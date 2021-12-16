#' Eliminate Multicollinearity in Multiple Regression
#'
#' This function loads data in two pieces, X and y. Then it selects the subset
#' of variables that minimize multicollinearity using a mixed-integer optimization
#' framework.
#'
#' @param X Dataframe containing potential covariates and their values
#' @param y Dataframe with response variables of interest
#' @param alpha Optional- standard constraint > 1
#' @return A list containing indicators for the selected variables ("Selected") and the
#' corresponding total number of covariates selected ("S"), R^2 for the linear multiple
#' regression model ("Rsq"), and maximum remaining VIF
#' value of the selected predictors ("VIF_max")
#'
#' @import Matrix gurobi
#' @export
vifmiqo <- function(X, y, alpha=5){

  # Reformatting data
  X <- apply(X,2,function(a) (a-mean(a))/sqrt(sum((a-mean(a))^2)))
  y <-(y-mean(y))/sqrt(sum((y-mean(y))^2))

  p=ncol(X)
  R = t(X)%*%X
  M=1e4

  #when constraint does not consider a or z
  zero_p <- Matrix(0,nrow=p,ncol=p,sparse = TRUE)
  #when constraint does not consider Q
  zero_p_q <- Matrix(matrix(rep(zero_p,p),nrow=p),nrow=p,sparse = TRUE)
  #when constraint does not consider U
  zero_p_u <- Matrix(matrix(rep(zero_p,p),nrow=p),nrow=p,sparse = TRUE)
  #pxp Identity matrix
  eye = diag(p)

  #================FORMULATE CONSTRAINTS ONE-AT-A-TIME============================

  #27
  #a_j - M*z_j <= 0
  L = matrix(cbind(eye,zero_p_q,zero_p_u,(eye*-M)),nrow=p)
  L = Matrix(L,sparse=TRUE)
  dir = rep('<=',p)
  sense = rep('<',p)
  rhs = rep(0,p)

  #a_j + M*z_j >= 0
  L = rbind(L,matrix(cbind(eye,zero_p_q,zero_p_u,(eye*M)),nrow=p))
  dir = c(dir,rep('>=',p))
  sense = c(sense,rep('>',p))
  rhs = c(rhs,rep(0,p))

  #28
  #q_ll <= alpha
  Q_28 <- matrix(0,nrow=p,ncol=p*p)
  for (i in 1:p){
    Q_28[i,(p*(i-1))+i] <- 1
  }

  L = rbind(L,matrix(cbind(zero_p,Q_28,zero_p_u,zero_p),nrow=p))
  dir = c(dir,rep('<=',p))
  sense = c(sense,rep('<',p))
  rhs = c(rhs,rep(alpha,p))

  #29
  #Q*R_p + U  = I_p
  a_29= matrix(0,nrow=p*p,ncol=p)
  Q_29 = matrix(0,nrow=p*p,ncol=p*p)
  U_29 = matrix(0,nrow=p*p,ncol=p*p)
  for (i in 1:p) {
    Q_29[((p*(i-1))+1):(p*i),((p*(i-1))+1):(p*i)] <- t(R)
    U_29[((p*(i-1))+1):(p*i),((p*(i-1))+1):(p*i)] <- eye
  }
  z_29= matrix(0,nrow=p*p,ncol=p)

  L29 = Matrix(cbind(a_29,Q_29,U_29,z_29),nrow=p*p,ncol=ncol(L),sparse = TRUE)
  L = rbind(L,L29)

  dir = c(dir,rep('==',p*p))
  sense = c(sense,rep('=',p*p))
  rhs = c(rhs,as.vector(t(diag(p))))

  #31
  U_31a <- Matrix(0,nrow=p*p,ncol=p*p,sparse = TRUE)
  Z_31a <- Matrix(0,nrow=p*p,ncol=p,sparse = TRUE)
  for (j in 1:p) {
    for (i in 1:p){
      U_31a[((p*(j-1))+i),(p*(i-1))+j] <- 1
      Z_31a[((p*(j-1))+i),j] <- M
    }
  }
  A_31a <- Matrix(0,nrow=p*p,ncol=p,sparse = TRUE)
  Q_31a <- Matrix(0,nrow=p*p,ncol=p*p,sparse = TRUE)
  L31a = Matrix(cbind(A_31a,Q_31a,U_31a,Z_31a),nrow=p*p,ncol=ncol(L),sparse=TRUE)

  dir = c(dir,rep('<=',p*p))
  sense = c(sense,rep('<',p*p))
  rhs = c(rhs,rep(M,p*p))

  L31b = Matrix(cbind(A_31a,Q_31a,U_31a,(Z_31a*-1)),nrow=p*p,ncol=ncol(L),sparse=TRUE)
  L31 <- rbind(L31a,L31b)
  L  = rbind(L,L31)

  dir = c(dir,rep('>=',p*p))
  sense = c(sense,rep('>',p*p))
  rhs = c(rhs,rep(-M,p*p))

  #32
  Q_32a <- Matrix(0,nrow=p*p,ncol=p*p,sparse = TRUE)
  Z_32a <- Matrix(0,nrow=p*p,ncol=p,sparse = TRUE)
  for (j in 1:p) {
    for (i in 1:p){
      Q_32a[((p*(j-1))+i),(p*(i-1))+j] <- 1
      Z_32a[((p*(j-1))+i),j] <- -M
    }
  }
  A_32a <- Matrix(0,nrow=p*p,ncol=p,sparse = TRUE)
  U_32a <- Matrix(0,nrow=p*p,ncol=p*p,sparse = TRUE)
  L32a = Matrix(cbind(A_32a,Q_32a,U_32a,Z_32a),nrow=p*p,sparse=TRUE)

  dir = c(dir,rep('<=',p*p))
  sense = c(sense,rep('<',p*p))
  rhs = c(rhs,rep(0,p*p))

  L32b = Matrix(cbind(A_32a,Q_32a,U_32a,(Z_32a*-1)),nrow=p*p,sparse=TRUE)
  L32a <- Matrix(L32a,sparse = TRUE)
  L32b <- Matrix(L32b,sparse = TRUE)
  L32 <- rbind(L32a,L32b)
  L  = rbind(L,L32)

  dir = c(dir,rep('>=',p*p))
  sense = c(sense,rep('>',p*p))
  rhs = c(rhs,rep(0,p*p))

  # Assigning values to the objective function
  Q_s <- R
  Q = Matrix(0,nrow=ncol(L),ncol=ncol(L),sparse = TRUE)
  Q[(1:p),(1:p)] <- Q_s
  L_s <- (-2*t(X)%*%y)
  C = rep(0,ncol(L))
  C[1:p] <- L_s

  # Model building
  model<- list()
  model$Q = Q
  model$lb <- rep(-Inf,ncol(L))
  model$A          <- L
  model$obj        <- C
  model$modelsense <- 'min'
  model$rhs        <- rhs
  model$sense      <- sense
  model$objcon      <- t(y)%*%y
  model$vtype <- c(rep("C",(ncol(L)-p)),rep("B",p))

  # Solving for solution
  res_gur <- gurobi(model)

  # Formatting the results of the model
  a <- res_gur$x[1:p]
  Q_mat <- matrix(res_gur$x[(p+1):((p*p)+(p))],ncol=p,nrow=p)
  U <- matrix(res_gur$x[((p*p)+(p+1)):((length(res_gur$x)-p))],ncol=p,nrow=p)
  # U = matrix(res_gur$x[381:741],nrow=p,ncol=p)
  z <- res_gur$x[((length(res_gur$x)-p)+1):length(res_gur$x)]

  num_selected = sum(z) # |S|
  vif_max = max(diag(Q_mat)) # VIF_{max}

  # To obtain R^2
  linear_trial <- X%*%a
  final_model = lm(y~linear_trial)
  Rsq = summary(final_model)$r.squared

  return(list(Selected=z, Rsq = Rsq, S = num_selected, VIF_max = vif_max))

}


