# Rvifmiqo
The 'Rvivmiqo' package and 'vifmiqo' funciton eliminate multicollinearity in multiple regression using a mixed-integer optimization framework. The 'vifmiqo' function returns a list of covariates to minimize multicollinearity and associated multiple regression attributes such as R-squared, maximum VIF value, and number of covariates.

# Tutorial
A detailed tutorial on how to use the 'vifmiqo' function along with a test example can be found in vignettes/Tutorial.Rmd.

# Installation
Prior to installation, visit https://www.gurobi.com/ to download the necessary 'gurobi' package.

```{r}
devtoools::install_git("egfried/Rvifmiqo")
```
