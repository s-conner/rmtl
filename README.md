# Estimators of the restricted mean time lost (RMTL) in the presence of competing risks
Conner and Trinquart (2021) Statistics in Medicine

The folders here include R functions for methods introduced in Conner and Trinquart (2021) Statistics in Medicine: 
1. estimation of the unadjusted RMTL
2. inverse probability weighted (IPW)-adjusted RMTL
3. RMTL regression model estimated with inverse probabilty of censoring weighting (IPCW). 

Folders 'simulation ipcw' and 'simulation ipw' include codes to replicate our simulation studies.

I created an R package, which you can install with the code below!

library(devtools)
install_github("s-conner/rmtl/rmtl package")
