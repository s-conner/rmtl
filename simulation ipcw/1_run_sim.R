### Simulate competing risks data and assess IPCW regression method
### Sarah Conner, Updated October 22 2020

# Load IPCW function
source("//restricted//projectnb//conner-thesis//rmtl//rmtl_functions//rmtl_ipcw_function.R")

param <- read.csv("sim_reg_scenarios.csv", header=TRUE, sep=",")

nsim <- 6000
beta1 <- matrix(NA, nrow=nsim, ncol=6)
beta2 <- matrix(NA, nrow=nsim, ncol=6)

p <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
if (is.na(p)) p <- 1 

simtime <- function(n, cens, gamma, rho, z1, z2){
  
  # Event type
  p2 <- exp(gamma/rho)
  cause <- 1 + rbinom(n, 1, p2)
  ev1 <- 2 - cause
  
  # Event times, obtained from conditional CIFs (on event type, exposure)
  u <- runif(n, 0, 1)
  aa <- 1 - exp(gamma/rho)
  bb <- (log(1 - u*aa) * rho)/gamma
  t <- ev1 * (log(1 - bb))/rho + (1-ev1)*-(log(1 - u))
  
  # Add censoring
  c <- runif(n, 0, cens)
  i.cens <- ifelse(t<=c, 0, 1)
  event <- cause*(1-i.cens)
  x <- pmin(t, c)
  
  # Lo and behold
  dat <- cbind(z1, z2, t, c, x, cause, i.cens, event)
}

simdata <- function(n, p.z, cens, gamma00, rho00, gamma10, rho10, gamma01, rho01, gamma11, rho11){

  dat00 <- simtime(n*(1-p.z)^2, cens, gamma00, rho00, 0, 0)
  dat10 <- simtime(n*p.z*(1-p.z), cens, gamma10, rho10, 1, 0)
  dat01 <- simtime(n*p.z*(1-p.z), cens, gamma01, rho01, 0, 1)
  dat11 <- simtime(n*p.z^2, cens, gamma11, rho11, 1, 1)
  
  dat <- data.frame(rbind(dat00, dat01, dat10, dat11))
  colnames(dat) <- c('z1', 'z2', 't', 'c', 'x', 'cause', 'i.cens', 'event')
  return(dat)
}


for(i in 1:nsim){
  
  # Generate data
  dat <- simdata(n=param$n[p], p.z=param$p.z[p], cens=param$cens[p], 
                 gamma00=param$gamma00[p], rho00=param$rho00[p], gamma10=param$gamma10[p], rho10=param$rho10[p], 
                 gamma01=param$gamma01[p], rho01=param$rho01[p], gamma11=param$gamma11[p], rho11=param$rho11[p])
  
  # IPCW model, overall weights
  res.ipcw.i <- rmtl.ipcw(times=dat$x, event=dat$event, eoi=1, cov=as.matrix(dat[,c('z1', 'z2')]), tau=1, group=NULL, strata=FALSE)
  beta1[i,] <- c(i, res.ipcw.i$res[2, c(1,2,6,3,4)])
  beta2[i,] <- c(i, res.ipcw.i$res[3, c(1,2,6,3,4)])
}


# Export all iterations
colnames(beta1) <- c('simid', 'rmtldiff', 'rmtldiff.se', 'rmtldiff.p', 'rmtldiff.cil', 'rmtldiff.ciu')
colnames(beta2) <- c('simid', 'rmtldiff', 'rmtldiff.se', 'rmtldiff.p', 'rmtldiff.cil', 'rmtldiff.ciu')
write.csv(beta1, paste0("output//beta1_", p, ".csv"), row.names=FALSE)
write.csv(beta2, paste0("output//beta2_", p, ".csv"), row.names=FALSE)



