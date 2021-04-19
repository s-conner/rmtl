### Code to simulate data and apply IPW and IPCW methods
### for adjusted difference in RMTL with competing risks

### Last updated by Sarah Conner on October 22 2020

library(survival)
library(dichromat)
#source("functions//rmtl_ipw_function.R")
#source("functions//rmtl_ipcw_function.R")

library(rmtl)

# --- Function to simulate competing risks data ---

simdata <- function(n, alpha0, gamma, rho, psi11, psi21, cens){

  psi31 <- -0.4
  psi1 <- c(psi11, log(0.9), log(0.8), log(0.7), log(0.6), log(0.5))
  psi2 <- c(psi21, 0, 0, 0, 0, 0)
  psi3 <- c(psi31, log(0.5), log(0.6), log(0.7), log(0.8), log(0.9))

  # Covariates
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  z3 <- rnorm(n)
  z4 <- rnorm(n)
  z5 <- rnorm(n)
  z <- cbind(rep(1, n), z1, z2, z3, z4, z5)

  # Exposure
  alpha <- c(alpha0, log(0.7), log(0.6), log(0.5), log(0.9), log(0.8))
  lp.a <- z %*% alpha
  p.a <- 1/(1 + exp(-lp.a))
  a <- rbinom(n, 1, p.a)
  zz <- cbind(a, z1, z2, z3, z4, z5)

  # Event type
  p2 <- exp((gamma * exp(zz %*% psi1))/(rho + zz %*% psi2))
  cause <- 1 + rbinom(n, 1, p2)
  ev1 <- 2 - cause

  # Event times, obtained from conditional CIFs (on event type, covariatse, exposure)
  u <- runif(n, 0, 1)
  aa <- 1 - exp(gamma * exp(zz %*% psi1)/(rho + zz %*% psi2))
  bb <- (log(1 - u*aa) * (rho + zz %*% psi2))/(gamma * exp(zz %*% psi1))
  t <- (ev1) * (log(1 - bb))/(rho + zz %*% psi2) + (1-ev1)*-(log(1 - u))/exp(zz %*% psi3)

  # Add censoring
  c <- runif(n, 0, cens)
  i.cens <- ifelse(t<=c, 0, 1)
  event <- cause*(1-i.cens)
  x <- pmin(t, c)

  # Lo and behold
  dat <- cbind(zz, t, c, x, cause, i.cens, event)
  dat <- data.frame(dat)
  colnames(dat)[c(7, 9, 11, 12)] <- c('t', 'x', 'i.cens', 'event')
  return(dat)
}




# --- Generate data ---
dat <- simdata(n=500, gamma=1.5, rho=-1.5, psi11=0.5, psi21=-2, alpha0=log(.25/(1-.25)), cens=2.18)

dat$i.event1 <- ifelse(dat$event==1, 1, 0)
dat$i.event2 <- ifelse(dat$event==2, 1, 0)




# --- Apply IPW estimator for marginal difference in RMTL ---

# Estimate IPW weights
ipw.mod <- glm(a ~ z1 + z2 + z3 + z4 + z5, data=dat, family='binomial')
p.a <- predict(ipw.mod, dat, type='response')
dat$weight <- (dat$a/p.a) + (1-dat$a)/(1-p.a)

# Estimate RMTLs and export results
ipw_res <- rmtl(times=dat$x, event=dat$event, eoi=1, group=as.factor(dat$a), weight=dat$weight, tau=1)

# Plot CIFs from RMTL output
plot(ipw_res$cif[[1]]$tj, ipw_res$cif[[1]]$cif, type='s', xlab='Time', ylab='Cumulative incidence')
lines(ipw_res$cif[[2]]$tj, ipw_res$cif[[2]]$cif, type='s', col=2)
abline(v=1, col=1, lty=3, lwd=2)



# --- Apply IPCW regression model for difference in RMTL conditional on covariates ---

# Be sure to check dataset has non-missing values for exposure, covariates, event time, and event type
covariates <- dat[, c('a', 'z1', 'z2', 'z3', 'z4', 'z5')]

# Censoring weights estimated from whole sample (default option)
# Can use strata=TRUE to derive censoring weights in each exposure arm
ipcw_res <- rmtl_mod(times=dat$x, event=dat$event, eoi=1, tau=1, cov=as.matrix(covariates), strata=FALSE)


# Compare IPW marginal difference in RMTL to conditional difference in RMTL; both ~0.06
ipw_res$rmtl.diff
ipcw_res$res[2,]





# --- Plot IPW estimator --- #
# Plot with and without accounting for competing risks
# i.e. treat competing event as censored (1- Kaplan Meier method, Conner Stat Med 2019)
# Notice the dotted 1-KM curves overestimate the CIF for each event and in each exposure group


# Subset for plots
dat.a0 <- dat[dat$a==0,]
dat.a1 <- dat[dat$a==1,]


# Function for IPW CIF and 1-KM

cif <- function(times, event, weight=NULL, tau){

  if(is.null(weight)){weight <- rep(1, length(times))}
  entry=rep(0, length(times))

  data <-  data.frame(entry, times, event, weight)
  data$event[data$times>tau] <- 0
  data$times[data$times>tau] <- tau

  tj <- data$times[data$event!=0]
  tj <- unique(tj[order(tj)])
  num.tj <- length(tj)

  num.atrisk <- sapply(tj, function(tj) sum(data$weight[data$entry<tj & data$times>=tj]))
  num.ev1 <- sapply(tj, function(tj) sum(data$weight[data$event==1 & data$times==tj]))
  num.ev2 <- sapply(tj, function(tj) sum(data$weight[data$event==2 & data$times==tj]))
  num.ev <- num.ev1 + num.ev2

  h1 <- num.ev1/num.atrisk
  h2 <- num.ev2/num.atrisk
  h <- num.ev/num.atrisk

  s <- cumprod(c(1, 1 - h))
  s <- s[1:length(s)-1]

  theta1 <- s * h1
  theta2 <- s * h2

  cif1 <- cumsum(theta1)
  cif2 <- cumsum(theta2)

  cr1 <- 1 - cumprod(c(1, 1 - h1[1:num.tj-1]))
  cr2 <- 1 - cumprod(c(1, 1 - h2[1:num.tj-1]))

  return(list(tj=tj, cif1=cif1, cif2=cif2, cr1=cr1, cr2=cr2))
}

#  Derive CIFs and 1-KMs in each exposure group

a0 <- cif(dat.a0$x, dat.a0$event, dat.a0$weight, tau=1)
a1 <- cif(dat.a1$x, dat.a1$event, dat.a1$weight, tau=1)


# Plot

palette <- colorschemes$BrowntoBlue.10[c(2,9)]

par(mfrow=c(1,2))

plot(a0$tj, a0$cif1, type='s', col=palette[1], lty=1, ylim=c(0, 1), main='Event 1',
     xlab='Time (years)', ylab='Cumulative incidence', lwd=2)
lines(a0$tj, a0$cr1, type='s', col=palette[1], lty=2, lwd=2)
lines(a1$tj, a1$cif1, type='s', col=palette[2], lty=1, lwd=2)
lines(a1$tj, a1$cr1, type='s', col=palette[2], lty=2, lwd=2)
legend('topleft', legend=c('A=0, IPW CIF', 'A=0, IPW 1-KM', 'A=1, IPW CIF', 'A=1, IPW 1-KM'),
       lty=c(1,2,1,2), col=c(palette[1],palette[1],palette[2],palette[2]), lwd=rep(2,4))

plot(a0$tj, a0$cif2, type='s', col=palette[1], lty=1, ylim=c(0, 1), main='Event 2',
     xlab='Time (years)', ylab='Cumulative incidence', lwd=2)
lines(a0$tj, a0$cr2, type='s', col=palette[1], lty=2, lwd=2)
lines(a1$tj, a1$cif2, type='s', col=palette[2], lty=1, lwd=2)
lines(a1$tj, a1$cr2, type='s', col=palette[2], lty=2, lwd=2)
legend('topleft', legend=c('A=0, IPW CIF', 'A=0, IPW 1-KM', 'A=1, IPW CIF', 'A=1, IPW 1-KM'),
       lty=c(1,2,1,2), col=c(palette[1],palette[1],palette[2],palette[2]), lwd=rep(2,4))





# --- IPCW model ignoring competing risks ----
# Notice betas differ from model accounting for competing risk
# This result is equivalent to using rmst2() in package survRM2
ipcw_res_ignore <- rmtl.ipcw(times=dat$x, event=dat$i.event1, tau=1, cov=as.matrix(covariates), strata=FALSE)




