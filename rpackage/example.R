### Code to simulate data and apply IPW and IPCW methods
### for adjusted difference in RMTL with competing risks

### Last updated by Sarah Conner on April 19 2021

library(survival)
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

# Indicators treating the competing event as censored
dat$i.event1 <- ifelse(dat$event==1, 1, 0)
dat$i.event2 <- ifelse(dat$event==2, 1, 0)



# --- Apply IPW estimator for marginal difference in RMTL ---

# Estimate IPW weights
ipw.mod <- glm(a ~ z1 + z2 + z3 + z4 + z5, data=dat, family='binomial')
p.a <- predict(ipw.mod, dat, type='response')
dat$weight <- (dat$a/p.a) + (1-dat$a)/(1-p.a)

# Estimate RMTLs and export results
ipw.res <- rmtl(times=dat$x, event=dat$event, eoi=1, group=as.factor(dat$a), weight=dat$weight, tau=1)

# Plot CIFs from RMTL output
plot(ipw.res$cif[[1]]$tj, ipw.res$cif[[1]]$cif, type='s', xlab='Time', ylab='Cumulative incidence')
lines(ipw.res$cif[[2]]$tj, ipw.res$cif[[2]]$cif, type='s', col=2)
abline(v=1, col=1, lty=3, lwd=2)

# Alternatively, use the CIF function
ipw.cif <- cif(times=dat$x, event=dat$event, eoi=1, group=as.factor(dat$a), weight=dat$weight)
plot(ipw.cif[[1]]$tj, ipw.cif[[1]]$cif, type='s', xlab='Time', ylab='Cumulative incidence', ylim=c(0,1))
lines(ipw.cif[[1]]$tj, ipw.cif[[1]]$cif.cil, type='s', lty=2)
lines(ipw.cif[[1]]$tj, ipw.cif[[1]]$cif.ciu, type='s', lty=2)
lines(ipw.cif[[2]]$tj, ipw.cif[[2]]$cif, type='s', col=2)
lines(ipw.cif[[2]]$tj, ipw.cif[[2]]$cif.cil, type='s', col=2, lty=2)
lines(ipw.cif[[2]]$tj, ipw.cif[[2]]$cif.ciu, type='s', col=2, lty=2)



# --- Apply IPCW regression model for difference in RMTL conditional on covariates ---

# Be sure to check dataset has non-missing values for exposure, covariates, event time, and event type
covariates <- dat[, c('a', 'z1', 'z2', 'z3', 'z4', 'z5')]

# Censoring weights estimated from whole sample (default option)
mod.res <- rmtl_mod(times=dat$x, event=dat$event, eoi=1, tau=1, cov=as.matrix(covariates), strata=FALSE)

# Can use strata=TRUE to derive censoring weights in each exposure arm
mod.res2 <- rmtl_mod(times=dat$x, event=dat$event, eoi=1, tau=1, cov=as.matrix(covariates), strata=TRUE, group=as.factor(dat$a))

# Compare IPW marginal difference in RMTL to conditional difference in RMTL (beta coefficient of a)
ipw.res$rmtl.diff
mod.res$res[2,]
mod.res2$res[2,]




# --- More figures: plot IPW CIF and 1-KM (without competing risks) --- #

# Plot with and without accounting for competing risks
# i.e. treat competing event as censored (1- Kaplan Meier method, Conner Stat Med 2019)
# Notice the dotted 1-KM curves overestimate the CIF for each event and in each exposure group

ipw.cif.naive <- cif(times=dat$x, event=dat$i.event1, eoi=1, group=as.factor(dat$a), weight=dat$weight)

plot(ipw.cif.naive[[1]]$tj, ipw.cif.naive[[1]]$cif, type='s', lty=2, xlab='Time', ylab='Cumulative incidence')
lines(ipw.cif.naive[[2]]$tj, ipw.cif.naive[[2]]$cif, type='s', col=2, lty=2)
lines(ipw.cif[[1]]$tj, ipw.cif[[1]]$cif, type='s')
lines(ipw.cif[[2]]$tj, ipw.cif[[2]]$cif, type='s', col=2)
legend('bottomright', legend=c('A=0, 1-KM', 'A=0, CIF', 'A=1, 1-KM', 'A=1, CIF'), col=c(1,1,2,2), lty=c(2,1,2,1))



# --- IPCW model ignoring competing risks ----

# Notice betas differ from model accounting for competing risk
# This result is equivalent to using rmst2() in package survRM2, but with flipped magnitude (difference in RMTL instead of RMTL)
mod.res_ignore <- rmtl.ipcw(times=dat$x, event=dat$i.event1, tau=1, cov=as.matrix(covariates), strata=FALSE)




