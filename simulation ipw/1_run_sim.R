### Code to replicate simulation study to assess IPW estimator 
### of adjusted difference in RMTL with competing risks

### Last updated by Sarah Conner on October 22 2020


# --- Load parameters for each scenario ---

# Each row contains parameters for data generation
param <- read.csv("sim_scenarios.csv", header=TRUE, sep=",")

nsim <- 6000
simresults <- matrix(NA, nrow=nsim, ncol=10)

# If using parallel jobs in Linux environment, can run each scenario in parallel
# Otherwise, just set p to value 1-16 for desired scenario
p <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
if (is.na(p)) p <- 1 


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


# --- Simplified RMTL IPW function for reduced computation time and extracting results ---

rmtl.sim <- function(entry=NULL, times, event, eoi=1, group=NULL, weight=NULL, tau=NULL, alpha=.05){  
  
  if(sum(times<0)>0){print("Error: times must be positive.")
  }else{
    if(sum(weight<=1)>0){print("Error: weights must be greater than 1.")
    }else{
      
      #--- Prep input data ---
      if(is.null(entry)){entry <- rep(0, length(times))}
      if(is.null(weight)){weight <- rep(1, length(times))}
      alldat <- data.frame(entry, times, event, group, weight)
      alldat <- alldat[!is.na(alldat$group) & !is.na(alldat$times),]
      alldat <- alldat[order(group),] 
      
      
      #--- If tau not specified, use minimum time from all groups. Check if provided tau is appropriate. ---
      gg <- length(levels(alldat$group))
      
      if(is.null(tau)){
        
        taui <- rep(NA, gg)
        
        for (i in (1:gg)){
          groupval <- (levels(alldat$group)[i])
          dat_group <- alldat[which(alldat$group==(groupval)),]
          taui[i] <- max(dat_group$times)
        }
        tau <- min(taui)
        
      } else {
        
        tau.error <- rep(0, gg)
        
        for (i in (1:gg)){
          groupval <- (levels(alldat$group)[i])
          dat_group <- alldat[which(alldat$group==(groupval)),]
          tau.error[i] <- ifelse(max(dat_group$times)<tau, 1, 0)
        }
      }
      
      if(sum(tau.error)>0){
        print("Error: observed times do not reach tau in each exposure group. Choose a different tau or leave unspecified for default value of tau.")
      }else{
        
        
        #--- Proceed with RMTL ----
        
        alldat$event[alldat$times>tau] <- 0
        alldat$times[alldat$times>tau] <- tau
        
        rmtl <- rep(NA, length(1:gg))
        groupval <- rep(NA, length(1:gg))
        rmtl.se <- rep(NA, length(1:gg))
        
        for (g in 1:gg){
          
          #--- Derive CIF and related quantities (theta, hazards, etc) ---
          
          groupval[g] <- (levels(alldat$group)[g])
          data <- alldat[which(alldat$group==(groupval[g])),]
          
          tj <- data$times[data$event!=0]
          tj <- unique(tj[order(tj)])
          num.tj <- length(tj)
          
          num.atrisk <- sapply(tj, function(x) sum(data$weight[data$entry<x & data$times>=x]))
          num.ev1 <- sapply(tj, function(x) sum(data$weight[data$event==eoi & data$times==x]))
          num.ev2 <- sapply(tj, function(x) sum(data$weight[data$event!=eoi & data$event!=0 & data$times==x]))
          num.ev <- num.ev1 + num.ev2
          
          m <- sapply(tj, function(x){sum((data$weight[data$entry<x & data$times>=x])^2)})
          mg <- ((num.atrisk^2)/m)
          
          h1 <- num.ev1/num.atrisk
          h <- num.ev/num.atrisk
          
          s <- cumprod(c(1, 1 - h))
          s <- s[1:length(s)-1]
          
          theta <- s * h1
          cif1 <- cumsum(theta)
          
          
          #---  Variance of each theta --- 
          
          a <- c(0,cumsum(num.ev/(mg * (num.atrisk - num.ev))))
          a <- a[1:num.tj]
          
          var.theta <- ((theta)^2) * (((num.atrisk - num.ev1)/(mg * num.ev1)) + a)
          var.theta[is.nan(var.theta)] <- 0
          
          
          #---  Covariance of thetas --- 
          
          cov.theta <- matrix(NA, nrow=num.tj, ncol=num.tj)
          b <- c(0,cumsum(num.ev/(mg * (num.atrisk - num.ev))))
          
          for(j in 1:(num.tj-1)){
            for(k in (j+1):num.tj){
              cov.theta[k,j] <- cov.theta[j,k] <- (theta[j]) * (theta[k]) * (-1/mg[j] + b[j])
            }
          }
          
          # Diagonal is variance of thetas
          diag(cov.theta) <- var.theta
          
          
          #---  Covariances of CIF --- 
          
          cov.f10 <- apply(cov.theta, 2, function(x){x[is.na(x)] <- 0;  cumsum(x)})
          cov.f1 <- apply(cov.f10, 1, function(x){x[is.na(x)] <- 0;  cumsum(x)})
          
          var.f1 <- diag(cov.f1) # not sure if this is needed, but for sanity check
          
          
          #--- RMTL and variance ---
          
          areas <- c(tj[2:num.tj], tau)-tj
          rmtl[g] <- sum(areas*cif1)
          
          cov.weights <- outer(areas,areas)
          cov.f1.weight <- cov.weights * cov.f1
          
          rmtl.var <- sum(cov.f1.weight)
          rmtl.se[g] <- sqrt(rmtl.var) 
          
        }
        
        
        #--- Compare RMTL between groups and compile output---
        
        z <- qnorm(1-alpha/2)
        rmtldiff <- rmtl[2]-rmtl[1]
        rmtldiff.se <- sqrt(rmtl.se[2]^2 + rmtl.se[1]^2)
        
        results <- c(rmtl0=rmtl[1], rmtl0.se=rmtl.se[1], rmtl1=rmtl[2], rmtl1.se=rmtl.se[2],
                     rmtldiff=rmtldiff, rmtldiff.se=rmtldiff.se, rmtldiff.p=2*(1-pnorm(abs(rmtldiff/rmtldiff.se))),
                     rmtl.cil=rmtldiff-(z*rmtldiff.se), rmtl.ciu=rmtldiff+(z*rmtldiff.se))
        
        return(results)
        
      }  
    }
  }
}


# --- For each iteration, generate data, apply IPW RMTL method, and save result ---

for(i in 1:nsim){
  
  # Generate data
  dat <- simdata(n=param$n[p], gamma=param$gamma[p], rho=param$rho[p], psi11=param$psi11[p], psi21=param$psi21[p],
                 alpha0=param$alpha0[p], cens=param$cens[p])
  
  # Estimate IPW weights
  ipw.mod <- glm(a ~ z1 + z2 + z3 + z4 + z5, data=dat, family='binomial')
  p.a <- predict(ipw.mod, dat, type='response')
  dat$weight <- (dat$a/p.a) + (1-dat$a)/(1-p.a)
  
  # Estimate RMTLs and export results
  res.i <- rmtl.sim(times=dat$x, event=dat$event, eoi=1, group=as.factor(dat$a), weight=dat$weight, tau=1)  
  simresults[i,] <- c(i, res.i)
}


# --- Export all iterations ---

colnames(simresults) <- c('simid', 'rmtl0', 'rmtl0.se', 'rmtl1', 'rmtl1.se', 'rmtldiff', 'rmtldiff.se', 'rmtldiff.p', 'rmtldiff.cil', 'rmtldiff.ciu')
write.csv(simresults, paste0("simresults_", p, ".csv"), row.names=FALSE)

