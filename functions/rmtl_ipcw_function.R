### Function for ICPW regression model with competing risks data
### Last updated by Sarah Conner on October 22, 2020

### Note: This code is modified from the original 'rmst2reg function' of the {survRM2} package,
### which was authored by Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, and Miki Horiguchi,
### in order to account for competing risks.

library(survival)

rmtl.ipcw <- function(times, event, eoi=1, tau, cov=NULL, strata=FALSE, group=NULL){
  
  if(is.null(group) & strata==TRUE){stop('Please specify a factor variable to statify weights.')}
  if(is.null(cov)){print('Warning: Fitting intercept-only model.')}
  
  # Round event times to avoid issues with survival() package rounding differently
  y <- round(times,4)
  id <- 1:length(y)
  
  # Recode so delta1 reflects event of interest, delta2 reflects all competing events. Assumes 0=censoring.
  delta1 <- ifelse(event==eoi, 1, 0)
  delta2 <- ifelse(event!=0 & event!=eoi, 1, 0)
  
  # Overall quantities
  x <- cbind(int=rep(1, length(y)), cov)
  p <- length(x[1,])
  if(is.null(group)){group <- as.factor(rep(1, length(y)))}
  
  # Recode event indicators to reflect status at chosen tau
  delta1[y>tau] <- 0
  delta2[y>tau] <- 0
  
  y <- pmin(y, tau)
  y1 <- y*delta1
  
  d0 <- 1 - (delta1 + delta2) # censoring indicator
  d0[y==tau] <- 0  # If follow-up lasts til tau, the event will not count as 'censored' in IPCW weights
  weights <- NULL
  
  ## Calculate IPCW weights (option to stratify by group) ## 
  
  if(strata==TRUE){
    for(aa in 1:length(unique(group))){
      # Subset the group
      a <- unique(group)[aa]
      d0.a <- d0[group==a]
      delta1.a <- delta1[group==a]
      y.a <- y[group==a]
      x.a <- x[group==a,]
      n.a <- length(d0.a)
      orig.id.a0 <- orig.id.a <- id[group==a]
      
      # Order the event times
      id.a <- order(y.a)
      y.a <- y.a[id.a]
      d0.a <- d0.a[id.a]
      delta1.a <- delta1.a[id.a]
      x.a <- x.a[id.a,]
      orig.id.a <- orig.id.a[id.a]
      
      # Derive IPCW
      fit <- survfit(Surv(y.a, d0.a) ~ 1)
      weights.a <- (1-d0.a)/rep(fit$surv, table(y.a))
      
      # Need to assign weights accordig to original ID, not ordered by event time
      linked.weights.a <- cbind(orig.id.a, weights.a, delta1.a, d0.a, y.a)
      weights <- rbind(weights, linked.weights.a)
    }
  } else {
    
    # Order the event times
    id.a <- order(y)
    y.a <- y[id.a]
    d0.a <- d0[id.a]
    delta1.a <- delta1[id.a]
    x.a <- x[id.a,]
    orig.id.a <- id[id.a]
    
    # Derive IPCW
    fit <- survfit(Surv(y.a, d0.a) ~ 1)
    weights.a <- (1-d0.a)/rep(fit$surv, table(y.a))
    
    # Need to assign weights accordig to original ID, not ordered by event time
    linked.weights.a <- cbind(orig.id.a, weights.a, delta1.a, d0.a, y.a)
    weights <- rbind(weights, linked.weights.a)
  }
  
  
  ## Fit linear model ## 
  
  # Link weights to original data frame
  #colnames(weights) <- c('id', 'weights')
  #data <- merge(data0, weights, by='id')
  #summary(lm(tau-y ~ x-1, weights=weights, data=data))
  
  # Or, sort weights and use vectors
  w <- weights[order(weights[, 1]),2]
  lm.fit <- lm(delta1*(tau-y) ~ x-1, weights=w)
  
  
  ## Derive SE ##
  
  beta0 <- lm.fit$coef
  error <- tau - y - as.vector(x %*% beta0)
  score <- x * w * error
  
  # Kappa (sandwich variance components) stratified by group
  kappa <- NULL
  
  for(aa in 1:length(unique(group))){
    
    # Subset the group
    a <- unique(group)[aa]
    d0.a <- d0[group==a]
    delta1.a <- delta1[group==a]
    y.a <- y[group==a]
    x.a <- x[group==a,]
    n.a <- length(d0.a)
    orig.id.a0 <- orig.id.a <- id[group==a]
    score.a <- score[group==a,]
    
    # Kappa calculations for sandwich variance
    kappa.a <- matrix(0, n.a, p)
    
    for(i in 1:n.a){
      kappa1 <- score.a[i,]
    
      kappa2 <- apply(score.a[y.a>=y.a[i],,drop=F], 2, sum)*(d0.a[i])/sum(y.a>=y.a[i])
    
      kappa3 <- rep(0, p)
    
      for(k in 1:n.a){
        if(y.a[k]<=y.a[i]){
          kappa3 <- kappa3+apply(score.a[y.a>=y.a[k],,drop=F], 2, sum)*(d0.a[k])/(sum(y.a>=y.a[k]))^2
        }
      }
  
      kappa.a[i,] <- kappa1+kappa2-kappa3
    }
    kappa <- rbind(kappa, kappa.a)
  }
  
  # Transpose the kappas rbinded from each group gives pxp matrix
  gamma <- t(kappa) %*% kappa
  
  A <- t(x) %*% x
  varbeta <- solve(A) %*% gamma %*% solve(A)
  se <- sqrt(diag(varbeta))
  
  
  #--- Return results ---
  
  res <- cbind(beta=lm.fit$coef, se=se, cil=lm.fit$coef-(1.96*se), ciu=lm.fit$coef+(1.96*se), 
               z=lm.fit$coef/se, p=2*(1-pnorm(abs(lm.fit$coef/se))))
  rownames(res) <- c("Intercept", colnames(x[,-1]))
  
  allres <- list(res=res, varbeta=varbeta)
  print(round(res, 3))
  invisible(allres)
}  
  