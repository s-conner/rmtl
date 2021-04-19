### Function to estimate the unadjusted or IPW-adjusted CIF, allowing delayed entry
### Null weights yields the unadjusted CIF
### Sarah Conner, Updated April 19 2021



#' Cumulative incidence function (CIF) with inverse probability weighting (IPW)
#'
#' @param times time to event (any-cause) or censoring
#' @param event integer denoting the event type, where 0=censored. If there is only 1 event type (values take 0,1), results are equivalent to a non-competing risk setting
#' @param eoi the event of interest, which should be an event value (defaults to 1)
#' @param group exposure/treatment group. If left blank, function calculates the CIF in the overall sample.
#' @param weight user can specify weights to adjust for confounding, i.e. inverse probability weights, overlap weights, or survey weights. If left blank, function calculates unadjusted results.
#' @param alpha level of significance (defaults to 0.05).
#' @param entry optional argument to allow delayed entry times (i.e., age as the time scale)
#'
#' @return IPW (or unadjusted) CIF in each group, confidence intervals, and timepoints for plotting cumulative incidence curves.
#' @export
#'
#' @examples
cif <- function(times, event, eoi=1, group=NULL, weight=NULL, alpha=.05, entry=NULL){

  #--- Prep input data ---
  if(is.null(entry)){entry <- rep(0, length(times))}
  if(is.null(group)){group <- as.factor(rep(1, length(times)))}
  if(is.null(weight)){weight <- rep(1, length(times))}
  alldat <- data.frame(entry, times, event, group, weight)
  alldat <- alldat[!is.na(alldat$group) & !is.na(alldat$times),]
  alldat <- alldat[order(group),]

  z <- stats::qnorm(1-alpha/2)
  gg <- length(levels(alldat$group))
  res.cif <- list()

  for (g in 1:gg){

    #--- Derive CIF and related quantities (theta, hazards, etc) ---

    groupval <- (levels(alldat$group)[g])
    data <- alldat[which(alldat$group==(groupval)),]

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

        # Modified formula
        cov.theta[k,j] <- cov.theta[j,k] <- (theta[j]) * (theta[k]) * (-1/mg[j] + b[j])
      }
    }

    # Diagonal is variance of thetas
    diag(cov.theta) <- var.theta


    #---  Covariances of CIF ---

    cov.f10 <- apply(cov.theta, 2, function(x){x[is.na(x)] <- 0;  cumsum(x)})
    cov.f1 <- apply(cov.f10, 1, function(x){x[is.na(x)] <- 0;  cumsum(x)})
    var.f1 <- diag(cov.f1)


    #---  Export CIF and variance ---

    res.cif.g[[length(res.cif.g)+1]] <- cif1
    res.cif.g[[length(res.cif.g)+1]] <- sqrt(var.f1)
    res.cif.g[[length(res.cif.g)+1]] <- tj
    res.cif.g[[length(res.cif.g)+1]] <- cif1 - z*sqrt(var.f1)
    res.cif.g[[length(res.cif.g)+1]] <- cif1 + z*sqrt(var.f1)
    names(res.cif.g) <- c('cif', 'se.cif', 'tj', 'cif.cil', 'cif.ciu')
    res.cif[[length(res.cif)+1]] <- res.cif.g
  }

  names(res.cif) <- groupval
  return(res.cif)

}



