
### Function to estimate the unadjusted or IPW-adjusted RMTL
### Null weights yields the unadjusted RMTL
### Sarah Conner, package created March 12 2021


#' Restricted mean time lost (RMTL) with competing risks
#'
#' @param times time to event (any-cause) or censoring
#' @param event integer denoting the event type, where 0=censored. If there is only 1 event type (values take 0,1), results are equivalent to a non-competing risk setting
#' @param eoi the event of interest, which should be an event value. Defaults to 1.
#' @param group exposure/treatment group. If left blank, function calculates the RMTL in the overall sample.
#' @param weight user can specify weights to adjust for confounding, i.e. inverse probability weights, overlap weights, or survey weights. If left blank, function calculates unadjusted results.
#' @param tau time horizon used to restrict event times. User can supply a clinically meaningful time horizon, otherwise the function calculates the minimum of the maximum event times in each exposure/treatment group
#' @param alpha level of significance (defaults to 0.05).
#' @param entry optional argument to allow delayed entry times (i.e., age as the time scale)
#' @param ...
#'
#' @return IPW (or unadjusted) RMTL in each group, differences between groups, confidence intervals, p-values, and IPW cumulative incidence curves.
#' @export
#'
#' @examples
rmtl <- function(times, event, eoi=1, group=NULL, weight=NULL, tau=NULL, alpha=.05, entry=NULL, ...){

  if(sum(times<0)>0){print("Error: times must be positive.")
  }else{

    #--- Prep input data ---
    if(is.null(entry)){entry <- rep(0, length(times))}
    if(is.null(group)){group <- as.factor(rep(1, length(times)))}
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
      graphics::plot(NULL, xlim=c(0, tau), ylim=c(0,1), xlab='Time', ylab='Cumulative incidence', ...)

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
        graphics::lines(c(tj, tau), c(cif1, cif1[num.tj]), type="s", col=(g+2), lwd=2)


        #---  Variance of each theta ---

        a <- c(0,cumsum(num.ev/(mg * (num.atrisk - num.ev))))
        a <- a[1:num.tj]

        var.theta <- ((theta)^2) * (((num.atrisk - num.ev1)/(mg * num.ev1)) + a)
        var.theta[is.nan(var.theta)] <- 0

        #sum.var.theta <- cumsum(var.theta)


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



      #--- Add legend and tau to plot ---
      graphics::abline(v=tau, col=1, lty=3, lwd=2)
      if(gg>1){
        graphics::legend('topleft', groupval, lty=rep(1, gg), lwd=rep(2, gg), col=3:(gg+2), cex=.75, bty ="n")
      }


      #--- Compare RMTL between groups and compile output---

      z <- stats::qnorm(1-alpha/2)
      res <- data.frame(groupval, rmtl, rmtl.se, cil=rmtl-(z*rmtl.se), ciu=rmtl+(z*rmtl.se))

      pwc <- ((gg^2)-gg)/2   #number of pairwise comparisons


      #--- Calculate difference in RMTL if group is not specified

      if(pwc>0){
        label.diff <- rep(NA,pwc)
        rmtl.diff <- rep(NA,pwc)
        rmtl.diff.se <- rep(NA,pwc)
        rmtl.diff.cil <- rep(NA,pwc)
        rmtl.diff.ciu <- rep(NA,pwc)
        rmtl.diff.p <- rep(NA,pwc)

        res.diff <- data.frame(label.diff, rmtl.diff, rmtl.diff.se, rmtl.diff.cil, rmtl.diff.ciu, rmtl.diff.p)
        l <- 1

        for (i in 1:(gg-1)){
          for (ii in (i+1):gg){

            #--- RMTL Difference ---
            res.diff[l,]$label.diff <- paste(res[ii,]$groupval, '-', res[i,]$groupval)
            res.diff[l,]$rmtl.diff <- (res[ii,]$rmtl - res[i,]$rmtl)
            res.diff[l,]$rmtl.diff.se <- sqrt(res[ii,]$rmtl.se^2 + res[i,]$rmtl.se^2)

            res.diff[l,]$rmtl.diff.cil <- res.diff[l,]$rmtl.diff - z*res.diff[l,]$rmtl.diff.se
            res.diff[l,]$rmtl.diff.ciu <- res.diff[l,]$rmtl.diff + z*res.diff[l,]$rmtl.diff.se
            res.diff[l,]$rmtl.diff.p <- 2*(1-stats::pnorm(abs(res.diff[l,]$rmtl.diff)/res.diff[l,]$rmtl.diff.se))

            l <- l+1
          }
        }

      }


      if(pwc>0){

        allres <- list(rmtl=res, rmtl.diff=res.diff)

        cat("RMTL per group, tau=", tau, "\n\n", sep="")
        rmtl.round <- round(allres$rmtl[,c(2:5)],3)
        colnames(rmtl.round) <- c("RMTL", "SE", "CIL", "CIU")
        rownames(rmtl.round) <- c(levels(res[,1]))
        print(rmtl.round)
        cat("\n\n")

        cat ("RMTL Differences, tau=", tau, "\n\n", sep="")
        rmtl.diff.round <- round(allres$rmtl.diff[c(2:6)],3)
        colnames(rmtl.diff.round) <- c("RMTL Diff.", "SE", "CIL", "CIU", "p")
        rownames(rmtl.diff.round) <- c(allres$rmtl.diff[,1])
        print(rmtl.diff.round)
        cat("\n\n")

        invisible(allres)

      } else { # No groups, thus no pairwise comparison

        cat("RMTL per group, tau=", tau, "\n\n", sep="")
        colnames(res) <- c("Group", "RMTL", "SE", "CIL", "CIU")
        rownames(res) <- c(paste("Group", res$Group,' '))
        print(round(res[c(2:5)],3))
        cat("\n\n")

        invisible(res)

      }
    }
  }
}
