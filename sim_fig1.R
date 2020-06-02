### Code to replicate Fig 1 (true CIFs in simulation study)
### Last updated by Sarah Conner on June 2 2020

library(dichromat)

# n.iter should be 1,000,000 to replicate the figures exactly, but a smaller number may be used to reduce computing time
n.iter <- 1000


#pdf(width=11, height=6, file = "cifs_sidebyside.pdf")

par(mfrow=c(1,2))


# Modified RMTL function (does not compute variance) to reduce computation time.

rmtl.plot <- function(entry=NULL, times, event, eoi=1, group=NULL, weight=NULL, tau=NULL, alpha=.05, yaxismax=1, title=NULL){  
  
  if(sum(times<0)>0){print("Error: times must be positive.")
  }else{
    if(sum(weight<=1)>0){print("Error: weights must be greater than 1.")
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
        
        
        palette <- colorschemes$BrowntoBlue.10[c(2,9)]
        plot(NULL, xlim=c(0, tau), ylim=c(0,yaxismax), xlab='Time (years)', ylab='Cumulative incidence',
             main=title)
        abline(v=tau, col='grey', lty=3, lwd=2)
        
        
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
          
          
          #--- Plot CIF ---
          lines(c(tj, tau), c(cif1, cif1[num.tj]), type="s", col=palette[g], lty=g, lwd=2)
          
          
          #--- RMTL ---
          areas <- c(tj[2:num.tj], tau)-tj
          rmtl[g] <- sum(areas*cif1)
        }
        
        
        #--- Add legend to plot ---
        legend('topleft', paste("A =", groupval), col=palette, lty=1:gg, lwd=rep(2, gg), 
               cex=.75, bty ="n", inset = c(0, 0))
        
        
        #--- Compare RMTL between groups and compile output---
        
        results <- data.frame(groupval, rmtl)
        pwc <- ((gg^2)-gg)/2   #number of pairwise comparisons
        
        label.diff <- rep(NA,(pwc))
        rmtl.diff <- rep(NA,(pwc))
        res.diff <- data.frame(label.diff, rmtl.diff)
        
        l <- 1
        
        for (i in 1:(gg-1)){
          for (ii in (i+1):gg){
            
            #--- RMTL Difference ---
            res.diff[l,]$label.diff <- paste('Groups',results[ii,]$groupval,'vs.',results[i,]$groupval,' ')
            res.diff[l,]$rmtl.diff <- (results[ii,]$rmtl - results[i,]$rmtl)
            l <- l+1
          }
        }
        
        cat ("RMTL per group \n\n")
        colnames(results) <- c("Group", "RMTL")
        rownames(results) <- c(paste("Group", results$Group,' '))
        print(round(results[2],3))
        cat("\n\n")
        
        cat ("RMTL Differences \n\n")
        colnames(res.diff) <- c("Groups", "RMTL Diff.")
        rownames(res.diff) <- c(res.diff$Groups)
        print(round(res.diff[2],3))
        cat("\n\n")
        
        
      }  
    }
  }
}


cif.plot <- function(n, gamma, rho, psi11, psi21, title){
  
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
  
  
  # ---- Exposed ---- 
  a.1 <- rep(1, n)
  zz.1 <- cbind(a.1, z1, z2, z3, z4, z5)
  
  # Event type
  p2.1 <- exp((gamma * exp(zz.1 %*% psi1))/(rho + zz.1 %*% psi2))
  cause.1 <- 1 + rbinom(n, 1, p2.1)
  ev1.1 <- 2 - cause.1
  
  # Event times, obtained from conditional CIFs (on event type, covariates, exposure)
  u.1 <- runif(n, 0, 1)
  aa.1 <- 1 - exp(gamma * exp(zz.1 %*% psi1)/(rho + zz.1 %*% psi2))
  bb.1 <- (log(1 - u.1*aa.1) * (rho + zz.1 %*% psi2))/(gamma * exp(zz.1 %*% psi1))
  t.1 <- (ev1.1) * (log(1 - bb.1))/(rho + zz.1 %*% psi2) + (1-ev1.1)*-(log(1 - u.1))/exp(zz.1 %*% psi3)
  
  
  # ---- Unexposed ---- 
  a.0 <- rep(0, n)
  zz.0 <- cbind(a.0, z1, z2, z3, z4, z5)
  
  # Event type
  p2.0 <- exp((gamma * exp(zz.0 %*% psi1))/(rho + zz.0 %*% psi2))
  cause.0 <- 1 + rbinom(n, 1, p2.0)
  ev1.0 <- 2 - cause.0
  
  # Event times, obtained from conditional CIFs (on event type, covariates, exposure)
  u.0 <- runif(n, 0, 1)
  aa.0 <- 1 - exp(gamma * exp(zz.0 %*% psi1)/(rho + zz.0 %*% psi2))
  bb.0 <- (log(1 - u.0*aa.0) * (rho + zz.0 %*% psi2))/(gamma * exp(zz.0 %*% psi1))
  t.0 <- (ev1.0) * (log(1 - bb.0))/(rho + zz.0 %*% psi2) + (1-ev1.0)*-(log(1 - u.0))/exp(zz.0 %*% psi3)
  
  
  # ---- Combine and derive RMTL, make figure ---- 
  dat.ate <- rbind(cbind(a.1, z1, z2, z3, z4, z5, t.1, cause.1), cbind(a.0, z1, z2, z3, z4, z5, t.0, cause.0))
  colnames(dat.ate)[c(1,7:8)] <- c('a', 't', 'cause')
  summary(dat.ate)
  
  dat.ate <- data.frame(dat.ate)
  
  rmtl.plot(times=dat.ate$t, event=dat.ate$cause, group=as.factor(dat.ate$a), tau=1, yaxismax=.6, title=title)  
}


cif.plot(n=n.iter, gamma=1.2, rho=-1.5, psi11=-0.3, psi21=0, title='Proportional subdistribution hazards')
cif.plot(n=n.iter, gamma=1.5, rho=-1.5, psi11=0.5, psi21=-2, title='Non-proportional subdistribution hazards')

# dev.off()


