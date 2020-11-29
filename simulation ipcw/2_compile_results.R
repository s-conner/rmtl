### Assess simulations of IPCW regression method and output summary measures
### Sarah Conner, Updated October 22 2020

n <- c(rep(500, 8), rep(1000, 8))
beta1 <- -.15
beta2 <--.1
ph <- rep(c(1, 1, 0, 0), 4)
prop.exposed <-  c(rep(.5, 4), rep(.25, 4), rep(.5, 4), rep(.25, 4))
cens <- rep(c(.1, .25), 8)
bias <- rel.bias <- rmse <- cov <- avg.se <- emp.se <- rel.se <- matrix(NA, nrow=16, ncol=2)

nsim <- 6000

for (i in 1:16){
  
  beta1.i <- read.csv(paste0('output//beta1_', i, '.csv'), sep=',')
  beta2.i <- read.csv(paste0('output//beta2_', i, '.csv'), sep=',')
  
  bias[i, 1] <- mean(beta1.i$rmtldiff) - beta1
  bias[i, 2] <- mean(beta2.i$rmtldiff) - beta2
  
  rel.bias[i, 1] <- bias[i, 1]/beta1
  rel.bias[i, 2] <- bias[i, 2]/beta2

  rmse[i, 1] <- sqrt(mean((beta1.i$rmtldiff - beta1)^2))
  rmse[i, 2] <- sqrt(mean((beta2.i$rmtldiff - beta2)^2))
  
  cov[i, 1] <- sum(ifelse(beta1.i$rmtldiff.cil<=beta1 & beta1.i$rmtldiff.ciu>=beta1, 1, 0))/nsim
  cov[i, 2] <- sum(ifelse(beta2.i$rmtldiff.cil<=beta2 & beta2.i$rmtldiff.ciu>=beta2, 1, 0))/nsim
  
  avg.se[i,1] <- sqrt(sum(beta1.i$rmtldiff.se^2)/nsim)
  avg.se[i,2] <- sqrt(sum(beta2.i$rmtldiff.se^2)/nsim)
  
  emp.se[i,1] <- sqrt(sum((beta1.i$rmtldiff - mean(beta1.i$rmtldiff))^2)/(nsim-1))
  emp.se[i,2] <- sqrt(sum((beta2.i$rmtldiff - mean(beta2.i$rmtldiff))^2)/(nsim-1))
  
  rel.se[i,1] <- avg.se[i,1]/emp.se[i,1]
  rel.se[i,2] <- avg.se[i,1]/emp.se[i,1]
  
}


allres <- cbind(bias,rel.bias,rmse,cov,rel.se)

colnames(allres) <- c('Bias Beta1', 'Bias Beta2', 
                      'RelBias Beta1', 'RelBias Beta2', 
                      'RMSE Beta1', 'RMSE Beta2', 
                      'Cov Beta1', 'Cov Beta2',
                      'Rel SE Beta1', 'Rel SE Beta2')

# Export
write.csv(cbind(n, ph, prop.exposed, cens, allres), file='allreg_results_6000.csv')


