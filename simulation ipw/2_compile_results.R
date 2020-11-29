### Assess simulations of IPW method and output summary measures
### Sarah Conner, Updated October 22 2020

# 1-4 are PH (RMTL1-RMTL0 = -0.062), 5-8 are nonPH (RMTL1-RMTL0 = 0.022)

n <- c(rep(500, 8), rep(1000, 8))
truth <- c(rep(-0.062, 4), rep(0.022, 4), rep(-0.062, 4), rep(0.022, 4))
ph <- c(rep(1, 4), rep(0, 4), rep(1, 4), rep(0, 4))
prop.exposed <-  rep(c(.25, .25, .5, .5), 4)
cens <- rep(c(.1, .25), 8)
bias <- rel.bias <- rmse <- cov <- avg.se <- emp.se <- rel.se <- matrix(NA, nrow=16, ncol=1)

for (i in 1:16){
  simres_est <- read.csv(paste0('sim_estandtrue_newvar_output//simresults_estweights_newvar_', i, '.csv'), sep=',')
  truth.i <- truth[i]
  
  bias[i] <- mean(simres_est$rmtldiff*365) - truth.i*365
  rel.bias[i] <- bias[i]/truth.i
  rmse[i] <- sqrt(mean((simres_est$rmtldiff*365 - truth.i*365)^2))
  cov[i] <- sum(ifelse(simres_est$rmtldiff.cil<=truth.i & simres_est$rmtldiff.ciu>=truth.i, 1, 0))/6000
  avg.se[i] <- sqrt(sum(simres_est$rmtldiff.se^2)/6000)
  emp.se[i] <- sqrt(sum((simres_est$rmtldiff - mean(simres_est$rmtldiff))^2)/5999)
  rel.se[i] <- avg.se[i,1]/emp.se[i,1]

}

all.res <- cbind(bias, rel.bias, rmse, cov, avg.se, emp.se, rel.se)

colnames(all.res) <- c('Bias est', 
                       'Rel bias est',
                        'RMSE est', 
                        'Cov est', 
                        'Avg SE est', 
                        'Emp SE est',
                        'Rel SE est')
summary(all.res)

cbind(n, ph, prop.exposed, cens, all.res)


# Export
write.csv(cbind(n, ph, prop.exposed, cens, all.res), file='all_ipw_results.csv')

