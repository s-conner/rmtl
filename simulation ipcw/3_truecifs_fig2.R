### Code to replicate Fig 2 (true CIFs in simulation study of IPCW model)
### Last updated by Sarah Conner on October 22 2020

library(dichromat)

rmtl <- function(gamma, rho){
  integrand <- function(t) {1 - exp(gamma * ((1 - exp(rho*t))/(rho)))}
  rmtl <- integrate(integrand, lower = 0, upper = 1)
  return(cbind(gamma, rho, rmtl=rmtl$value))
}

cif <- function(gamma, rho){
  cif <- 1 - exp(gamma * ((1 - exp(rho*t))/(rho)))
  return(cif)
}

t <- seq(0, 1, by=.01)

palette <- colorschemes$BrowntoBlue.10[c(2,4,9,10)]

pdf('cifs_sim2_days.pdf', width=12, height=6)
par(oma = c(2,1,1,1), mfrow = c(1,2), mar = c(4,4,4,1))

plot(t, cif(1.213442, -1.7), type='l', ylim=c(0,1), lwd=3, col=palette[1], lty=1, xaxt='n', 
     ylab='Cumulative incidence', xlab='Time (days)', main='Proportional subdistribution hazards')
lines(t, cif(1.770851, -1.7), type='l', lwd=3, col=palette[2], lty=2)
lines(t, cif(2.096025, -1.7), type='l', lwd=3, col=palette[3], lty=4)
lines(t, cif(2.879494, -1.7), type='l', lwd=3, col=palette[4], lty=3)
axis(1, at=c(0,60,120,180,240,300,365)/365, labels=c(0,60,120,180,240,300,365))
abline(v=1, col='grey', lty=3, lwd=3)

plot(t, cif(1.546588, -2.8), type='l', ylim=c(0,1), lwd=3, col=palette[1], lty=1, xaxt='n', 
     ylab='Cumulative incidence', xlab='Time (days)', main='Non-proportional subdistribution hazards')
lines(t, cif(2.291237, -2.9), type='l', lwd=3, col=palette[2], lty=2)
lines(t, cif(1.952647, -1.4), type='l', lwd=3, col=palette[3], lty=4)
lines(t, cif(2.879494, -1.7), type='l', lwd=3, col=palette[4], lty=3)
axis(1, at=c(0,60,120,180,240,300,365)/365, labels=c(0,60,120,180,240,300,365))
abline(v=1, col='grey', lty=3, lwd=3)

par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
plot(0, 0.2, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom', horiz=TRUE, col=c(palette[4],palette[3],palette[2],palette[1]), lty=c(3,4,2,1), lwd=rep(2,4), cex=.8,
       legend=c('Z1=0, Z2=0', 'Z1=0, Z2=1', 'Z1=1, Z2=0','Z1=1, Z2=1'))

dev.off()




