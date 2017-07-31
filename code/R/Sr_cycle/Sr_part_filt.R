library("MASS")
library(mvtnorm)
library(parallel)
library(scales)
library(readxl)
require(stats)

graphics.off()
rm(list = ls())

set.seed(1976)

source("./code/R/Sr_cycle/Sr_init.R")
source("./code/R/C_cycle/C_Cycle_init.R")
source("./code/R/General/PF.R")


plot(q$Sr_data$AGE, q$Sr_data$RSr, 'p', xlab = 'Age [Ma]', ylab = expression(Sr^87/Sr^86 ) )
lines(q$Sr_spline$Age, q$Sr_spline$RSr, col='red')    
legend('bottomright', c('Data', 'Smoothing Spline'), lty = c(NA, 1), pch = c(1, NA), col = c('black', 'red'))
title('Sr data')


state_eqs0 <- function (x, k){
    
    
    MC <- x[1] 
    MSr <- x[2]
    RSr <- x[3]
    F_forcing <- x[4]
    
    
    RCO2 = (MC/p$MCss)^2
    Fvolc = p$Fv*q$f(q$tt[k])
    Fwsil = p$Fws*(log2(RCO2)+1)
    Fborg = p$kbo*p$MPss
    Fwcarb = p$Fwcarb
    
    dMC = Fvolc + p$Fwo - Fwsil - Fborg
    dMSr = q$F_Sr_hydro_in + Fwcarb*q$k_Sr_Fwcarb + Fwsil*q$k_Sr_Fwsil
              - q$F_Sr_hydro_out - q$k_Sr_bcarb*MSr
    
    dRSr = (q$F_Sr_hydro_in*(q$R_Sr_hydro_in - RSr) + p$Fwcarb*q$k_Sr_Fwcarb*(q$R_Sr_wcarb  - RSr) +
             Fwsil*(q$R_Sr_wsil - RSr))/MSr
    
    
    return(
        c(
        (MC + q$dt*dMC),
        (MSr + q$dt*dMSr),
        (RSr + q$dt*dRSr),
        1)        )}


#define equations that determine observability
obs_eqs <- function (x, k){
    MC <- x[1] 
    MSr <- x[2]
    RSr <- x[3]
    F_forcing <- x[4]
    
    RCO2 = (MC/p$MCss)^2
    return(RSr) 
}



#noise
W = function(){mvrnorm(n = 1, mu = matrix(0, 1, dim(q$Q)[2]), Sigma = q$Q_sample)}
V = function(){mvrnorm(n = 1, mu = matrix(0, 1, dim(q$R)[2]), Sigma = q$R_sample)}

#preallocate
x = matrix(0, ncol=dim(p$Q)[2], nrow=p$nt)
y = matrix(0, ncol=dim(p$R)[2], nrow=p$nt)

#initial values
x[1,] = c(p$MCss, q$MSr,  q$RSr, 1)
y[1,] = q$RSr

#make some synthetic data
for (k in 1:(q$nt-1)){
    x[k+1,] = state_eqs0(x[k,], k) #+ W()
    y[k+1] = obs_eqs(x[k+1,], k)  #+ V()
}


plot(q$tt, y)
plot(q$tt, x[,1], 'l')
plot(q$tt, x[,2], 'l')
plot(q$tt, x[,3], 'l')
plot(q$tt, x[,4], 'l')

# dC <- f.wcarb - k*(log2(y[1]/C.i)+1)
# dSr <- f.w.hyd.sr + k.sr.w*f.wcarb - f.b.hyd.sr - k.sr.b*k*(log2(y[1]/C.i)+1)
# dSr87 <- (f.w.hyd.sr*(R.Sr.hyd - Sr87) + k.sr.w*k*(log2(y[1]/C.i)+1)*(R.Sr.w - y[3]))/y[2]




