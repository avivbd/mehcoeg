library("MASS")
library(mvtnorm)
library(parallel)
library(scales)
library(readxl)
require(stats)
library(parallel)

dev.off()
rm(list = ls())

set.seed(1976)

source("./code/R/Sr_cycle/Sr_init.R")
source("./code/R/C_cycle/C_Cycle_init.R")
source("./code/R/General/PF.R")



#make synthetic data
state_eqs0 <- function (x, k, Forcing, dt){

    MC <- x[1]
    MSr <- x[2]
    RSr <- x[3]

    Fvolc = p$Fv*Forcing[k]
    Fwsil = p$kws*MC

    dMC = Fvolc + p$Fwo - Fwsil - p$Fbo
    
    dMSr = q$F_Sr_hydro_in + p$Fwcarb*q$k_Sr_Fwcarb + Fwsil*q$k_Sr_Fwsil - q$F_Sr_hydro_out - q$k_Sr_bcarb*MSr
    
    print(q$k_Sr_bcarb*MSr)
    
    dRSr = (    q$F_Sr_hydro_in*(q$R_Sr_hydro_in - RSr) + 
                p$Fwcarb*q$k_Sr_Fwcarb*(q$R_Sr_wcarb  - RSr) +
                Fwsil*q$k_Sr_Fwsil*(q$R_Sr_wsil - RSr))/MSr

    # print(dMC)
    # print(dMSr)
    return(c((MC + dt*dMC),
             (MSr + dt*dMSr),
             (RSr + dt*dRSr)) )}


#define equations that determine observability
obs_eqs <- function (x, k){
    MC <- x[1]
    MSr <- x[2]
    RSr <- x[3]
    return(RSr) 
}

l = length(q$Sr_spline$Age)
x = matrix(ncol=3, nrow=l)
x[1,] = c(p$MCss, q$MSr,  q$Sr_spline$RSr[1])
#make some synthetic data
z = 1 + 0.25*(q$Sr_spline$RSr - min(q$Sr_spline$RSr))/sd(q$Sr_spline$RSr - min(q$Sr_spline$RSr))
for (k in 1:(l-1) ){
    x[k+1,] = state_eqs0(x[k,], k, z, mean(diff(q$Sr_spline$Age))*1e6 ) 
    
}


par(mfrow=c(2,2))
plot(q$Sr_spline$Age, z, ylab = 'Forcing', xlab='Age')
plot(q$Sr_spline$Age, x[, 1], ylab = 'MC', xlab='Age')
plot(q$Sr_spline$Age, x[, 2], ylab = 'MSr', xlab='Age')
plot(q$Sr_data$AGE, q$Sr_data$RSr, ylab = 'RSr', xlab='Age')
lines(q$Sr_spline$Age, x[, 3])















# obs_eqs <- function (x, k){
#     MC <- x[1] 
#     MSr <- x[2]
#     RSr <- x[3]
#     F_forcing <- x[4]
#     
#     RCO2 = (MC/p$MCss)^2
#     return(RSr) 
# }


q$M0 = c(q$MSr,  q$RSr, 1)
q$Q = diag(c(1e-15, 1e-15, 1e-2))
q$R = 1e-4

mod_pars <- list(m0=q$M0, #initial state estimates
                 C0=q$Q,   #uncertainty of initial state estimates
                 V=q$R,    #measurement noise
                 W=q$Q)    #process noise

#about 1 minute
pf_results = PF(y=q$Sr_interp[,2], mod=mod_pars, N=5000, Nthreshold=4000,
                roughening=TRUE, GGfunction=state_eqs, FFfunction=obs_eqs, MCparticles=TRUE)

# pf_smooth = PFsmooth(pf_results)

# start.time <- Sys.time()
# pf_smooth = replicate(5, PFsmooth(pf_results) )
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# print(time.taken)

# clusterExport(cl,c("pf_results","p", "PFsmooth"))
# 
# pf_smooth = parReplicate(cl, 5, PFsmooth(pf_results) )
# 
# 
# MC_smooth = matrix(0, ncol=p$n_trajectories, nrow=p$nt+1)
# MSr_smooth = matrix(0, ncol=p$n_trajectories, nrow=p$nt+1)
# RSr_smooth = matrix(0, ncol=p$n_trajectories, nrow=p$nt+1)
# F_forcing_smooth = matrix(0, ncol=p$n_trajectories, nrow=p$nt+1)
# 
# for (i in 1:p$n_trajectories){
#     MC_smooth[,i] = pf_smooth[i]$s[,1]
#     MSr_smooth[,i] = pf_smooth[i]$s[,2]
#     RSr_smooth[,i] = pf_smooth[i]$s[,3]
#     F_forcing_smooth[,i] = pf_smooth[i]$s[,4]
# }
# 


# par(mfcol=c(1,2))
plot(q$Sr_data$AGE, q$Sr_data$RSr, 'p', xlab = 'Age [Ma]', ylab = expression(Sr^87/Sr^86 ) )
# plot(q$Temp_data$AgeTemp, q$Temp_data$Temp, 'p', xlab = 'Age [Ma]', 
     # ylab = ('degrees C' ), xlim = par("usr")[1:2] )
t = c(q$Sr_interp[1,1], q$Sr_interp[,1]) 
points(t, pf_results$m[,2], col="red", pch=4)
# lines(t, pf_smooth, lty = 2, col='red')


plot(t, pf_results$m[,1], lty=2, col="blue", lwd=3,
     ylab = expression(MSr), xlab='Age')

plot(t, pf_results$m[,3], lty=2, col="blue", lwd=3, 
     ylab = expression(Forcing), xlab='Age')





state_eqs1 <- function (x, k, F_forcing){
    
    MSr <- x[1]
    RSr <- x[2]
    
    Fwsil = p$Fws*F_forcing
    
    # dMC = Fvolc + p$Fwo - Fwsil - Fborg
    dMSr = q$F_Sr_hydro_in + p$Fwcarb*q$k_Sr_Fwcarb + 
        Fwsil*q$k_Sr_Fwsil - q$F_Sr_hydro_out - q$k_Sr_bcarb*MSr
    
    dRSr = (q$F_Sr_hydro_in*(q$R_Sr_hydro_in - RSr) + p$Fwcarb*q$k_Sr_Fwcarb*(q$R_Sr_wcarb  - RSr) +
                Fwsil*q$k_Sr_Fwsil*(q$R_Sr_wsil - RSr))/MSr
    
    print(dMSr)
    
    return(
        c(  (MSr + q$dt*dMSr),
            (RSr + q$dt*dRSr))        )}

l = length(pf_results$m[,3])
           
xx = matrix(0, ncol=2, nrow=l)
xx[1,] = q$M0[1:2]

for (k in 1:(l-1)){
    xx[k+1,] = state_eqs1(xx[k,], k, pf_results$m[k,3]) 
}




# plot(t, pf_results$m[,4], lty=2, col="blue", lwd=3, 
#      ylab = expression(Forcing), xlab='Age')

# legend('bottomright', c('Data', 'Smoothing Spline'), lty = c(NA, 1), pch = c(1, NA), col = c('black', 'red'))
# title('Sr data')

# graphics.off()
# plot(q$tt, x[,1], 'l')
# plot(q$tt, x[,2], 'l')
# plot(q$tt, x[,3], 'l')
# plot(q$tt, x[,4], 'l')
# plot(q$tt, y, 'l')

