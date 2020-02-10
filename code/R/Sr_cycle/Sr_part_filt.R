
library("MASS")
library(mvtnorm)
library(parallel)
library(scales)
library(readr)
require(stats)
library(parallel)
library(magrittr)

dev.off()
rm(list = ls())
# set.seed(1976)

source("./code/R/Sr_cycle/Sr_init.R")
source("./code/R/C_cycle/C_Cycle_init.R")
source("./code/R/General/PF.R")

p = C_Cycle_init()
q = Sr_init('y')

state_eqs <- function (x, k)
    t = q$t[1] + k*q$dt*1e-6
    MC <- x[1]
    MSr <- x[2]
    RSr <- x[3]
    F_forcing_v = x[4]
    F_forcing_sw = x[5]
    
    # F_forcing_sw = predict(q$ssp, t)$y
    Fvolc = p$Fv*F_forcing_v
    Fwsil = p$kws*MC*F_forcing_sw
    
    
    # print(Fwsil)

    dMC = Fvolc + p$Fwo - Fwsil - p$Fbo
    
    dMSr = q$F_Sr_hydro_in + p$Fwcarb*q$k_Sr_Fwcarb + 
            Fwsil*q$k_Sr_Fwsil - q$F_Sr_hydro_out - q$k_Sr_bcarb*MSr
    
    dRSr = (    q$F_Sr_hydro_in*(q$R_Sr_hydro_in - RSr) +
                p$Fwcarb*q$k_Sr_Fwcarb*(q$R_Sr_wcarb  - RSr) +
                Fwsil*q$k_Sr_Fwsil*(q$R_Sr_wsil - RSr))/MSr
    
    # dMC = 0
    # dMSr = 0
    # dRSr = 0
    dForcing_v = 0
    dForcing_sw = 0
    
    return(c((MC + q$dt*dMC), (MSr + q$dt*dMSr), (RSr + q$dt*dRSr), 
             (F_forcing_v + q$dt*dForcing_v), (F_forcing_sw + q$dt*dForcing_sw)) )}

obs_eqs <- function (x, k){
    MC <- x[1]
    RSr <- x[3]
    F_forcing_sw <- x[5]
    RCO2 = q$co2_conv(MC)
    return(cbind(RSr, RCO2, F_forcing_sw))
}


#TODO put in a try catch box and restart if fails. 
# zz = log10(c(1e7, 1e7, 1e-7, .05))
# best <- optim(zz, optim_fun, method = 'SANN',
               # control = list(maxit = 20000, temp = 20))


# zz = c(1e10, 5e6, 1e-7, 1e-6, 1e-1)
zz = c(1e10, 5e6, 1e-7, 1e-2, 5e-1)
mod_pars <- list(m0=c(p$MCss, q$MSr,  q$RSr, 1, 1), #initial state estimates
                 C0=diag(1*zz),   #uncertainty of initial state estimates
                 V=diag(c(1e-3, 5e-1, 1e-3)),    #measurement noise
                 W=diag(zz))    #process noise

N_particles = 25000

# #about 1 minute
pf_results = PF(y=cbind(q$y_synth, q$pCO2_noisy, q$F_Forcing_ws_nsy), 
                mod=mod_pars, N=N_particles, 
                Nthreshold=floor(0.75*N_particles),
                roughening=TRUE, GGfunction=state_eqs, 
                FFfunction=obs_eqs, MCparticles=TRUE)



tt = seq(q$t[1] - q$dt*1e-6, q$t %>% tail(.,1), by = q$dt*1e-6 )
l = length(tt)

results_to_mat <- function(pf_results, i){
    matrix(unlist(lapply(pf_results$xpt, function(x) x[,i])), nrow=l, byrow=TRUE)}

MC_filt = results_to_mat(pf_results, 1)
MSr_filt = results_to_mat(pf_results, 2)
RSr_filt = results_to_mat(pf_results, 3)
F_forcing_v_filt = results_to_mat(pf_results, 4)
F_forcing_ws_filt = results_to_mat(pf_results, 5)


#95% ci for x
ci <- function(pf_results, i){
    se = sqrt(unlist(lapply(pf_results$C, function(x) diag(x)[i])))
    l = vector('list')
    l$ci95 = pf_results$m[,i] + qnorm(.05/2)*se%o%c(1, -1)
    l$ci65 = pf_results$m[,i] + qnorm(.35/2)*se%o%c(1, -1)
    return(l)
}

ci_MC <- ci(pf_results, 1)
ci_MSr <- ci(pf_results, 2)
ci_RSr <- ci(pf_results, 3)
ci_FFv <- ci(pf_results, 4)
ci_FFws <- ci(pf_results, 5)


# layout(rbind(c(1,1, 2,2, 3,3), c(4,4,4, 5,5,5)))
par(mfrow=c(2,3))

plot(q$t, q$x_synth[, 1], ylab = 'MC', xlab='Age', ylim=(range(ci_MC$ci95)))
lines(tt, pf_results$m[,1], col='blue')
matlines(tt, ci_MC$ci95, lty=3, col="gray")


plot(q$t, q$x_synth[, 2], ylab = 'MSr', xlab='Age', ylim=range(ci_MSr$ci95))
lines(tt, pf_results$m[,2], col='blue')
matlines(tt, ci_MSr$ci95, lty = 3, col = 'gray')


plot(q$t, q$F_Forcing_v, ylab = 'F_Forcing_v', xlab='Age', ylim = range(ci_FFv$ci95))
lines(tt, pf_results$m[,4], col='blue')
matlines(tt, ci_FFv$ci95, lty = 3, col = 'gray')

plot(q$t, q$F_Forcing_ws_nsy, ylab = 'F_Forcing_ws', xlab='Age', ylim = range(ci_FFws$ci95))
lines(tt, pf_results$m[,5], col='blue')
matlines(tt, ci_FFws$ci95, lty = 3, col = 'gray')

plot(q$t, q$y_synth, ylab = 'RSr', ylim=(range(ci_RSr$ci95)))
lines(tt, pf_results$m[,3], col='blue')
matlines(tt, ci_RSr$ci95, lty=3, col="gray")


plot(q$t, q$pCO2_noisy, ylab = 'Temp', xlab='Age', ylim = range(q$co2_conv(ci_MC$ci95)))
lines(tt, q$co2_conv(pf_results$m[,1]), col='blue')
matlines(tt, q$co2_conv(ci_MC$ci95), lty=3, col="gray")



#---------

#plot full posterior distributions

n_trajectories = 20
pf_smooth = replicate(n_trajectories, PFsmooth(pf_results) )

MC_smooth = matrix(0, ncol=n_trajectories, nrow=l)
MSr_smooth = matrix(0, ncol=n_trajectories, nrow=l)
RSr_smooth = matrix(0, ncol=n_trajectories, nrow=l)
F_forcing_smooth = matrix(0, ncol=n_trajectories, nrow=l)

for (i in 1:n_trajectories){
    MC_smooth[,i] = pf_smooth[i]$s[,1]
    MSr_smooth[,i] = pf_smooth[i]$s[,2]
    RSr_smooth[,i] = pf_smooth[i]$s[,3]
    F_forcing_smooth[,i] = pf_smooth[i]$s[,4]
}


layout(rbind(c(1,1, 2,2, 3,3), c(4,4,4, 5,5,5)))

plot(q$t, q$x_synth[, 1], ylab = 'MC', xlab='Age')
# lines(tt, pf_results$m[,1], col='blue')
lines(tt, MC_smooth, col='red')


plot(q$t, q$x_synth[, 2], ylab = 'MSr', xlab='Age')
# lines(tt, pf_results$m[,2], col='blue')
lines(tt, MSr_smooth, col='red')

plot(q$t, q$F_Forcing_v, ylab = 'F_Forcing_v', xlab='Age')
# lines(tt, pf_results$m[,4], col='blue')
lines(tt, F_forcing_smooth, col='red')

plot(q$t, q$y_synth, ylab = 'RSr')
# lines(tt, pf_results$m[,3], col='blue')
lines(tt, RSr_smooth, col='red')

plot(q$t, q$pCO2_noisy, ylab = 'Temp', xlab='Age', ylim = c(16,40))
# lines(tt, q$TD*log2((pf_results$m[,1]/p$MCss)^2) + 20 , col='blue')
lines(tt, q$TD*log2((MC_smooth/p$MCss)^2) + 20, col='red' )









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





