library("MASS")
library(mvtnorm)
library(parallel)
library(scales)

graphics.off()
rm(list = ls())

set.seed(1976)


source("./code/R/C_cycle/C_Cycle_init.R")
source("./code/R/General/PF.R")


# preliminaries ----

#set up for parallel computation
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)
#source the aux code for the cluster
clusterCall(cl, function() { source(("./PF.R")) })

#define model and observation equations with a known forcing

state_eqs0 <- function (x, k){
  MP <- x[1] 
  MC <- x[2]
  delC <- x[3]
  F_forcing <- x[4]
  c(
    MP + p$dt*(p$kwp*MC - p$kbp*MP),
    MC + p$dt*(p$Fv + p$Fwo - p$kws*MC - p$kbo*p$F_forcing(k*p$dt)*MP),
    delC + p$dt*((p$Fv*(p$d13C_volc - delC) + p$Fwo*(p$delwo - delC) + 
                    p$Fwcarb*(p$delC - delC) + p$kbo*MP*p$eps )/MC),
    1
  )}

#define equations that determine observability
obs_eqs <- function (x, k){
  MP = x[1]
  MC = x[2]
  delC = x[3]
  F_forcing = x[4]
  return(delC) 
}


W = function(){mvrnorm(n = 1, mu = matrix(0,1,dim(p$Q)[2]), Sigma = p$Q_sample)}
V = function(){mvrnorm(n = 1, mu = matrix(0,1,dim(p$R)[2]), Sigma = p$R_sample)}

x = matrix(0, ncol=dim(p$Q)[2], nrow=p$nt)
y = matrix(0, ncol=dim(p$R)[2], nrow=p$nt)

x[1,] = p$M0

#make some synthetic data
for (k in 1:(p$nt-1)){
  x[k+1,] = state_eqs0(x[k,], k) + W()
  y[k+1] = obs_eqs(x[k+1,], k)  + V()
}


# #visualize states and observations
# plot(p$tt/1e6, y, type = "p", xlab = "Time [Myr]", ylab = expression(paste("δ"^"13"*"C" )) )
# plot(p$tt/1e6, x[,1], "b", xlab = "Time [Myr]", ylab = "MP")
# plot(p$tt/1e6, x[,2], "b", xlab = "Time [Myr]", ylab = "MC")
# plot(p$tt, p$F_forcing(p$tt), "b", xlab = "Time [Myr]", ylab = "Forcing")


#apply the particle filter ----

#redefine the state equations and now let the forcing be an unkown
state_eqs <- function (x, k){
  MP <- x[1] 
  MC <- x[2]
  delC <- x[3]
  F_forcing <- x[4]
  c(
    MP + p$dt*(p$kwp*MC - p$kbp*MP),
    MC + p$dt*(p$Fv + p$Fwo - p$kws*MC - p$kbo*F_forcing*MP),
    delC + p$dt*((p$Fv*(p$d13C_volc - delC) + p$Fwo*(p$delwo - delC) + 
                    p$Fwcarb*(p$delC - delC) + p$kbo*MP*p$eps )/MC),
    F_forcing
  )}


#particle filter 
mod_pars <- list(m0=p$M0, #initial state estimates
                C0=p$Q,   #uncertainty of initial state estimates
                V=p$R,    #measurement noise
                W=p$Q)    #process noise




pf_results = PF(y=y, mod=mod_pars, N=p$n_ensemble, resampling="strat", Nthreshold=p$n_ensemble/2,
   roughening=TRUE, Grough=.2, 
   GGfunction=state_eqs, FFfunction=obs_eqs, MCparticles=TRUE)



# Apply particle smoother ----

#run smoother in serial
# pf_smooth <- replicate(p$n_trajectories, PFsmooth(pf_results))

#run smoother in parallel for a bit of speedup
parReplicate <- function(cl, n, expr, simplify=TRUE, USE.NAMES=TRUE){
  parSapply(cl, 1:n, function(i, ex) eval(ex, envir=.GlobalEnv),
            substitute(expr), simplify=simplify, USE.NAMES=USE.NAMES)
}

clusterExport(cl,c("pf_results","p", "PFsmooth"))

pf_smooth = parReplicate(cl, p$n_trajectories, PFsmooth(pf_results) )


MP_smooth = matrix(0, ncol=p$n_trajectories, nrow=p$nt+1)
MC_smooth = matrix(0, ncol=p$n_trajectories, nrow=p$nt+1)
delC_smooth = matrix(0, ncol=p$n_trajectories, nrow=p$nt+1)
F_forcing_smooth = matrix(0, ncol=p$n_trajectories, nrow=p$nt+1)

for (i in 1:p$n_trajectories){
  MP_smooth[,i] = pf_smooth[i]$s[,1]
  MC_smooth[,i] = pf_smooth[i]$s[,2]
  delC_smooth[,i] = pf_smooth[i]$s[,3]
  F_forcing_smooth[,i] = pf_smooth[i]$s[,4]
}



# plotting ----

plot(p$tt/1e6, x[,1], type="p", pch=16, col="black", xlab = "Time [Myr]", ylab = "MP")
matlines(c(0, p$tt)/1e6, MP_smooth, col= gray(0.5, 0.8), lty=3 )
lines(c(0,p$tt)/1e6, pf_results$m[,1], lty=2, col="blue", lwd=3)
lines(c(0, p$tt)/1e6, apply(MP_smooth, 1, median), lty=2, col="red", lwd=3)
legend("topleft", col=c("black", "blue", "red", "gray"),
       lty=c(NA, 2, 2, 2), pch=c(16, NA, NA, NA), lwd=c(NA, 2, 2, 1),
       legend=c("true state", "filtered PF", "smoothed PF", "particle trajectories"),
       bty="n", y.intersp=1.2, cex=.7)



plot(p$tt/1e6, x[,2], type="p", pch=16, col="black", xlab = "Time [Myr]", ylab = "MC")
lines(c(0,p$tt)/1e6, pf_results$m[,2], lty=2, col="blue", lwd=3)
matlines(c(0, p$tt)/1e6, MC_smooth, col= gray(0.5, 0.8), lty=3 )
lines(c(0, p$tt)/1e6, apply(MC_smooth, 1, median), lty=2, col="red", lwd=3)
legend("topleft", col=c("black", "blue", "red", "gray"),
       lty=c(NA, 2, 2, 2), pch=c(16, NA, NA, NA), lwd=c(NA, 2, 2, 1),
       legend=c("true state", "filtered PF", "smoothed PF", "particle trajectories"),
       bty="n", y.intersp=1.2, cex=.7)


plot(p$tt/1e6, p$F_forcing(p$tt), type="p", pch=16, col="black", xlab = "Time [Myr]", ylab = "Forcing")
lines(c(0,p$tt)/1e6, pf_results$m[,4], lty=2, col="blue", lwd=3)
matlines(c(0, p$tt)/1e6, F_forcing_smooth, col= gray(0.5, 0.8), lty=3 )
lines(c(0, p$tt)/1e6, apply(F_forcing_smooth, 1, median), lty=2, col="red", lwd=3)
legend("topleft", col=c("black", "blue", "red", "gray"),
       lty=c(NA, 2, 2, 2), pch=c(16, NA, NA, NA), lwd=c(NA, 2, 2, 1),
       legend=c("true state", "filtered PF", "smoothed PF", "particle trajectories"),
       bty="n", y.intersp=1.2, cex=.7)


plot(p$tt/1e6, y, type="p", pch=16, col="black", xlab = "Time [Myr]", ylab = "Forcing")
lines(c(0,p$tt)/1e6, pf_results$m[,3], lty=2, col="blue", lwd=3)
matlines(c(0, p$tt)/1e6, delC_smooth, col= gray(0.5, 0.8), lty=3 )
lines(c(0, p$tt)/1e6, apply(delC_smooth, 1, median), lty=2, col="red", lwd=3)
legend("topleft", col=c("black", "blue", "red", "gray"),
       lty=c(NA, 2, 2, 2), pch=c(16, NA, NA, NA), lwd=c(NA, 2, 2, 1),
       legend=c("true state", "filtered PF", "smoothed PF", "particle trajectories"),
       bty="n", y.intersp=1.2, cex=.7)







