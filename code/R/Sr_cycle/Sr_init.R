
Sr_init = function(plot_='n'){
    #list of carbon cycle parameters
    q = list()
    
    #data
    q$Sr_data <- read_csv("./Data/Sr_data.csv", col_types = 'dd')
    q$Sr_spline <- read_csv("./Data/Sr_spline.csv", col_types = 'dd')
    q$Temp_data <- read_csv("./Data/Temp_data.csv", col_types = 'dd')
    
    #C-cycle params
    source("./code/R/C_cycle/C_Cycle_init.R")
    p = C_Cycle_init()
    
    #reservoirs
    q$MSr = 1.27e5 # Tmol, from Kump 87.
    # results in a residence time of 2.2 Myr
    
    #Sr fluxes
    q$F_Sr_hydro_in = 2.5e-2# From Kump and Arthur 1997 
    q$F_Sr_hydro_out = 2.5e-2
    q$F_Sr_wcarb = 2.5e-2 #from K+A 1997
    q$F_Sr_wsil  = 0.25*q$F_Sr_wcarb #from K+A 1997
    q$F_Sr_bcarb = q$F_Sr_wcarb + q$F_Sr_wsil
    
    #Sr values
    q$R_Sr_hydro_in = 0.7035
    q$R_Sr_wsil = 0.7164 #
    q$R_Sr_wcarb = 0.7080 # carbonate weathering from Kump and Arthur 1997
    
    q$RSr = (q$F_Sr_hydro_in*q$R_Sr_hydro_in + q$F_Sr_wsil*q$R_Sr_wsil + q$F_Sr_wcarb*q$R_Sr_wcarb)/
        (q$F_Sr_hydro_in+q$F_Sr_wcarb+q$F_Sr_wsil)
    
    #sensitivities
    q$k_Sr_Fwcarb = q$F_Sr_wcarb/p$Fwcarb
    q$k_Sr_Fwsil = q$F_Sr_wsil/p$Fws
    q$k_Sr_bcarb = q$F_Sr_bcarb/q$MSr
    
    
    #observation errors
    q$sigma_RSr = 1e-7
    q$R = matrix(q$sigma_RSr^2)
    
    
    #make synthetic data
    state_eqs0 <- function (x, k, F_Forcing_v, F_Forcing_ws, dt){
        
        MC <- x[1]
        MSr <- x[2]
        RSr <- x[3]
        
        Fvolc = p$Fv*F_Forcing_v
        Fwsil = p$kws*MC*F_Forcing_ws
        
        dMC = Fvolc + p$Fwo - Fwsil - p$Fbo
        
        dMSr = q$F_Sr_hydro_in + p$Fwcarb*q$k_Sr_Fwcarb + Fwsil*q$k_Sr_Fwsil - q$F_Sr_hydro_out - q$k_Sr_bcarb*MSr
        
        dRSr = (        q$F_Sr_hydro_in*(q$R_Sr_hydro_in - RSr) + 
                        p$Fwcarb*q$k_Sr_Fwcarb*(q$R_Sr_wcarb  - RSr) +
                        Fwsil*q$k_Sr_Fwsil*(q$R_Sr_wsil - RSr))/MSr
        
        
        return(c((MC + dt*dMC),
                 (MSr + dt*dMSr),
                 (RSr + dt*dRSr)) )}
    
    #make forcing for synthetic data
    q$n = 100
    q$t = seq(-257, -235, length.out = q$n )
    midp_v = 20
    q$TD = 3
    q$F_Forcing_v = #matrix(1, q$n)
        c(1 + 1.5/(1+exp(-10*(q$t[1:midp_v] + 254))),
          1.5 + 1/(1+exp(1*(q$t[(midp_v+1):q$n] + 249))))
    
    
    
    midp_ws = 20
    q$F_Forcing_ws = #matrix(1, q$n)
    c(1 - 1/(1.5 + exp(-10*(q$t[1:midp_ws] + 254))),
      1 - 1/(1.5 + exp(1*(q$t[(midp_ws+1):q$n] + 249))))
    
    q$F_Forcing_ws_nsy = q$F_Forcing_ws + rnorm(q$n, 0, 1e-1)
    q$ssp = smooth.spline(x=q$t, y=q$F_Forcing_ws)
    
    q$dt = mean(diff(q$t))*1e6
    q$x_synth = matrix(ncol=3, nrow=q$n)
    q$x_synth[1,] = c(p$MCss, q$MSr,  q$RSr[1])
    
    for (k in 1:(q$n-1) ){
        q$x_synth[k+1,] = state_eqs0(q$x_synth[k,], k, q$F_Forcing_v[k], q$F_Forcing_ws[k], q$dt) 
    }
    
    #produce noisy observations
    q$y_synth = matrix(q$x_synth[, 3] + rnorm(q$n, 0, 1e-4))
    q$co2_conv <- function(MC){q$TD*log2((MC/p$MCss)^2) + 20}
    q$pCO2 = q$co2_conv(q$x_synth[,1])
    q$pCO2_noisy = q$pCO2 + rnorm(q$n, 0, 2e0)
    
    if (plot_=='y'){
        
        #plot synthetic data and the states that result in it
        # layout(rbind(c(1,1, 2,2, 3,3), c(4,4,4, 5,5,5)))
        par(mfrow=c(2,3))
        
        plot(q$t, q$x_synth[, 1], ylab = 'MC', xlab='Age')
        plot(q$t, q$x_synth[, 2], ylab = 'MSr', xlab='Age')
        
        plot(q$t, q$F_Forcing_v, ylab = 'F_Forcing_v', xlab='Age')
        plot(q$t, q$F_Forcing_ws_nsy, ylab = 'F_Forcing_ws', xlab='Age')
        
        plot(q$t, q$x_synth[, 3], ylab = 'RSr', xlab='Age', type = 'l', ylim=c(0.7065, 0.7085))
        points(q$t, q$y_synth , cex=0.5, col='black')
        points(q$Sr_data$AGE, q$Sr_data$RSr, col='red')
    
        plot(q$t, q$pCO2_noisy, 'p', cex=0.5, col='black', ylim = c(16,40))
        lines(q$t, q$pCO2, col='black')
        points(q$Temp_data$AgeTemp, q$Temp_data$Temp, col='red')
        
        
    }
    
    
    
    return(q)
    
}










#set up for parallel computation

# # Calculate the number of cores
# no_cores <- detectCores() - 1
# 
# # Initiate cluster
# cl <- makeCluster(no_cores)
# 
# #source the aux code for the cluster
# clusterCall(cl, function() { source(("./PF.R")) })
# 
# #parallel replicate function
# parReplicate <- function(cl, n, expr, simplify=TRUE, USE.NAMES=TRUE){
#     parSapply(cl, integer(n), function(i, ex) eval(ex, envir=.GlobalEnv),
#               substitute(expr), simplify=simplify, USE.NAMES=USE.NAMES)
# }
# 

