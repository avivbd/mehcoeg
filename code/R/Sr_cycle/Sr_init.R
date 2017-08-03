
#list of carbon cycle parameters
q = list()

#data
q$Sr_data <- read_csv("~/Google Drive/Research/StanfordPostdoc/mehcoeg/Data/Sr_data.csv", col_types = 'dd')
q$Sr_spline <- read_csv("~/Google Drive/Research/StanfordPostdoc/mehcoeg/Data/Sr_spline.csv", col_types = 'dd')
q$Temp_data <- read_csv("~/Google Drive/Research/StanfordPostdoc/mehcoeg/Data/Temp_data.csv", col_types = 'dd')

#C-cycle params
source("./code/R/C_cycle/C_Cycle_init.R")

#time series attributes
q$t0 = q$Sr_spline$Age[1]
q$tf = tail(q$Sr_spline$Age,n = 1)
q$nt = 500
q$dt = (q$tf - q$t0)/q$nt
q$tt = matrix(seq(from=q$t0, to = q$tf, length.out = q$nt))
colnames(q$tt) = 'Age'

q$Sr_interp = matrix(nrow = q$nt, ncol = 2)
colnames(q$Sr_interp) = c('Age', 'RSr')
inds = findInterval(q$Sr_data$AGE, q$tt)
q$Sr_interp[inds, 'Age'] = q$Sr_data$AGE
q$Sr_interp[inds, 'RSr'] = q$Sr_data$RSr

q$Temp_interp = matrix(nrow = q$nt, ncol = 2)
colnames(q$Temp_interp) = c('Age', 'Temp')
inds = findInterval(q$Temp_data$AgeTemp, q$tt)
q$Temp_interp[inds, 'Age'] = q$Temp_data$AgeTemp
q$Temp_interp[inds, 'Temp'] = q$Temp_data$Temp

#reservoirs
q$MSr = 1.27e5 # Tmol, from Kump 87.
# results in a residence time of 2.2 Myr

#Sr fluxes
q$F_Sr_hydro_in = 2e-2# From Kump and Arthur 1997 
q$F_Sr_hydro_out = 2e-2
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

q$M0 = c(p$MCss, q$MSr,  q$RSr, 1)

#forcing
# q$f <- splinefun(q$Sr_spline$Age, q$Sr_spline$RSr)
# q$f = function(x){1}

#State variances
q$Q = diag(c(0, 0, 1e-9, 1e-1))


#observation errors
q$sigma_RSr = 1e-7
q$R = matrix(q$sigma_RSr^2)


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
