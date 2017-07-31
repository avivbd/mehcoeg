
#list of carbon cycle parameters
q = list()

#data
q$Sr_data <- read_csv("~/Google Drive/Research/StanfordPostdoc/mehcoeg/Data/Sr_data.csv", col_types = 'dd')
q$Sr_spline <- read_csv("~/Google Drive/Research/StanfordPostdoc/mehcoeg/Data/Sr_spline.csv", col_types = 'dd')

#C-cycle params
source("./code/R/C_cycle/C_Cycle_init.R")

#time series attributes
q$t0 = q$Sr_spline$Age[1]
q$tf = tail(q$Sr_spline$Age,n = 1)
q$dt = mean(diff(q$Sr_spline$Age))
q$tt = q$Sr_spline$Age
q$nt = length(q$tt)

#reservoirs
q$MSr = 1.27e5 # Tmol, from Schaal disertation based on Kump 87.
# results in a residence time of 2.6 Myr

q$RSr = 0.7072

#Sr fluxes
q$F_Sr_hydro_in = 1e-2 # Note that is larger by a factor of 2.8 than Kump and Arthur 1997 (3.5e-3 Tmol/yr)
q$F_Sr_hydro_out = 1e-2
q$F_Sr_bcarb = 3.5e-2 #from Kump 87, compare to 2.2e-2 Tmol/yr @K+A
q$F_Sr_wcarb = 0.5e-2 #from Kump 87
q$F_Sr_wsil = 3e-2 #from SchKump 87aal

#Sr values
q$R_Sr_hydro_in = 0.7035
q$R_Sr_wsil = 0.712 #a bit higher than K+A which give 0.7164 for silicate weathering 
q$R_Sr_wcarb = 0.7080 # carbonate weathering from Kump and Arthur 1997

#sensitivities
q$k_Sr_Fwcarb = q$F_Sr_wcarb/p$Fwcarb
q$k_Sr_Fwsil = q$F_Sr_wsil/p$Fws
q$k_Sr_bcarb = q$F_Sr_bcarb/q$MSr

#forcing
q$f <- splinefun(q$Sr_spline$Age, q$Sr_spline$RSr - min(q$Sr_spline$RSr) + 1) 

#State variances
q$sigma_MSr = 0.075*q$MSr
q$sigma_RSr = 0.075*q$RSr

p$Q = diag(c(p$sigma_MC, q$sigma_MSr^2, q$sigma_RSr^2, 0.05))
p$Q_sample = diag(c(p$sigma_MC, q$sigma_MSr^2, q$sigma_RSr^2, 0.05))

#observation errors
q$sigma_RSr = 1*q$RSr
q$R = matrix(q$sigma_RSr^2)



