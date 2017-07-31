
#list of carbon cycle parameters
p = list()

#time series attributes
p$t0 = 0
p$tf = 6e6
p$dt = 5e4
p$tt = seq(from=p$t0, to=p$tf, by=p$dt )
p$nt = length(p$tt)

#reservoirs
p$MCss = 3.8e6 #Tmol
p$MPss = 2e3

#fluxes
p$total_input = 50 #Tmol/yr
p$fworg = 0.2
p$forg = 0.2
p$Fv = 5
p$Fwo = 9
p$Fwcarb = 36
p$Fwp = 3.6e-2
p$Fbp = 3.6e-2
p$pCO2 = 540 #ppm
p$P = 2 #muM
p$Fws = (1-p$forg)*p$Fv
p$Fbo = p$Fwo + p$forg*p$Fv
p$Fbcarb = p$Fwcarb + p$Fws

#isotope values
p$eps = 28
p$d13C_volc = -5
p$delC = p$d13C_volc + p$forg*p$eps
p$delwo = p$delC - p$eps

#initial states
p$M0 = c(p$MPss, p$MCss, p$delC, 1)

#sensitivities
p$kws = p$Fws/p$MCss
p$kbc = p$Fbcarb/p$MCss
p$kbo = p$Fbo/p$MPss
p$kbp = p$Fbp/p$MPss
p$kwp = p$Fwp/p$MCss
p$kCO2 = p$pCO2/p$MCss
p$kP = p$P/p$MPss


#state variances
p$sigma_MP = 0.075*p$MPss
p$sigma_MC = 0.075*p$MCss
# p$sigma_delC = 0.5*p$delC
# p$sigma_Fv = 0.01*p$Fv

# p$Q = diag(c(p$sigma_MP^2, p$sigma_MC^2, 0.001, 0.005))
p$Q = diag(c(p$sigma_MP^2, p$sigma_MC^2, 1, 0.05))

#observation errors
p$sigma_d13C = 5*p$delC
# p$sigma_pCO2 = 0.5*p$pCO2
# p$sigma_P = 0.5*p$P

p$R = matrix(p$sigma_d13C^2)

p$W = function(){mvrnorm(n = 1, mu = matrix(0,1,dim(p$Q)[2]), Sigma = p$Q)}
p$V = function(){mvrnorm(n = 1, mu = matrix(0,1,dim(p$R)[2]), Sigma = p$R)}

p$Q_sample = diag(c(p$sigma_MP^2, p$sigma_MC^2, 0, 0.1))
p$R_sample = 0.5*p$delC

#forcing functions for synthetic data
# p$F_fun_d13C = function(t) (1 + 3*sin(1*pi*t/1e6) )
p$F_forcing  = function(t) (1 + 10e-1*sin(-1*pi*t/1e6) )

#ensemble characteristics
p$n_ensemble = 1000
p$n_trajectories = 50




