function [ p ] = KalFilt_init(  )
%UNTITLED initialization of p used for carbon cycle modeling
%   Detailed explanation goes here

% need to scale

p.t0 = 0;
p.tf = 1e6;
p.n_steps = 150;
p.dt = (p.tf - p.t0)/p.n_steps;
p.tt = p.t0:p.dt:p.tf;

%reservoirs
p.MCss = 3.8e18;
p.MPss = 2e15;

%fluxes
% total_input = 50e12;
% fworg = 0.2;
forg = 0.2;
p.Fv = 5e12; 
p.Fwo = 9e12;
p.Fwcarb = 36e12;

p.Fws = (1-forg)*p.Fv;
p.Fbo = p.Fwo + forg*p.Fv;
p.Fbcarb = p.Fwcarb + p.Fws;


p.Fwp = 3.6e10;
p.Fbp = p.Fwp;

p.pCO2 = 540; %ppm
p.P = 2; %muM

%sensitivities
p.kws = p.Fws/p.MCss;
p.kbo = p.Fbo/p.MPss;
p.kbp = p.Fbp/p.MPss;
p.kwp = p.Fwp/p.MCss;
p.kCO2 = p.pCO2/p.MCss;
p.kP = p.P/p.MPss;

%isotope values
p.eps = 28;
p.d13C_volc = -5;
p.delC = p.d13C_volc + forg*p.eps;
p.delwo = p.delC - p.eps;

%states
p.M = [p.MPss; p.MCss; p.delC];

%states variances
p.sigma_MP = 0.5*p.MPss;
p.sigma_MC = 0.05*p.MCss;
p.sigma_delC = 0.5*p.delC;
% p.sigma_Fv = 0.01*p.Fv;

%observation errors
p.sigma_d13C = 0.1*p.delC;
p.sigma_pCO2 = 0.5*p.pCO2;
p.sigma_P = 0.5*p.P;

%forcing functions
p.F_fun_d13C = @(t) (1 + 3*sin(1*pi*t/1e6) );
p.F_fun_pCO2 = @(t) (1 - 1e-1*t/1e6 );
p.F_fun_P  = @(t) (1 + 1e-1*sin(-2*pi*t/1e6) );

%ensemble characteristics
p.n_ensemble = 100;
% p.w = [p.sigma_MP; p.sigma_MC; p.sigma_delC; p.sigma_Fv];  %process noise standard deviation
p.w = [p.sigma_MP; p.sigma_MC; p.sigma_delC];  %process noise standard deviation
% p.z = [p.sigma_d13C; p.sigma_pCO2; p.sigma_P];  %measurement noise standard deviation
% p.z = [p.sigma_d13C; p.sigma_pCO2];
p.z = [p.sigma_d13C];

end

