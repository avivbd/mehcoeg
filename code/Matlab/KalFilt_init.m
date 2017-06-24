function [ p ] = KalFilt_init(  )
%UNTITLED initialization of p used for carbon cycle modeling
%   Detailed explanation goes here

% need to scale

p.t0 = 0;
p.tf = 10e6;
p.n_steps = 100;
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

%sensitivities
p.kws = p.Fws/p.MCss;
p.kbo = p.Fbo/p.MPss;
p.kbp = p.Fbp/p.MPss;
p.kwp = p.Fwp/p.MCss;

%isotope values
p.eps = 28;
p.d13C_volc = -5;
p.delC = p.d13C_volc + forg*p.eps;
p.delwo = p.delC - p.eps;
end

