function [ params ] = KalFilt_init(  )
%UNTITLED initialization of params used for carbon cycle modeling
%   Detailed explanation goes here

% need to scale

params.t0 = 0;
params.tf = 10e6;
params.n_steps = 100;
params.dt = (params.tf - params.t0)/params.n_steps;
params.tt = params.t0:params.dt:params.tf;

%reservoirs
params.MCss = 3.8e18;
params.MPss = 2e15;
%fluxes
params.Fbo = 10e12;
params.Fwo = params.Fbo;
params.Fwcarb = 34e12;
params.Fv = 6e12;
params.Fws = params.Fv;
params.Fbcarb = params.Fwcarb + params.Fws;
params.Fwp = 3.6e10;
params.Fbp = params.Fwp;
%sensitivities
params.kws = params.Fws/params.MCss;
params.kbo = params.Fbo/params.MPss;
params.kbp = params.Fbp/params.MPss;
params.kwp = params.Fwp/params.MCss;
%isotope values
params.eps = 28;
params.d13C_volc = -5;
params.delwo = -25;
end

