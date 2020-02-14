function [ p ] = params( )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

p = struct();

% All values from Kump and Arthur 1997
p.R_Sr_hydro = 0.7035;
p.R_Sr_riv = 0.7101;
p.R_Sr_wsil = 0.7164;
p.R_Sr_wcarb = 0.7080;

p.F_Sr_w = 0.022; %10^18 / 1 Ma = Tmol/yr
p.F_Sr_b = p.F_Sr_w;
p.F_Sr_hydro_in = 0.0035;
p.F_Sr_hydro_out = p.F_Sr_hydro_in;


p.F_C_worg = 5;
p.F_C_wcarb = 20;
p.F_C_volc_0 = 10;

p.F_C_bcarb = 30;
p.F_C_borg = 5; 

p.F_Ca_wsil = 10;

p.delC_volc = -5;
p.delC_wcarb = +2;
p.delC_worg = -22;
p.delC_Fbcarb = 0;
p.delC_Fbo = -22;
p.epsilon = -22;

p.M_Sr = 0.12; % 10^18 Mol
p.M_C = 7.42; % 10^18 Mol
p.R_Sr = 0.7092;
p.delC = 0; 

p.M_ALK = 8.06;
p.Voc = 1.32e21; %[L]
[p.pCO2_i, p.CO3_i] = c_sys(p.M_ALK, p.M_C, p);


end

