function forward_model()
clear
clc
clf

p = param_vals();
p.tspan = linspace(-253, -246, 100);

A1 = 0.75;
A2 = 0.25;
omega1 = 0.5;
omega1 = 1.3;
phi1 = pi/2;
phi2 = pi/3;
p.F_C_volc_t = p.F_C_volc_0*( 1 + A1*sin(omega1*p.tspan + phi1) + ...
                                  A2*sin(omega1*p.tspan + phi2)); 

M0 = [p.M_C , p.delC, p.M_ALK, p.M_Sr, p.R_Sr];
[T, Y] = ode45( @(t, x) ode_eqs(x, t, p), p.tspan, M0);

plot_ode(T,Y, p)

end


function p = param_vals()

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

function [dm] = ode_eqs(M, t, p)
% [p.M_C , p.M_Sr, p.delC, p.R_Sr]
M_C = M(1);
del_C = M(2);
M_Alk = M(3);
M_Sr = M(4);
R_Sr = M(5);

[pCO2, CO3] = c_sys(M_Alk, M_C, p);
RCO2 = pCO2/p.pCO2_i;

[RWsil, RWcarb] = weathering(RCO2);
omega = CO3/p.CO3_i;
% F_C_bcarb = omega/omega_i add a log?
% F_b_org = Fpb*CP

F_C_volc_t = interp1(p.tspan, p.F_C_volc_t, t);

dMC = F_C_volc_t + p.F_C_worg + p.F_C_wcarb*RWcarb ...
      - p.F_C_bcarb*omega - p.F_C_borg;

d_delC_MC = F_C_volc_t*(p.delC_volc - del_C) + ...
            p.F_C_wcarb*RWcarb*(p.delC_wcarb - del_C) + ...
            p.F_C_worg*(p.delC_worg - del_C) - ...
            p.F_C_borg*(p.epsilon);

dMALK = 2*p.F_Ca_wsil*RWsil + 2*p.F_C_wcarb*RWcarb - 2*p.F_C_bcarb;
% 
% dMP = F_P_w - F_P_b;
% 
dMSr = p.F_Sr_hydro_in + p.F_Sr_w - p.F_Sr_b - p.F_Sr_hydro_out;
% 
dRSr = p.F_Sr_hydro_in*p.R_Sr_hydro + p.F_Sr_w*p.R_Sr_riv -...
        p.F_Sr_hydro_out*p.R_Sr_hydro - p.F_Sr_b*p.R_Sr_riv;

dm = [dMC; d_delC_MC ; dMALK; dMSr; dRSr];

end

function plot_ode(T,Y, p)

ylab = {'MC', 'delC', 'MALK','M_Sr', 'RSr'};

for i = 1:length(ylab)
    subplot(3,2,i)
    plot(T, Y(:, i))
    ylabel(ylab(i))
    grid on
    
end

subplot(3, 2, i+1)
plot(p.tspan, p.F_C_volc_t)
ylabel('Fvolc')

xlabel('time [Ma]')
grid on


end


function [RWsil, RWcarb] = weathering(RCO2)

G = 3.3;
Z_s = 0.09;
Z_c = 0.07;

RWsil = (RCO2^(G*Z_s))*(1 + (G*Z_c)*log(RCO2))^0.65;
RWcarb = (RCO2^(G*Z_s))*(1 + (G*Z_s)*log(RCO2));

end

function [pCO2, CO3] = c_sys(M_ALK, M_C, p)

DIC = M_C*1e21/p.Voc;
ALK = M_ALK*1e21/p.Voc;
alpha = 3.4e-5; % mM/ppm
k1 = 1e-6;
k2 = 1e-9;
k = k1/k2;
bar = (4*ALK + DIC*k - ALK*k)^2 + 4*(k-4)*ALK^2;
foo = ALK*k - DIC*k - 4*ALK + sqrt(bar);
H2CO3 = ((DIC - ALK) + foo/(2*(k-4)));
pCO2 = H2CO3/alpha;
HCO3 = H2CO3*k1;
CO3 = HCO3*k2;
end