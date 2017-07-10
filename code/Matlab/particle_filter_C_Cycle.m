function particle_filter_C_Cycle()
%% initialize

close all 
clear
clc

rng('default')


%load init file
p = KalFilt_init;

%make some synthetic data 
[Mdl, Data] = make_data(p);

% data = [d13C_noisy; pCO2_noisy; P_noisy];
data = [Data.d13C_noisy; Data.pCO2_noisy];

p.odeFun = @odeFun;
p.obsFun = @obsFun;


% filter it
[S, E] = part_filt(data, p);


%plot results
plotting_fun(Mdl, Data, E, S, p)




end

function [Mdl, Data] = make_data(p)

% rng(1980)

M0 = [p.MPss; p.MCss; p.delC; p.Fv];

[T,Y] = ode45( @(t,y) odeFun(t,y,p), p.tt, M0); 

Mdl.Mp_out = Y(:, 1);
Mdl.Mc_out = Y(:, 2);
Mdl.delC_out = Y(:, 3);
Mdl.Fv_out = Y(:,4);
                        
Z = obsFun([Mdl.Mp_out'; Mdl.Mc_out'; Mdl.delC_out'; Mdl.Fv_out'], p);

d13C_obs = Z(1,:);
pCO2_obs = Z(2,:);
% P_obs = Z(3,:);


%add trend to d13C
Data.d13C_noisy = d13C_obs.*p.F_fun_d13C(T') + p.sigma_d13C*randn(size(d13C_obs));
Data.pCO2_noisy = pCO2_obs.*p.F_fun_pCO2(T') + p.sigma_pCO2*randn(size(pCO2_obs));
% P_noisy = P_obs.*p.F_fun_P(T') + p.sigma_P*randn(size(P_obs));

end

function dy = odeFun(t,y,p)

    Mp = y(1,:);
    Mc = y(2,:);
    delC = y(3,:);
    Fv = y(4,:);
    
    dMp = p.kwp*Mc - p.kbp*Mp;
    
    dMc = Fv + p.Fwo - p.kws*Mc - p.kbo*Mp;
    
    ddelta = ( p.Fv*(p.d13C_volc - delC) + p.Fwo*(p.delwo - delC) + ...
               p.Fwcarb*(p.delC - delC) + p.kbo*Mp*p.eps )./Mc;
    
    dy = [dMp; dMc; ddelta; zeros(size(ddelta))];

end

function Z = obsFun(M,p)
    
pCO2 = p.kCO2.*M(2,:);
delC = M(3,:);
% P = p.kP.*M(1,:);
% Z = [delC; pCO2; P];
Z = [delC; pCO2];

end

function [] = plotting_fun(Mdl, Data, E, S, p)

figure('Position',[360 7 882 689])
axf = subplot(221);
hold on
plot(p.tt, Mdl.Fv_out,'k')
plot(p.tt, S.F_forcing_filt', 'b') 
% axis manual
he = plot(p.tt, E.F_forcing, '.', 'color', [0,0,0]+0.5);
xlabel('Time [My]')
ylabel('F_{volc}')
legend(axf, {'model input', 'filtered output', 'ensemble'})
uistack(he, 'down',3)



ax1 = subplot(222);
hold on
plot(p.tt, p.kP.*Mdl.Mp_out,'k')
plot(p.tt, p.kP.*S.Mp_filt', 'b')
% hd = plot(p.tt, Data.P_noisy,'ro' , 'MarkerSize', 3);
axis manual
he = plot(p.tt, p.kP.*E.Mp, '.', 'color', [0,0,0]+0.5);
ylabel('[P]')
legend(ax1, {'model input'; 'filtered output'; 'data'; 'ensemble'})
uistack(he, 'down',3)
% uistack(hd, 'down',2)


ax2 = subplot(223);
hold on
plot(p.tt, p.kCO2.*Mdl.Mc_out,'k')
plot(p.tt, p.kCO2.*S.Mc_filt','b')
hd = plot(p.tt, Data.pCO2_noisy, 'ro', 'MarkerSize', 3);
axis manual
he = plot(p.tt, p.kCO2.*E.Mc, '.', 'color', [0,0,0]+0.5);
ylabel('pCO_2')
legend(ax2, {'model input'; 'filtered output'; 'data'; 'ensemble'})
uistack(he, 'down',3)
uistack(hd, 'down',2)

ax3 = subplot(224);
hold on
plot(p.tt, Mdl.delC_out,'k')
plot(p.tt, S.delC_filt','b')
hd = plot(p.tt, Data.d13C_noisy','ro', 'MarkerSize', 3);
axis manual
he = plot(p.tt', E.delC, '.', 'color', [0,0,0]+0.5);
ylabel(['\delta^{13}C ', char(8240)])
xlabel('Time [yr]')
legend(ax3, {'model input'; 'filtered output'; 'data'; 'ensemble'})
uistack(he, 'down',3)
uistack(hd, 'down',2)
end

