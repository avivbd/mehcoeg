function [] = KalFilt_C_Cycle

%todo: estimate changing sensitivity

close all 
clear
clc

% rng(11051980)


params = KalFilt_init;
KalFiltCCyle(params);

% KalFiltCCyle_ParamEstimation()
% UKF_CC_ParamEst()

end


function KalFiltCCyle()
%%
%reservoirs
MCss = 3.8e18;
MPss = 2e15;
%fluxes
Fbo = 10e12;
Fwo = Fbo;
Fwcarb = 34e12;
Fv = 6e12;
Fws = Fv;
Fbcarb = Fwcarb+Fws;
Fwp = 3.6e10;
Fbp = Fwp;
%sensitivities
kws = Fws/MCss;
kbo = Fbo/MPss;
kbp = Fbp/MPss;
kwp = Fwp/MCss;
%isotope values
eps = 28;
d13C_volc = -5;
delwo = -25;


d13C_carb = d13C_volc + kbo*MPss/(kbo*MPss + kws*MCss + Fwcarb)*eps;



%%
T = 10e6;
a = -0.25;
F_forcing = @(t) Fv*(1+a*(sin(2*pi*t/T).*(1 - cos(pi*t/(0.5*T)))));

%%
tt = linspace(0,10e6,30);
figure('Position',[360 7 882 689])
axf = subplot(221);
plot(tt,F_forcing(tt))
xlabel('Time [My]')
ylabel('F_{volc}')

%%

M0 = [MPss;MCss];
TT = linspace(0,10e6,100);
params = struct('F_forcing',F_forcing,'kwp', kwp,...
    'kbp', kbp, 'Fwo', Fwo, 'kws', kws,'kbo', kbo);

[T,Y] = ode45(@(t,y) odeFun(t,y,params), TT, M0); 
Mpout = Y(:,1);
Mcout = Y(:,2);

%%

ax1 = subplot(222);
plot(T, Mpout,'k')
hold on
ylabel('M_P')

ax2 = subplot(223);
plot(T, Mcout,'k')
ylabel('M_C')
hold on
%%
for i = 1:10
    Q = diag([1e-4*MPss,1e-4*MCss]);
    options = sdeset('DiagonalNoise','no');

    YOUT = sde_euler(@(t,y) odeFun(t,y,params), Q ,TT, M0 + randn(2,1),options); 
    Mp_stoch = YOUT(:,1);
    Mc_stoch = YOUT(:,2);

    plot(T, Mp_stoch, 'LineStyle', '-','Color',[0.5 0.5 0.5],'Parent',ax1)
    plot(T, Mc_stoch, 'LineStyle','-','Color',[0.5 0.5 0.5],'Parent',ax2)
end

%%
                        
d13C_out = obsfun(Mpout, Mcout);
ax3 = subplot(224);
plot(T, d13C_out,'k')
hold on
ylabel(['\delta^{13}C ', char(8240)])
xlabel('Time [yr]')

%%
%make observations by adding noise to the model d13C output

sigma = 0.2;
% z = d13C_out + sigma*randn(size(d13C_out));
z = d13C_out + sigma*randn(size(d13C_out));
plot(T, z,'o','Parent',ax3)

%%

M = M0;
P = eye(length(M));
F = [-kbp kwp; -kbo, -kws];
dt = TT(2) - TT(1);
A = eye(2) + dt*F;
R = 1e-15;
% [A,Q] = lti_disc(F,[],Qc,dt);

MM2 = zeros(size(M,1),size(Y,2));
% PP2 = zeros(size(M,1),size(M,1),size(Y,2));
PP2 = zeros(size(M,1),size(Y,2));

for k=1:length(z)
    U = [0;F_forcing(TT(k))+Fwo];
    H = obsfun_deriv(M(1),M(2));
    
    M = A*M + dt*U;
    P = A*P*A' + Q;

%     IM = H*M;
    IM = obsfun(M(1),M(2));
    S = (R + H*P*H');
    K = P*H'/S;
    M = M + K * (z(k)-IM);
    P = P - K*S*K';
    
    MM2(:,k) = M; 
%     PP2(:,:,k) = P;
    PP2(:,k) = diag(chol(P));
end

Mp_filt = MM2(1,:);
Mc_filt = MM2(2,:);
d13C_filt = obsfun(Mp_filt, Mc_filt);

sigma_Mp_filt = PP2(1,:);
sigma_Mc_filt = PP2(2,:);
%%

% figure
plot(TT,Mp_filt, 'b','Parent',ax1)
% hold on 
% plot(TT, Mp_filt-sigma_Mp_filt,'--r',...
%      TT, Mp_filt+sigma_Mp_filt,'--r')

% subplot 312
plot(TT,Mc_filt,'g','Parent',ax2)
% hold on 
% plot(TT, Mc_filt-sigma_Mc_filt,'--r',...
%      TT, Mc_filt+sigma_Mc_filt,'--r')


% subplot 313
plot(TT, d13C_filt,'r','Parent',ax3)


end

function dy = odeFun(t,y,p)

    Mp = y(1);
    Mc = y(2);
    dMp = p.kwp*Mc - p.kbp*Mp;
    dMc = p.F_forcing(t) + p.Fwo - p.kws*Mc - p.kbo*Mp;
    dy = [dMp;dMc];

end

function delta = obsfun(M,p) 
    Mp = M(:,1);
    Mc = M(:,2);
    delta = p.d13C_volc + p.kbo * Mp./(p.kbo*Mp + p.kws*Mc + p.Fwcarb)*p.eps;
end

function [J] = obsfun_deriv(M,p)
    Mp = M(:,1);
    Mc = M(:,2);
    
    J = [ (p.eps*p.kbo)/(p.Fwcarb + Mp*p.kbo + Mc*p.kws) - ...
      (Mp*p.eps*p.kbo^2)/(p.Fwcarb + Mp*p.kbo + Mc*p.kws)^2, ...
      -(Mp*p.eps*p.kbo*p.kws)/(p.Fwcarb + Mp*p.kbo + Mc*p.kws)^2];
end                        




