function [] = KalFilt_C_Cycle

%todo: estimate changing sensitivity

close all 
clear
clc

params = KalFilt_init;
KalFiltCCyle_ParamEstimation(params)

end




function KalFiltCCyle_ParamEstimation(params)
%%
params.F_forcing = @(t) params.Fv*ones(size(t)); %flat forcing
tt = linspace(0,10e6,30);

figure('Position',[360 7 882 689])
axf = subplot(221);
plot(tt,params.F_forcing(tt),'k')
hold on
xlabel('Time [My]')
ylabel('F_{volc}')

%%

M0 = [params.MPss; params.MCss];
TT = linspace(0,10e6,100);
% params = struct('F_forcing',F_forcing,'kwp', kwp,...
%     'kbp', kbp, 'Fwo', Fwo, 'kws', kws,'kbo', kbo);

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


% obsfun = @(Mp, Mc) d13C_volc + kbo*Mp./(kbo*Mp + kws*Mc + Fwcarb)*eps;

% %analytical gradient of the observation function
% obsfun_deriv = @(Mp, Mc) [ (eps*kbo)/(Fwcarb + Mp*kbo + Mc*kws) - ...
%                            (Mp*eps*kbo^2)/(Fwcarb + Mp*kbo + Mc*kws)^2, ...
%                             -(Mp*eps*kbo*kws)/(Fwcarb + Mp*kbo + Mc*kws)^2];

% obsfun_deriv_numeric = @(Mp, Mc) jacobianest(obsfun(Mp,Mc),[Mp; Mc]);                    
                        
d13C_out = obsfun([Mpout, Mcout], params);
ax3 = subplot(224);
plot(T, d13C_out,'k')
hold on
ylabel(['\delta^{13}C ', char(8240)])
xlabel('Time [yr]')


%% 
%make observations by adding noise to the model d13C output
rng(1980)
sigma = 0.2;
% z = d13C_out + sigma*randn(size(d13C_out));
zz = d13C_out.*(-1e-7*T+1 + 1*sin(2*pi*T/10e6) ); 
z = zz + sigma*randn(size(d13C_out));
plot(T, zz, '-','Parent',ax3)
plot(T, z,'o','Parent',ax3)



%% Kalman filter it

M = [M0;6e12]; %add the forcing function to the filter
P = eye(length(M));
F = [-params.kbp, params.kwp, 0; -params.kbo, -params.kws, 1; 0, 0, 0];
Q = diag([0.001*params.MPss, 0.001*params.MCss, 0.01*params.Fv]);
dt = TT(2) - TT(1);
A = eye(3) + dt*F;
R = 1e-13;
% [A,Q] = lti_disc(F,[],Qc,dt);

% MM2 = zeros(length(M),length(Y));


MM = zeros(size(M,1), size(Y,2));
PP = zeros(size(M,1), size(M,1), size(Y,2));

for k=1:length(z)
    U = [0;params.Fwo;0];
%     H = [obsfun_deriv([M(1),M(2)],params),0];
    H = obsfun_deriv_numeric([M(1),M(2)], @obsfun,params);
    H = [H, 0];
    M = A*M + dt*U;
    P = A*P*A' + Q;

%     IM = H*M;
    IM = obsfun([M(1),M(2)], params);
    S = (R + H*P*H');
    K = P*H'/S;
    M = M + K * (z(k)-IM);
    P = P - K*S*K';
    
    MM(:,k) = M; 
    PP(:,:,k) = P;

end


Mp_filt = MM(1,:);
Mc_filt = MM(2,:);
% F_forcing_filt = MM(3,:);
% d13C_filt = obsfun(Mp_filt, Mc_filt);

%%

% plot(TT,Mp_filt, 'b','Parent',ax1)
% plot(TT,Mc_filt,'g','Parent',ax2)
% plot(TT, d13C_filt,'r','Parent',ax3)
% plot(TT, F_forcing_filt, 'r', 'Parent',axf) 


%% run the RTS smoother
D = zeros(size(M,1),size(M,1),size(M,2));
M = MM;
P = PP;
A = repmat(A,[1 1 size(M,2)]);
Q = repmat(Q,[1 1 size(M,2)]);

PP2 = zeros(size(M,1),length(Y));

for k=(size(M,2)-1):-1:1
    P_pred   = A(:,:,k) * P(:,:,k) * A(:,:,k)' + Q(:,:,k);
    D(:,:,k) = P(:,:,k) * A(:,:,k)' / P_pred;
    M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - (A(:,:,k) * M(:,k) + dt*U));
    P(:,:,k) = P(:,:,k) + D(:,:,k) * (P(:,:,k+1) - P_pred) * D(:,:,k)';
    
    PP2(:,k) = diag(chol(P(:,:,k)));
    
end

sigma_Mp_filt = PP2(1,:);
sigma_Mc_filt = PP(2,:);


Mp_filt = M(1,:);
Mc_filt = M(2,:);
F_forcing_filt = M(3,:);
d13C_filt = obsfun([Mp_filt', Mc_filt'],params);




plot(TT',Mp_filt', 'b','Parent',ax1)
legend(ax1, {'model input', 'reconstructed estimate'})
plot(TT',Mc_filt','g','Parent',ax2)
legend(ax2, {'model input', 'reconstructed estimate'})
plot(TT', d13C_filt','r','Parent',ax3)
legend(ax3,{'forward model output','true curve',...
    'noisy observations','filtered observations'})
plot(TT', F_forcing_filt', 'r', 'Parent',axf) 
legend(axf, {'model input', 'reconstructed estimate'})

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

function [J] = obsfun_deriv_numeric(M, obsfun, params) 
    Mp = M(:,1);
    Mc = M(:,2);
    J = jacobianest( @(M) obsfun(M, params), [Mp, Mc]);                    
end


%%%%%%%%%%%%%%%%%%%%

% KalFiltCCyle();
% UKF_CC_ParamEst()



% function KalFiltCCyle()
% %%
% %reservoirs
% MCss = 3.8e18;
% MPss = 2e15;
% %fluxes
% Fbo = 10e12;
% Fwo = Fbo;
% Fwcarb = 34e12;
% Fv = 6e12;
% Fws = Fv;
% Fbcarb = Fwcarb+Fws;
% Fwp = 3.6e10;
% Fbp = Fwp;
% %sensitivities
% kws = Fws/MCss;
% kbo = Fbo/MPss;
% kbp = Fbp/MPss;
% kwp = Fwp/MCss;
% %isotope values
% eps = 28;
% d13C_volc = -5;
% delwo = -25;
% 
% 
% d13C_carb = d13C_volc + kbo*MPss/(kbo*MPss + kws*MCss + Fwcarb)*eps;
% 
% 
% 
% %%
% T = 10e6;
% a = -0.25;
% F_forcing = @(t) Fv*(1+a*(sin(2*pi*t/T).*(1 - cos(pi*t/(0.5*T)))));
% 
% %%
% tt = linspace(0,10e6,30);
% figure('Position',[360 7 882 689])
% axf = subplot(221);
% plot(tt,F_forcing(tt))
% xlabel('Time [My]')
% ylabel('F_{volc}')
% 
% %%
% 
% M0 = [MPss;MCss];
% TT = linspace(0,10e6,100);
% params = struct('F_forcing',F_forcing,'kwp', kwp,...
%     'kbp', kbp, 'Fwo', Fwo, 'kws', kws,'kbo', kbo);
% 
% [T,Y] = ode45(@(t,y) odeFun(t,y,params), TT, M0); 
% Mpout = Y(:,1);
% Mcout = Y(:,2);
% 
% %%
% 
% ax1 = subplot(222);
% plot(T, Mpout,'k')
% hold on
% ylabel('M_P')
% 
% ax2 = subplot(223);
% plot(T, Mcout,'k')
% ylabel('M_C')
% hold on
% %%
% for i = 1:10
%     Q = diag([1e-4*MPss,1e-4*MCss]);
%     options = sdeset('DiagonalNoise','no');
% 
%     YOUT = sde_euler(@(t,y) odeFun(t,y,params), Q ,TT, M0 + randn(2,1),options); 
%     Mp_stoch = YOUT(:,1);
%     Mc_stoch = YOUT(:,2);
% 
%     plot(T, Mp_stoch, 'LineStyle', '-','Color',[0.5 0.5 0.5],'Parent',ax1)
%     plot(T, Mc_stoch, 'LineStyle','-','Color',[0.5 0.5 0.5],'Parent',ax2)
% end
% 
% %%
%                         
% d13C_out = obsfun(Mpout, Mcout);
% ax3 = subplot(224);
% plot(T, d13C_out,'k')
% hold on
% ylabel(['\delta^{13}C ', char(8240)])
% xlabel('Time [yr]')
% 
% %%
% %make observations by adding noise to the model d13C output
% 
% sigma = 0.2;
% % z = d13C_out + sigma*randn(size(d13C_out));
% z = d13C_out + sigma*randn(size(d13C_out));
% plot(T, z,'o','Parent',ax3)
% 
% %%
% 
% M = M0;
% P = eye(length(M));
% F = [-kbp kwp; -kbo, -kws];
% dt = TT(2) - TT(1);
% A = eye(2) + dt*F;
% R = 1e-15;
% % [A,Q] = lti_disc(F,[],Qc,dt);
% 
% MM2 = zeros(size(M,1),size(Y,2));
% % PP2 = zeros(size(M,1),size(M,1),size(Y,2));
% PP2 = zeros(size(M,1),size(Y,2));
% 
% for k=1:length(z)
%     U = [0;F_forcing(TT(k))+Fwo];
%     H = obsfun_deriv(M(1),M(2));
%     
%     M = A*M + dt*U;
%     P = A*P*A' + Q;
% 
% %     IM = H*M;
%     IM = obsfun(M(1),M(2));
%     S = (R + H*P*H');
%     K = P*H'/S;
%     M = M + K * (z(k)-IM);
%     P = P - K*S*K';
%     
%     MM2(:,k) = M; 
% %     PP2(:,:,k) = P;
%     PP2(:,k) = diag(chol(P));
% end
% 
% Mp_filt = MM2(1,:);
% Mc_filt = MM2(2,:);
% d13C_filt = obsfun(Mp_filt, Mc_filt);
% 
% sigma_Mp_filt = PP2(1,:);
% sigma_Mc_filt = PP2(2,:);
% %%
% 
% % figure
% plot(TT,Mp_filt, 'b','Parent',ax1)
% % hold on 
% % plot(TT, Mp_filt-sigma_Mp_filt,'--r',...
% %      TT, Mp_filt+sigma_Mp_filt,'--r')
% 
% % subplot 312
% plot(TT,Mc_filt,'g','Parent',ax2)
% % hold on 
% % plot(TT, Mc_filt-sigma_Mc_filt,'--r',...
% %      TT, Mc_filt+sigma_Mc_filt,'--r')
% 
% 
% % subplot 313
% plot(TT, d13C_filt,'r','Parent',ax3)
% 
% 
% end



% function [J] = obsfun_deriv(M,p)
%     Mp = M(:,1);
%     Mc = M(:,2);
%     
%     J = [ (p.eps*p.kbo)/(p.Fwcarb + Mp*p.kbo + Mc*p.kws) - ...
%       (Mp*p.eps*p.kbo^2)/(p.Fwcarb + Mp*p.kbo + Mc*p.kws)^2, ...
%       -(Mp*p.eps*p.kbo*p.kws)/(p.Fwcarb + Mp*p.kbo + Mc*p.kws)^2];
% end                        
