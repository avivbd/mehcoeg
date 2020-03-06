function EnKF_C_Cycle()
%% initialize

close all 
clear
clc

rng('default')


%load init file
p = KalFilt_init;

%make some synthetic data 
[Mp_out, Mc_out, d13C_out, d13C_noisy, pCO2_noisy, P_noisy] = make_data(p);


%Kalman filter it
[Mp_filt, Mc_filt, d13C_filt, F_forcing_filt] = enkf([d13C_noisy; pCO2_noisy; P_noisy], p);


%plot results
plotting_fun(Mp_out, Mc_out, d13C_out, d13C_noisy, Mp_filt,... 
             Mc_filt, d13C_filt, F_forcing_filt, p)




end

function [Mp_out, Mc_out, delC_out, d13C_noisy, pCO2_noisy, P_noisy] = make_data(p)

% rng(1980)

M0 = [p.MPss; p.MCss; p.delC; p.Fv];

[T,Y] = ode45( @(t,y) odeFun(t,y,p), p.tt, M0); 

Mp_out = Y(:, 1);
Mc_out = Y(:, 2);
delC_out = Y(:, 3);
Fv_out = Y(:,4);
                        
Z = obsFun([Mp_out'; Mc_out'; delC_out'; Fv_out'], p);

d13C_obs = Z(1,:);
pCO2_obs = Z(2,:);
P_obs = Z(3,:);


%add trend to d13C
d13C_noisy = d13C_obs.*p.F_fun_d13C(T') + p.sigma_d13C*randn(size(d13C_obs));
pCO2_noisy = pCO2_obs.*p.F_fun_pCO2(T') + p.sigma_pCO2*randn(size(pCO2_obs));
P_noisy = P_obs.*p.F_fun_P(T') + p.sigma_P*randn(size(P_obs));

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
P = p.kP.*M(1,:);
Z = [delC; pCO2; P];

end


function [Mp_filt, Mc_filt, delC_filt, F_forcing_filt, p] = enkf(y, p)

muM = [p.MPss; p.MCss; p.delC; p.Fv]; %augemented state vector
Z_cov = diag(p.z.^2); %create measurement noise covariance matrix


for k=1:(length(y))

    W = repmat(p.w', 1, p.n_ensemble) .* randn(length(muM), p.n_ensemble);            %create process noise
    V = repmat(p.z, 1, p.n_ensemble) .* randn(size(y,1), p.n_ensemble);     %create measurement noise
    M = repmat(muM, 1, p.n_ensemble);     %state ensemble
    M = M + p.dt*odeFun(p.dt*k, M, p) + W;     %Update states. 
    Y = obsFun(M, p) + V; %make ensemble of predictions. 
    y_obs = repmat(y(:,k),1, p.n_ensemble) + V;  %make ensemble of observations
    
    muM = mean(M, 2);
    muY = mean(Y, 2);
    
    EM = M - repmat(muM, 1, length(M));
    EY = Y - repmat(muY, 1, length(Y));
    
    Pxy=EM*EY'/(length(EY)-1);
    Pyy=EY*EY'/(length(EY)-1) + Z_cov;
    
    K=Pxy*inv(Pyy);
    
    M = M + K*(y_obs - Y);
    
    
    MM(:,k) = muM; 

end

Mp_filt = MM(1,:);
Mc_filt = MM(2,:);
delC_filt = MM(3,:);
F_forcing_filt = MM(4,:);




end


function [] = plotting_fun(Mp_out, Mc_out, d13C_out, d13C_noisy, Mp_filt,... 
             Mc_filt, d13C_filt, F_forcing_filt, p)

figure('Position',[360 7 882 689])
axf = subplot(221);
plot(p.tt,p.Fv*ones(size(p.tt)),'k')
hold on
plot(p.tt, F_forcing_filt', 'r', 'Parent',axf) 
xlabel('Time [My]')
ylabel('F_{volc}')
legend(axf, {'model input', 'reconstructed estimate'})


ax1 = subplot(222);
plot(p.tt, Mp_out,'k')
hold on
plot(p.tt, Mp_filt', 'b','Parent',ax1)
ylabel('M_P')
legend(ax1, {'model input', 'reconstructed estimate'})

ax2 = subplot(223);
plot(p.tt, Mc_out,'k')
ylabel('M_C')
hold on
plot(p.tt, Mc_filt','g','Parent',ax2)
legend(ax2, {'model input', 'reconstructed estimate'})

ax3 = subplot(224);
plot(p.tt, d13C_out,'k')
hold on
% plot(p.tt, d13C_forced','b','Parent',ax3)
plot(p.tt, d13C_noisy','ob','Parent',ax3)
plot(p.tt, d13C_filt','r','Parent',ax3)
ylabel(['\delta^{13}C ', char(8240)])
xlabel('Time [yr]')
legend(ax3,{'true curve',...
    'noisy observations','filtered observations'})


end

