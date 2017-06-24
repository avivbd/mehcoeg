function EKF_ParamEstimation()
%% initialize

close all 
clear
clc

%load init file
p = KalFilt_init;

%make some synthetic data 
[Mp_out, Mc_out, d13C_out, d13C_forced, d13C_noisy] = make_data(p);

%Kalman filter it
[MM, PP, AA, Q, U] = kalf(d13C_noisy, p);

%RTS smooth it
[Mp_filt, Mc_filt, d13C_filt, F_forcing_filt] = RTS(MM, PP, AA, Q, U, p);

%plot results
plotting_fun(Mp_out, Mc_out, d13C_out, d13C_forced, d13C_noisy, Mp_filt,... 
             Mc_filt, d13C_filt, F_forcing_filt, p)



end


function dy = odeFun(~,y,p)

    Mp = y(1);
    Mc = y(2);
    delC = y(3);
    
    dMp = p.kwp*Mc - p.kbp*Mp;
    dMc = p.Fv + p.Fwo - p.kws*Mc - p.kbo*Mp;
    ddelta = ( p.Fv*(p.d13C_volc - delC) + p.Fwo*(p.delwo - delC) + ...
               p.Fwcarb*(p.delC - delC) + p.kbo*Mp*p.eps )/Mc;
    
    dy = [dMp; dMc; ddelta];

end

function [J] = obsfun_deriv_numeric(M, obsfun, params) 
    Mp = M(:,1);
    Mc = M(:,2);
    delC = M(:,3);
    J = jacobianest( @(M) obsfun(M, params), [Mp, Mc, delC]);                    
end

function [J] = odefun_deriv_numeric(M, odeFun, params) 
    Mp = M(:,1);
    Mc = M(:,2);
    delC = M(:,3);
    J = jacobianest( @(M) odeFun([], M, params), [Mp, Mc, delC]);                    
end


function delta = obsfun(M,p)
%     Mp = M(:,1);
%     Mc = M(:,2);
    delC = M(:,3);
%     delta = p.d13C_volc + p.kbo * Mp./(p.kbo*Mp + p.kws*Mc + p.Fwcarb)*p.eps;
    delta = delC;
end

function [Mp_out, Mc_out, delC_out, d13C_forced, d13C_noisy] = make_data(p)

% rng(1980)

M0 = [p.MPss; p.MCss; p.delC];

[T,Y] = ode45( @(t,y) odeFun(t,y,p), p.tt, M0); 

Mp_out = Y(:, 1);
Mc_out = Y(:, 2);
delC_out = Y(:, 3);
                        
d13C_obs = obsfun([Mp_out, Mc_out, delC_out], p);

%add trend to d13C
F_fun = @(t) (-1e-7*t+1 + 1*sin(2*pi*t/10e6) );
sigma = 0.2;
d13C_forced = d13C_obs.*F_fun(T); 
d13C_noisy = d13C_forced + sigma*randn(size(d13C_obs));
end

function [MM, PP, AA, Q, U] = kalf(d13C_obs, p)

M = [p.MPss; p.MCss; p.delC; p.Fv]; %augemented state vector
P = eye(length(M));
% F = [-p.kbp, p.kwp, 0; -p.kbo, -p.kws, 1; 0, 0, 0];
F = odefun_deriv_numeric([M(1), M(2), M(3)], @odeFun, p);
F = padmat(F);

Q = diag([0.01*p.MPss, 0.01*p.MCss, 0.001*p.delC, 0.01*p.Fv]);
A = eye(length(M)) + p.dt*F;
R = 1e-13;

MM = zeros(size(M,1), length(d13C_obs));
PP = zeros(size(M,1), size(M,1), length(d13C_obs));
AA = zeros(size(M,1), size(M,1), length(d13C_obs));

for k=1:(length(d13C_obs))
    U = [0; p.Fwo; 0; 0];
    H = obsfun_deriv_numeric([M(1),M(2), M(3)], @obsfun, p);
    H = [H, 0];
    A = odefun_deriv_numeric([M(1), M(2), M(3)], @odeFun, p);
    A = padmat(A);
    M = M + (A*M + U)*p.dt;%A*M + p.dt*U;
    P = A*P*A' + Q;

    IM = obsfun([M(1), M(2), M(3)], p);
    S = (R + H*P*H');
    K = P*H'/S;
    M = M + K * (d13C_obs(k)-IM);
    P = P - K*S*K';
    
    MM(:,k) = M; 
    PP(:,:,k) = P;
    AA(:,:,k) = A;

end

% Mp_filt = MM(1,:);
% Mc_filt = MM(2,:);
% F_forcing_filt = MM(3,:);

end

function [Mp_filt, Mc_filt, d13C_filt, F_forcing_filt]= RTS(M, P, A, Q, U, p)

% M = [p.MPss; p.MCss; p.Fv]; %augemented states vector
D = zeros(size(M,1),size(M,1),size(M,2));
% M = MM;
% P = PP;
Q = repmat(Q,[1 1 size(M,2)]);

% PP2 = zeros(size(M,1),length(Y));

for k=(size(M,2)-1):-1:1
    P_pred   = A(:,:,k) * P(:,:,k) * A(:,:,k)' + Q(:,:,k);
    D(:,:,k) = P(:,:,k) * A(:,:,k)' / P_pred;
    M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - (A(:,:,k) * M(:,k) + p.dt*U));
    P(:,:,k) = P(:,:,k) + D(:,:,k) * (P(:,:,k+1) - P_pred) * D(:,:,k)';
    
%     PP2(:,k) = diag(chol(P(:,:,k)));
    
end

% sigma_Mp_filt = PP2(1,:);
% sigma_Mc_filt = PP(2,:);


Mp_filt = M(1,:);
Mc_filt = M(2,:);
delC_filt = M(3,:);
F_forcing_filt = M(4,:);
d13C_filt = obsfun([Mp_filt', Mc_filt', delC_filt'], p);

end

function [] = plotting_fun(Mp_out, Mc_out, d13C_out, d13C_forced, d13C_noisy, Mp_filt,... 
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
plot(p.tt, d13C_forced','b','Parent',ax3)
plot(p.tt, d13C_noisy','ob','Parent',ax3)
plot(p.tt, d13C_filt','r','Parent',ax3)
ylabel(['\delta^{13}C ', char(8240)])
xlabel('Time [yr]')
legend(ax3,{'forward model output','true curve',...
    'noisy observations','filtered observations'})


end

function F = padmat(F)

F = [F; zeros(1, size(F,2))];
F = [F , [0; 1; 0; 0]];

end
