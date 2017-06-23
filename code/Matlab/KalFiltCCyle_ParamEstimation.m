function KalFiltCCyle_ParamEstimation()
%% initialize

close all 
clear
clc

%load init file
p = KalFilt_init;

%% make some synthetic data 
[t_obs, d13C_obs] = make_data(p);


%% Kalman filter it
[t_filt, Mp_filt, Mc_filt] = kalf(t_obs, d13C_obs, p);

% plot(

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
d13C_filt = obsfun([Mp_filt', Mc_filt'],p);


plotting_fun(p)



end


function dy = odeFun(~,y,p)

    Mp = y(1);
    Mc = y(2);
    dMp = p.kwp*Mc - p.kbp*Mp;
    dMc = p.Fv + p.Fwo - p.kws*Mc - p.kbo*Mp;
    dy = [dMp;dMc];

end

function [J] = obsfun_deriv_numeric(M, obsfun, params) 
    Mp = M(:,1);
    Mc = M(:,2);
    J = jacobianest( @(M) obsfun(M, params), [Mp, Mc]);                    
end

function delta = obsfun(M,p) 
    Mp = M(:,1);
    Mc = M(:,2);
    delta = p.d13C_volc + p.kbo * Mp./(p.kbo*Mp + p.kws*Mc + p.Fwcarb)*p.eps;
end

function [] = plotting_fun(p)

figure('Position',[360 7 882 689])
axf = subplot(221);
plot(p.tt,p.F_forcing(tt),'k')
hold on
xlabel('Time [My]')
ylabel('F_{volc}')

ax1 = subplot(222);
plot(p.T, p.Mpout,'k')
hold on
ylabel('M_P')

ax2 = subplot(223);
plot(p.T, p.Mcout,'k')
ylabel('M_C')
hold on

ax3 = subplot(224);
plot(p.T, p.d13C_out,'k')
hold on
ylabel(['\delta^{13}C ', char(8240)])
xlabel('Time [yr]')

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

function [T, z] = make_data(p)

% rng(1980)

F_fun = @(t) (-1e-7*t+1 + 1*sin(2*pi*t/10e6) );

p.tt = linspace(p.t0,p.tf,30);

p.M0 = [p.MPss; p.MCss];

[T,Y] = ode45(@(t,y) odeFun(t,y,p),[p.t0, p.tf], p.M0); 

Mp_out = Y(:,1);
Mc_out = Y(:,2);
                        
d13C_out = obsfun([Mp_out, Mc_out], p);

sigma = 0.2;
zz = d13C_out.*F_fun(T); 
z = zz + sigma*randn(size(d13C_out));
end

function [t_filt, Mp_filt, Mc_filt] = kalf(t_obs, d13C_obs, p)

M = [p.MPss; p.MCss; p.Fv]; %augemented states vector
P = eye(length(M));
F = [-p.kbp, p.kwp, 0; -p.kbo, -p.kws, 1; 0, 0, 0];
Q = diag([0.001*p.MPss, 0.001*p.MCss, 0.01*p.Fv]);
dt = mean(diff(t_obs));
A = eye(3) + dt*F;
R = 1e-13;

MM = zeros(size(M,1), length(d13C_obs));
PP = zeros(size(M,1), size(M,1), length(d13C_obs));

for k=1:(length(d13C_obs)-1)
    U = [0; p.Fwo; 0];
    H = obsfun_deriv_numeric([M(1),M(2)], @obsfun,p);
    H = [H, 0];
    M = A*M + dt*U;
    P = A*P*A' + Q;

    IM = obsfun([M(1),M(2)], p);
    S = (R + H*P*H');
    K = P*H'/S;
    M = M + K * (d13C_obs(k)-IM);
    P = P - K*S*K';
    
    MM(:,k) = M; 
    PP(:,:,k) = P;

end

t_filt = [0 (1:(length(d13C_obs)-1))*dt];
Mp_filt = MM(1,:);
Mc_filt = MM(2,:);


end