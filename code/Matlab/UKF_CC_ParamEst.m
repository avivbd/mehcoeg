function UKF_CC_ParamEst()



%%
tstart = 0;
tend = 15e6;
% T = 10e6;
% a = -0.25;
% F_forcing = @(t) Fv*(1+a*(sin(2*pi*t/T).*(1 - cos(pi*t/(0.5*T)))));
F_forcing = @(t) Fv*ones(length(t));
tt = linspace(tstart,tend,30);

figure('Position',[360 7 882 689])
axf = subplot(221);
plot(tt,F_forcing(tt),'k')
hold on
xlabel('Time [My]')
ylabel('F_{volc}')

%%

M0 = [MPss;MCss];
TT = linspace(tstart,tend,100);
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

obsfun = @(Mp, Mc) d13C_volc + kbo*Mp./(kbo*Mp + kws*Mc + Fwcarb)*eps;


d13C_out = obsfun([Mpout, Mcout]);
ax3 = subplot(224);
plot(T, d13C_out,'k')
hold on
ylabel(['\delta^{13}C ', char(8240)])
xlabel('Time [yr]')


%%

rng(1980)
sigma = 0.2;
% z = d13C_out + sigma*randn(size(d13C_out));
zz = d13C_out.*(-1e-7*T+1 + 1*sin(2*pi*T/10e6) ); 
z = zz + sigma*randn(size(d13C_out));
plot(T, zz, '-','Parent',ax3)
plot(T, z,'o','Parent',ax3)


%% UKF that shit

  [WM,WC,c] = ut_weights(size(M,1),alpha,beta,kappa);
   X = ut_sigmas(M,P,c);
   w = {WM,WC,c};

  tr_param = {alpha beta kappa mat};
  [M,P,D] = ut_transform(M,P,f,f_param,tr_param);
  P = P + Q;




end