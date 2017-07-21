function [  ] = KalFilt(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

clc
close all
clear

plotdata = 'off'; %on/off

plothistdiffdata = 'off';

simsyst = 'yes';
%%
[t_PT, delC_PT] = plotfun(plotdata);

histdiffdata(plothistdiffdata,t_PT,delC_PT)




%model values
Fin = 50e12;
MP_0 = 2e15;
MC_0 = 3600e15;
Fwsil_0 = 4e12; 
Fborg_0 = 14e12;
Fwp_0 = 36e9; 
Fbp_0 = 36e9;
kwsil = Fwsil_0/MC_0;    
kborg = Fborg_0/MP_0;
kwp = Fwp_0/MC_0;   
kbp = Fbp_0/MP_0;  
delC_in = -5;
epsilon = 25;
% delC_0 = 0;
%%
% test values
% A = [1.1269   -0.4940    0.1129
%      1.0000         0         0
%           0    1.0000         0];
% 
% B = [-0.3832
%       0.5919
%       0.5191];
% 
%  C = [1 0 0]; 

%% simulate the system 

x_0 = [0;0;0];% start from steady state

A = @jacfun;
B = [0;1;0];
C = [0 0 1];
D = 0;
G = [0;1;0];
H = 1;

% Q = 0;
Q = (MC_0/100000)^2;
R = 0.05;

randn('seed',0)

switch simsyst
    case 'yes'
        tspan = [0 1e6];
        tinterp = 0:1e3:1e6;
        w_interp = sqrt(Q)*randn(length(tinterp),1);
        u_interp = 1e13*sin(2*2*pi*tinterp/1e6);
%         u_interp = zeros(size(tinterp));
        
        options = odeset('OutputFcn',@odeplot,'OutputSel', 3);
        
        [t,x] = ode15s(@(t,x) odefun(t,x,A,B,B,tinterp,w_interp,u_interp),...
                tspan,x_0,options);


        v = sqrt(R)*randn(length(t),1);
        
        for i = 1:length(x);
            yv(i,:) = C*(x(i,:))' + v(i);
        end
        
        figure
        subplot(211); 
        plotyy(t,x(:,1),t,x(:,2)); 
        subplot(212); 
        plotyy(t,x(:,3),t,yv)
        
end
% n = length(t);

% sys = ss(A,B,C,D);
% Q - covariance of process noise
% R - covariance of measurement noise

% R = 0;
% w = process noise
% v = measurement noise





% [y,t,x] = lsim(sys,w,t,x_0);

%%
% Q = 1; 
% R = 1;
% 
% t = (1:length(t_PT))';%[0:100]';
% u = zeros(size(t));%sin(t/5);
% 
% n = length(t);

% w = sqrt(Q)*randn(n,1);
% v = sqrt(R)*randn(n,1);


% %ss model 
% x_0 = [MP_0; MC_0; delC_0];
% u_0 = zeros(length(t_PT),1);

%    B = [1;1;1];
%   C = [1 1 1];

% %Generated input
% % Fextra = 1e12*sin(t);
n = length(t);
x = x_0;     % Initial condition on the state
P_0 = B*Q*B';         % Initial error covariance
u_interp = zeros(size(tinterp));
% u = u_interp(1);

ycov = zeros(n,1); 

global ye
ye = zeros(n,1);

options2 = odeset('OutputFcn',...
    @(t,y,flag) kalfun(t,y,flag,P_0,C,R,yv,tinterp,u_interp,B,Q,A) );

[t,x] = ode15s(@(t,x) odefun(t,x,A,B,B,tinterp,w_interp,u_interp),...
                t,x_0,options2);

% Use process noise w and measurement noise v generated above
% sys = ss(A,B,C,0,-1);
% y = lsim(sys,u+w);      % w = process noise
% yv = y + v;             % v = measurement noise
% y = delC_PT;
% yv = delC_PT;

% for i=1:n
 
%   u = interp1(tinterp,u_interp,t(i));
%     
%   % Measurement update
%   Mn = P*C'/(C*P*C'+R);
%   x = x + Mn*(yv(i)-C*x);   % x[n|n]
%   P = (eye(3)-Mn*C)*P;      % P[n|n]
% 
%   ye(i) = C*x;
%   errcov(i) = C*P*C';
% 
%   % Time update
%   x = A(x,B*u)*x + B*u;        % x[n+1|n]
%   P = A(x,B*u)*P*(A(x,B*u))' + B*Q*B';     % P[n+1|n]
% end



% subplot(211)
% plot(t_PT,delC_PT,t_PT,ye)
% subplot(211),
hold on
plot(t,ye,'r-')
title('Time-varying Kalman filter response')
xlabel('No. of samples'), ylabel('Output')
% subplot(212), plot(t,y-yv,'-.',t,y-ye,'-')
% xlabel('No. of samples'), ylabel('Output')

%%
    function A = jacfun(x,u)
        A = [-kbp      kwp   0 ; ...
             -kborg   -kwsil 0 ]; 
        
        A(3,1) = kborg*epsilon/(MC_0 + x(2));
        A(3,2) = -((Fin + u(2))*(delC_in - x(3)) + kborg*epsilon*(MP_0 + x(1)) )...
                 /(MC_0 + x(2))^2;
        A(3,3) = -(Fin + u(2) )/(MC_0 + x(2));
    end
 


end

function dx = odefun(t,x,A,B,G,tinterp,w_interp,u_interp)

w = interp1(tinterp,w_interp,t);
u = interp1(tinterp,u_interp,t);

dx = A(x,B*u)*x + B*u + G*w;

end

function status = kalfun(t,x,flag,P_0,C,R,yv,tinterp,u_interp,B,Q,A)

global ye

persistent i P

switch flag
    case 'init'        
        i = 1;
        t = 0;
        P = P_0;
end

for j = 1:length(t)

switch flag
    case 'done'
        return
    otherwise
        i = i + 1;
end

    
    
t_in = t(j);   
x_in = x(:,j);
    
u = interp1(tinterp,u_interp,t_in);
% yvi = interp1(tinterp,yv,t);

% Measurement update
Mn = P*C'/(C*P*C'+R);
x_in = x_in + Mn*(yv(i)-C*x_in);   % x[n|n]
P = (eye(3)-Mn*C)*P;      % P[n|n]

ye(i) = C*x_in;
errcov(i) = C*P*C';

% Time update
x_in = A(x_in,B*u)*x_in + B*u;        % x[n+1|n]
P = A(x_in,B*u)*P*(A(x_in,B*u))' + B*Q*B';     % P[n+1|n]

end

status = 0;

% P


end

function [t_PT, delC_PT] = plotfun(plotdata)
%%load P-T data and prep it for regression

load('/Users/avivbachan/My_Documents_Google_Drive/Research/PennStatePostdoc/Phanerozoic_d13C/Data/Jon_Data_Including_Neoprot/delC_full.mat')
load('/Users/avivbachan/My_Documents_Google_Drive/Research/PennStatePostdoc/Phanerozoic_d13C/Data/Jon_Data_Including_Neoprot/t_full.mat')

%remove NaN datapoints
nanind = isnan(delC_full);
delC_full(nanind) = [];
t_full(nanind) = [];

%replace duplicates with their average
[t_unique,~,idx] = unique(t_full);
delC_unique = accumarray(idx,delC_full,[],@mean);

%remove outliers (point to point difference larger than plus minus 2 permil)
ind_diff_delC_outliers = find(or(diff(delC_unique)>2,diff(delC_unique)<-2));
delC_no_outliers = delC_unique;
t_no_outliers = t_unique;
delC_no_outliers(ind_diff_delC_outliers) = [];
t_no_outliers(ind_diff_delC_outliers) = [];

%restrict data to P-T interval
idx_PT = find(and(t_no_outliers>245,t_no_outliers<255));
t_PT = t_no_outliers(idx_PT);
delC_PT = delC_no_outliers(idx_PT);

t_PT = -1e6*t_PT;

%plot data
switch plotdata
    case 'on'
        scrsz = get(0,'ScreenSize');

        h.fig = figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2]);

        h.lines = plot(t_PT,delC_PT,'o');

        hold on

        h.ax = gca;

        xlabel('Time [Ma]')

        ylabel('\delta^{13}C  [permil]')

        legend off

        box on

        grid on

        ylim([-2 +8])

        set(h.fig,'PaperPositionMode','auto');

%         title('\delta^{13}C data of the past 1 Ga')
end


end

function histdiffdata(plothistdiffdata,t_PT,delC_PT)

diffdelC = diff(delC_PT);

datarange = min(diffdelC):1:max(diffdelC);

[nelements,centers] = hist(diffdelC,datarange);

fitnorm = fitdist(diffdelC,'Normal');

datarangeinterp = min(diffdelC):.1:max(diffdelC);

evalfitnorm = pdf(fitnorm,datarangeinterp);


switch plothistdiffdata
    case 'on'
        figure
        bar(centers,nelements'/(sum(nelements)))
        hold on
        plot(datarangeinterp,evalfitnorm)
        
        annotation(gcf,'textbox',...
            [0.1414    0.7457    0.1759    0.1643],...
            'String',{...
            ['n = ' num2str(sum(nelements))],...
            ['\mu = ' num2str(fitnorm.mu)],...
            ['\sigma = ' num2str(fitnorm.sigma)], ...
            },...
            'FitBoxToText','on','EdgeColor','none');
end

end

