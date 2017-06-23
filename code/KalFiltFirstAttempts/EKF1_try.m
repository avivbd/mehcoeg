function EKF1_try

clear
clc
close all

plotdata = 'on'; %on/off

[t, Y] = plotfun(plotdata);

hold on

h_func = @ekf_sine_h;
dh_dx_func = @ekf_sine_dh_dx;

% Initial guesses for the state mean and covariance.
M = [0 0 0]';
P = diag([3 3 3]);

R = 0.05;

MM = zeros(size(M,1),size(Y,2)); 
PP = zeros(size(M,1),size(M,1),size(Y,2));

% Estimate with EKF
for k=1:length(Y)
    %prediction step
    %P = A * P * A' + W * Q * W';
    [M,P] = ekf_predict1(M,P,[],[]);
    
    %update steps
%   S = (V*R*V' + H*P*H');
%   K = P*H'/S;
%   M = M + K * (y-MU);
%   P = P - K*S*K';
    
    [M,P] = ekf_update1(M,P,Y(k),dh_dx_func,R*eye(1),h_func);
    
    %accumulate output
    MM(:,k) = M;
    PP(:,:,k) = P;
end

Y_m = feval(h_func, MM);
plot(t,Y_m,'r');

[SM1,SP1] = erts_smooth1(MM,PP,A,Q);
Y_s1 = feval(h_func, SM1);

plot(t,Y_s1,'g');

legend('Noisy measurements','model prediction','EKF','RTS')

function Y = ekf_sine_h(x,~)
%first order EKF demo

f = x(1,:);
a = x(3,:);
Y = a.*sin(f);

function dY = ekf_sine_dh_dx(x, ~)
f = x(1,:);
w = x(2,:);
a = x(3,:);
dY = [(a.*cos(f))' zeros(size(f,2),1) (sin(f))'];


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

