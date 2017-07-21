%% Kalman filtering example from Poor Mans Explanation
clc
clear
close all

cd('/Users/avivbachan/My_Documents_Google_Drive/Research/StanfordPostdoc/KalmanFiltering')

img = imread('Kim_U_paper_standards_fig.jpg');

img(1:4180,:,:) = [];

h.image = imshow(img,'InitialMagnification',33);
h.fig = gcf;
hold on

%DWP means
DWP.mu = [ ...
    0.0578   -0.3165
    0.1168   -0.3653
    0.1751   -0.1730
    0.2356   -0.3162
    0.2940   -0.1582
    0.3531   -0.2476
    0.4120   -0.2443
    0.4697   -0.2442
    0.5300   -0.3549
    0.5899   -0.3059
    0.6478   -0.3677
    0.7075   -0.2617
    0.7650   -0.1672
    0.8252   -0.2534
    0.8845   -0.3836
    0.9430   -0.2662 ];

%standard deviatiosn in the plus direction
DWP.sd_plus = [ ...
    0.0580   -0.2340
    0.1176   -0.3063
    0.1753   -0.1089
    0.2351   -0.2321
    0.2941   -0.1103
    0.3514   -0.1859
    0.4120   -0.1907
    0.4703   -0.1808
    0.5300   -0.3040
    0.5892   -0.1937
    0.6487   -0.2347
    0.7080   -0.1968
    0.7661   -0.1211
    0.8257   -0.1983
    0.8843   -0.3264
    0.9434   -0.1964];

%standard deviatiosn in the minus direction
DWP.sd_minus = [ ...
    0.0584   -0.3979
    0.1173   -0.4239
    0.1756   -0.2365
    0.2348   -0.3946
    0.2942   -0.2056
    0.3534   -0.3053
    0.4122   -0.2976
    0.4710   -0.3098
    0.5301   -0.4002
    0.5890   -0.4217
    0.6492   -0.4953
    0.7075   -0.3232
    0.7660   -0.2156
    0.8251   -0.3046
    0.8843   -0.4396
    0.9429   -0.3382];

h.ax = axes('Parent',h.fig, 'Position',[0.1259    0.1675    0.8124    0.7618]);
plot(h.ax,DWP.mu(:,1),DWP.mu(:,2),'.')
set(h.ax,'Color', 'none','XTick',[],'YTick',[])
ylim([-0.5275 -0.095])


n = length(DWP);
% initial estimate of measurement precision
% take the mean of the sample sd's
sigma_measurement_0 = mean(DWP.sd_plus(:,2) - DWP.sd_minus(:,2));
% initial estimate of measurement value
% take the mean of the sample means
mu_0 = mean(DWP.mu(:,2));
sigma_process_0 =  0.01;

%%

b = zeros(1,n); 
x_hat = zeros(1,n); 
sigma_xhat = zeros(1,n);

%steps are:
% 1. compute weighting coefficient
b(1) = sigma_x^2/(sigma_x^2 + sigma_eps^2);
% 2. use wighting coefficient to make new estimate 
x_hat(1) = mu_x + b(1)*(y(1)  - mu_x);
% 3. update mean square error in estimation
sigma_xhat(1) = realsqrt((1 - b(1))*sigma_x^2);


for i = 2:length(y), 
    b(i) = sigma_xhat(i-1)^2/(sigma_xhat(i-1)^2 + sigma_eps^2); 
    x_hat(i) = x_hat(i-1) +b(i)*(y(i) - x_hat(i-1)); 
    sigma_xhat(i) = realsqrt((1 - b(i))*sigma_xhat(i-1));
end

subplot 211
plot(1:length(x_hat),x_hat,1:length(y),y)
legend('Optimal estimate','measured value')

subplot 212
plot(1:length(sigma_xhat),sigma_xhat)

%%
% %% Matlab example
% 
% clear 
% 
% A = [1.1269   -0.4940    0.1129
%      1.0000         0         0
%           0    1.0000         0];
% 
% B = [-0.3832
%       0.5919
%       0.5191];
% 
% C = [1 0 0]; 
% 
% x = zeros(3,1);     % Initial condition on the state
% 
% t = [0:100]';
% % u = 10*sin(2*pi*t/50);
% u = sin(t/5);
% 
% n = length(t);
% randn('seed',0)
% 
% Q = 1; R = 1;
% 
% w = sqrt(Q)*randn(n,1);
% v = sqrt(R)*randn(n,1);
% 
% % Use process noise w and measurement noise v generated above
% sys = ss(A,B,C,0,-1);
% y = lsim(sys,u+w);      % w = process noise
% yv = y + v;             % v = measurement noise
% 
% P = B*Q*B';         % Initial error covariance
% ye = zeros(length(t),1);
% ycov = zeros(length(t),1); 
% 
% for i=1:length(t)
%   % Measurement update
%   Mn = P*C'/(C*P*C'+R);
%   x = x + Mn*(yv(i)-C*x);   % x[n|n]
%   P = (eye(length(x))-Mn*C)*P;      % P[n|n]
% 
%   ye(i) = C*x;
%   errcov(i) = C*P*C';
% 
%   % Time update
%   x = A*x + B*u(i);        % x[n+1|n]
%   P = A*P*A' + B*Q*B';     % P[n+1|n]
% end
% 
% figure
% subplot(211), 
% plot(t,y,'--',t,ye,'-')% ,t,yv,'-.')
% legend('true','estimated','measured')
% 
% subplot(212), plot(t,y-yv,'-.',t,y-ye,'-')
% xlabel('No. of samples'), ylabel('Output')
% legend('estimated','measured')
% 
% title('Time-varying Kalman filter response')
% xlabel('No. of samples'), ylabel('Output')
% % subplot(212), plot(t,y-yv,'-.',t,y-ye,'-')
% % xlabel('No. of samples'), ylabel('Output')
% 
% %% DSP Kalman object example
% 
% numSamples = 4000;
% R = 0.02;
% hSig = dsp.SignalSource;
% hSig.Signal = [ones(numSamples/4,1);   -3*ones(numSamples/4,1);...
%               4*ones(numSamples/4,1); -0.5*ones(numSamples/4,1)];
%       
% hTScope = dsp.TimeScope('NumInputPorts', 3, 'TimeSpan', numSamples, ...
%           'TimeUnits', 'Seconds', 'YLimits',[-5 5], ...
%           'ShowLegend', true); % Create the Time Scope
% hKalman = dsp.KalmanFilter('ProcessNoiseCovariance', 0.0001,...
%           'MeasurementNoiseCovariance', R,...
%           'InitialStateEstimate', -5,...
%           'InitialErrorCovarianceEstimate', 1,...
%           'ControlInputPort',false); %Create Kalman filter
%                       
% while(~isDone(hSig))
%    trueVal = step(hSig);
%    noisyVal = trueVal + sqrt(R)*randn;
%    estVal = step(hKalman, noisyVal);
%    step(hTScope,noisyVal,trueVal,estVal);
% end                      
                      
                      
