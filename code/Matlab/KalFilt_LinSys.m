function [] = KalFilt_LinSys
%Matlab code that demonstrates application of kalman filter to smoothing
%and parameter estimation in a simple linear system
%dy1/dt = -2y1 + 2y2 +1
%dy2/dt = -4y1 +2y2
%%
close all 
clear
addpath('/Users/avivbachan/Google_Drive/Software/MATLAB/ekfukfToolbox')

%% This system is simple enough that it has an analytical solution
%y1 = 1.5sin2t + 1.5cos2t ? 0.5
%y2 = 3cos2t -1

t = linspace(0,25,501);
dt = t(2) - t(1);

[t, yan1, yan2] = analyticalSolution(t);

%This solution can be verified using numerical integration
F = [-2, 2; -4, 2];
U = [1;0];
y0 = [1;2];


[t, ynum1, ynum2] = numericalSolution(t,F,U,y0);

%plot analytical and numerical solutions
figure('Name','Analytical Solution and Numerical Verification')
subplot 211
plot(t, yan1,'b-',t,ynum1,'k.')
xlabel('t')
ylabel('y1')
title('Analytical Solution and Numerical Verification')

subplot 212
plot(t,yan2,'g-',t,ynum2,'k.')
xlabel('t')
ylabel('y2')

%% Make some synthetic data by adding noise to the analytical solution

y1 = yan1 + randn(1,length(yan1));
y2 = yan2 + randn(1,length(yan2));
Y = [y1;y2];


%% Show stochastic runs of model

% addpath('/Users/avivbachan/Google_Drive/Software/MATLAB/SDETools-horchler-db4b64c/SDETools')
%discretize the equation for a given dt
[A,Q] = lti_disc(F,[],1,dt);
% 
% [ystoch1, ystoch2] = stochasticSolution(t,A,U,y0,Q);
% 
% figure('Name','Stochastic Model Simulation')
% subplot 211
% plot(t, yan1,'b-',t,ystoch1,'k-')
% xlabel('t')
% ylabel('y1')
% title('Stochastic Model Simulation')
% 
% subplot 212
% plot(t,yan2,'g-',t,ystoch2,'k-')
% xlabel('t')
% ylabel('y2')
% 

%% Now for kalman filtering

%initialize
m = y0;
P = eye(2);
H = eye(2);
R = 0.1*eye(2);

%preallocate means and variances
MM = zeros(size(m,1), size(Y,2));
PP = zeros(size(m,1), size(Y,2));


for i = 1:size(Y,2)
    [m,P] = kf_predict(m,P,A,Q,[],U*dt);
    [m,P] = kf_update(m,P,Y(:,i),H,R);
    MM(:,i) = m;
    PP(:,i) = diag(chol(P));
end

yfilt1 = MM(1,:);
yfilt2 = MM(2,:);
sigmaFilt1 = PP(1,:);
sigmaFilt2 = PP(2,:);


%plot the filter solution and the noisy data

figure('Name','Noisy Data and Filter Solution')
subplot 211
plot(t, y1,'k.',t,yfilt1,'b-')
hold on
plot(t, yfilt1-sigmaFilt1,'r--')
plot(t, yfilt1+sigmaFilt1,'r--')
xlabel('t')
ylabel('y1')
title('Noisy Data and Filter Solution')

subplot 212
plot(t,y2,'k.',t,yfilt2,'g-')
hold on
plot(t, yfilt2+sigmaFilt2,'r--')
plot(t, yfilt2-sigmaFilt2,'r--')
xlabel('t')
ylabel('y2')


%% Parameter estimation

F = [-2, 2, 1;...
    -4, 2, 0 ;...
     0, 0, 0 ];
[A,Q] = lti_disc(F,[],1,dt);
y0 = [1;2;5];
m = y0;
P = eye(3);
H = [1, 0, 0 ; 0, 1, 0];
R = 0.1*eye(2);

%preallocate means and variances
MM = zeros(size(m,1), size(Y,2));
PP = zeros(size(m,1), size(Y,2));


for i = 1:size(Y,2)
    [m,P] = kf_predict(m,P,A,Q,[],[]);
    [m,P] = kf_update(m,P,Y(:,i),H,R);
    MM(:,i) = m;
    PP(:,i) = diag(chol(P));
end


yfilt1 = MM(1,:);
yfilt2 = MM(2,:);
yfilt3 = MM(3,:);
sigmaFilt1 = PP(1,:);
sigmaFilt2 = PP(2,:);
sigmaFilt3 = PP(3,:);


figure('Name','Filter Solution and Parameter Estimation')
subplot 311
plot(t, y1,'k.',t,yfilt1,'b-')
hold on
plot(t, yfilt1-sigmaFilt1,'r--')
plot(t, yfilt1+sigmaFilt1,'r--')
xlabel('t')
ylabel('y1')
title('Filter Solution and Parameter Estimation')

subplot 312
plot(t,y2,'k.',t,yfilt2,'g-')
hold on
plot(t, yfilt2+sigmaFilt2,'r--')
plot(t, yfilt2-sigmaFilt2,'r--')
xlabel('t')
ylabel('y2')


subplot 313
plot(t,yfilt3,'k-')
hold on
plot(t, yfilt3+sigmaFilt3,'r--')
plot(t, yfilt3-sigmaFilt3,'r--')
xlabel('t')
ylabel('y2')


end

function [t, yan1, yan2] = analyticalSolution(t)
%% Analytical solution
yan1 = (3*sin(2*t) + 3*cos(2*t) -1)/2;
yan2 = 3*cos(2*t) - 1;

end



function [t, ynum1, ynum2] = numericalSolution(t,F,u,y0)
%% Numerical solution
odefun = @(t,y) F*y+u;
[~,Y] = ode45(odefun, t, y0);
ynum1 = Y(:,1);
ynum2 = Y(:,2);

end


function [t, ystoch1, ystoch2] = stochasticSolution(t,F,u,y0,Q)
%% Numerical solution

FFUN = @(t,y) F*y+u;
GFUN = @(t,y) Q;
options = sdeset('DiagonalNoise','no');
Y = sde_euler(FFUN,GFUN,t,y0,options);
ystoch1 = Y(:,1);
ystoch2 = Y(:,2);

end
