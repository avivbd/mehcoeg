function EKF1_demo

clear
clc
close all


h_func = @ekf_sine_h;
dh_dx_func = @ekf_sine_dh_dx;



% Initial values for the signal.
f = 0;
w = 10;
a = 1;
  
% Number of samples and stepsize.
d = 5;
n = 500;
dt = d/n;
x = 1:n;

% der_check(h_func, dh_dx_func, 1, [f w a]');


% Dynamic state transition matrix in continous-time domain.
F = [0 1 0;
     0 0 0;
     0 0 0];
  
% Noise effect matrix in continous-time domain.
L = [0 0;
     1 0;
     0 1];
  
% Spectral power density of the white noise.
q1 = 0.2;
q2 = 0.1;
Qc = diag([q1 q2]);
  
% Discretize the plant equation.
[A,Q] = lti_disc(F,L,Qc,dt);
  
% Generate the real signal.
X = zeros(3, n);
X(:,1) = [f w a]';
for i = 2:n
   X(:,i) = A*X(:,i-1) + gauss_rnd([0 0 0]', Q);
end  
  
% Generate the observations with Gaussian noise.
sd = 1;
R = sd^2;

Y = zeros(1,n);
Y_real = feval(h_func,X);     
Y = Y_real + gauss_rnd(0,R,n);

plot(x,Y,'.',x,Y_real)
hold on


% Initial guesses for the state mean and covariance.
M = [f w a]';
P = diag([3 3 3]);

MM = zeros(size(M,1),size(Y,2)); 
PP = zeros(size(M,1),size(M,1),size(Y,2));

% Estimate with EKF
for k=1:size(Y,2)
[M,P] = ekf_predict1(M,P,A,Q);
[M,P] = ekf_update1(M,P,Y(:,k),dh_dx_func,R*eye(1),h_func);
MM(:,k) = M;
PP(:,:,k) = P;
end

Y_m = feval(h_func, MM);
plot(x,Y_m,'r');

[SM1,SP1] = erts_smooth1(MM,PP,A,Q);
Y_s1 = feval(h_func, SM1);

plot(x,Y_s1,'g');

legend('Noisy measurements','model prediction','EKF','RTS')

function Y = ekf_sine_h(x,param)
%first order EKF demo

f = x(1,:);
a = x(3,:);
Y = a.*sin(f);
if size(x,1) == 7
Y = Y + x(7,:);
end

function dY = ekf_sine_dh_dx(x, param)
f = x(1,:);
w = x(2,:);
a = x(3,:);
dY = [(a.*cos(f))' zeros(size(f,2),1) (sin(f))'];