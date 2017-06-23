%%
n=4; %number of state 
q=0.1; %std of process 
r=0.1; %std of measurement 
Q=q^2*eye(n); % covariance of process 
% or Q=diag[Q 0]; % if no process noise is included in the parameter 
R=r^2; % covariance of measurement 
f=@(x)[x(2);x(3);x(4)*x(1)*(x(2)+x(3));x(4)]; % nonlinear state equations 
h=@(x)x(1); % measurement equation 
s=[0;0;1;0.1]; % initial state 
x=s+q*randn(4,1); %initial state % initial state with noise 
P = eye(n); % initial state covraiance 
N=100; % total dynamic steps 
xV = zeros(n,N); %estmate % allocate memory 
sV = zeros(n,N); %actual 
zV = zeros(1,N); 
for k=1:N 
z = h(s) + r*randn; % measurments 
sV(:,k)= s; % save actual state 
zV(k) = z; % save measurment 
[x, P] = ekf(f,x,P,h,z,Q,R); % ekf 
xV(:,k) = x; % save estimate 
s = f(s) + q*randn(4,1); % update process 
end;


for k=1:4                                 % plot results
  subplot(4,1,k)
  plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--')
end