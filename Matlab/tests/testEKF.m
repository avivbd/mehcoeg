classdef testEKF < matlab.unittest.TestCase
                  
    properties
        TestData
    end
    

    methods (TestMethodSetup)
        
        function setupOnce(testCase)
            %%
                        
            A = [   1.1269      -0.4940    0.1129;
                    1.0000      0          0;
                    0           1.0000     0];
                
            testCase.TestData.A = A;

            B = [-0.3832;
                  0.5919;
                  0.5191];
              
            testCase.TestData.B = B;  

            C = [1 0 0];
            
            testCase.TestData.C = C;  
            
            Q = 1; 
            R = 1;
            
            testCase.TestData.R = R;  
            testCase.TestData.Q = Q;              
            
            t = [0:100]';
            testCase.TestData.t = t;  
            
            yv = sim_sys(testCase);
            testCase.TestData.yv = yv;  
            
        end        
        
    end
    
    methods (Test)        
                
        function testKFDisc(testCase)
%             [f, x0, P0, h, Q, R, sV, xV, zV] = YI_CAO_TEST_DATA();

            [t, y, yv, u, x] = sim_sys(testCase);
            
            [xks_C, ~, yks_C] = control_sys_tutorial_kf(testCase, yv, u);
            
            [xks_Y, ~, yks_Y] = yckf(testCase, yv, u);
            
            [xks_M, ~, yks_M] = mykf(testCase, yv, u);
            
            plot_results=true;
            
            if plot_results
                
                close all
                figure()
                h1 = subplot(4, 1, 1);
                plot(t, y, '-b', ...
                     t, yks_C, '--r', ...
                     t, yks_Y, '.-m', ...
                     t, yks_M, 'x-k')

                title('Time-varying Kalman filter response')
                xlabel('No. of samples'), ylabel('Output')
                legend(h1, 'obs', 'CS filtered obs', ...
                       'YC filtered obs', 'My kf')
                   
                for i = 2:4
                    h = subplot(4, 1, i);
                    plot(t, x(i-1), '-b', ...
                         t, xks_C(i-1, :), '--r', ...
                         t, xks_Y(i-1, :), '.-m', ...
                         t, xks_M(i-1, :), 'x-k')

                    title('Time-varying Kalman filter response')
                    xlabel('No. of samples'), ylabel('Output')
                    legend(h, 'states', ...
                           'CS filtered obs', ...
                           'YC filtered obs', ...
                           'My kf')   
                end
                   
            end
            
            m1 = max(abs(yks_Y - yks_M));
            m2 = max(abs(yks_C - yks_M));
            % todo make different tests for the differt implementations
            testCase.assertLessThan(m1, 0.1)
            testCase.assertLessThan(m2, 0.68)
                        
        end
        
        function testKFCOnt(testCase)
            
            tspan = [0, 3000];
            tinterp = linspace(tspan(1), tspan(2), 100);
            y0 = [2; 0];
            Mu = 1000;
            
            vp_odefun = @(t,y) vanderpoldemo(t,y, Mu);
            obsfun = @(t, y) y(1, :);
            [t,y] = ode15s(vp_odefun, tinterp, y0);
            yv = y + .5*randn(size(y));
            
            P0 = diag([2, 0]);
            R = (0.1)^2;
            Q = (0.1)^2;
            
            dt = tinterp(2) - tinterp(1);
            sol = EKF(vp_odefun, ...
                obsfun, t, yv, y0, P0, R, Q, ...
                'cont-time', ...
                'implementation','regular',  ...
                'theta', 1);
            
            xks = sol.xks_post;
            Pks = sol.Pks_post;
            zks = sol.zks_post;

%             testCase.assertLessThan( y(:, 1), 3.2)
            
            plot_results = true;
%             close all 
            if plot_results
                figure()
                plot(t, y(:, 1), '-')
                hold on
                plot(t, yv(:, 1), '--')
                plot(t, zks(1, :), 'x-k')
                
                title([sprintf('van der Pol Equation, mu =  %i', Mu)])
                xlabel('t')
                ylabel('solution y')
                legend('solution', 'noisy solution', 'filtered solution')
            end

            
            
        end
        
    end
    
end 

function [xks, Pks, yks] = yckf(testCase, yv, u)

    [A, B, C, R, Q, t] = getParams(testCase);
            
    P0 = B*Q*B';         % Initial error covariance
    x0 = zeros(3,1);     % Initial condition on the state
    
    hmeas = @(x) C*x;
    fstate = @(k, x) odefun(k, x, A, B, t, u);
    
    yks = zeros(length(hmeas(x0)), length(t));
    xks = zeros(length(x0), length(t)); 
    Pks = zeros(length(x0), length(x0), length(t));             
    xk = x0;
    Pk = P0;

    for k=1:length(t)
        [xk, Pk, yk] = YI_CAO_EKF(@(x) fstate(k, x), xk, Pk, ...
                                  hmeas, yv(k), Q, R);
        xks(:, k) = xk;
        Pks(:,:, k) = Pk;
        yks(:, k) = yk;
    end
    
end

function [xks, Pks, zks] = mykf(testCase, yv, u)

[A, B, C, R, Q, t] = getParams(testCase);
odefun_ = @(k, x) odefun(k, x, A, B, t, u);
obsfun = @(t, x) C*x;

P0 = B*Q*B';         % Initial error covariance
x0 = zeros(3,1);     % Initial condition on the state
sol = EKF(odefun_, obsfun, t, ...
    yv, x0, P0, R, Q, 'disc-time', 'sqrt');

xks = sol.xks_post;
Pks = sol.Pks_post;
zks = sol.zks_post;
end

function xt_1 = odefun(t, x, A, B, tu, u)
uk = interp1(tu, u, t, 'pchip');
xt_1 = A*x + B*uk;
end

function [xks, Pks, yks] = control_sys_tutorial_kf(testCase, yv, u)
            %%
            [A, B, C, R, Q, t] = getParams(testCase);
            
            P = B*Q*B';         % Initial error covariance
            x = zeros(3,1);     % Initial condition on the state
            
            yks = zeros(1, length(t));
            xks = zeros(length(x), length(t)); 
            Pks = zeros(length(x), length(x), length(t));             

            
            for i = 1:length(t)
                
                % Time update
                x = A*x + B*u(i);        % x[n+1|n]
                P = A*P*A' + B*Q*B';     % P[n+1|n]
                
                % Measurement update
                Mn = P*C'/(C*P*C'+R);
                y = C*x;
                x = x + Mn*(yv(i)-y);   % x[n|n]
                P = (eye(3)-Mn*C)*P;      % P[n|n]
                
                yks(i) = C*x;
                xks(:, i) = x;
                Pks(:,:, i) = P;
                
            end

end

function [f, x0, P0, h, Q, R, sV, xV, zV] = YI_CAO_TEST_DATA()
rng default
n=3;      %number of state
q=0.1;    %std of process 
r=0.1;    %std of measurement
Q=q^2*eye(n); % covariance of process
R=r^2;        % covariance of measurement  
f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
h=@(x)x(1);                               % measurement equation
s=[0;0;1];                                % initial state
x0=s+q*randn(3,1); %initial state          % initial state with noise
P0 = eye(n);                               % initial state covraiance
N=20;                                     % total dynamic steps
xV = zeros(n,N);          %estmate        % allocate memory
sV = zeros(n,N);          %actual
zV = zeros(1,N);
x = x0;
P = P0;

for k=1:N
  z = h(s) + r*randn;                     % measurments
  sV(:,k)= s;                             % save actual state
  zV(k)  = z;                             % save measurment
  [x, P] = YI_CAO_EKF(f,x,P,h,z,Q,R);            % ekf 
  xV(:,k) = x;                            % save estimate
  s = f(s) + q*randn(3,1);                % update process 
end


for k=1:3                                 % plot results
  subplot(3,1,k)
  plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--')
end


end

function [xk, Pk, zk] = YI_CAO_EKF(fstate,xk,Pk,hmeas,zk,Q,R)
% EKF   Extended Kalman Filter for nonlinear dynamic systems
% [x, P] = ekf(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P 
% for nonlinear dynamic system:
%           x_k+1 = f(x_k) + w_k
%           z_k   = h(x_k) + v_k
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R
% Inputs:   f: function handle for f(x)
%           x: "a priori" state estimate
%           P: "a priori" estimated state covariance
%           h: fanction handle for h(x)
%           z: current measurement
%           Q: process noise covariance 
%           R: measurement noise covariance
% Output:   x: "a posteriori" state estimate
%           P: "a posteriori" state covariance
%

% By Yi Cao at Cranfield University, 02/01/2008
%
[x1,A] = jaccsd(fstate,xk);    %nonlinear update and linearization at current state
Pk = A*Pk*A'+Q;                 %partial update
[z1,H] = jaccsd(hmeas,x1);    %nonlinear measurement and linearization
P12 = Pk*H';                   %cross covariance
% K=P12*inv(H*P12+R);       %Kalman filter gain
% x=x1+K*(z-z1);            %state estimate
% P=P-K*P12';               %state covariance matrix
S = chol(H*P12+R);            %Cholesky factorization
U = P12/S;                    %K=U/R'; Faster because of back substitution
xk = x1 + U*(S'\(zk-z1));         %Back substitution to get state update
Pk = Pk - U*U';                   %Covariance update, U*U'=P12/R/R'*P12'=K*P12.
[zk, ~] = jaccsd(hmeas,xk);
end

function [z,A]=jaccsd(fun,x)
% JACCSD Jacobian through complex step differentiation
% [z J] = jaccsd(f,x)
% z = f(x)
% J = f'(x)
%
    z=fun(x);
    n=numel(x);
    m=numel(z);
    A=zeros(m,n);
    h=n*eps;
    for k=1:n
        x1=x;
        x1(k)=x1(k)+h*1i;
        A(:,k)=imag(fun(x1))/h;
    end

end

function [t, y, yv, u, x] = sim_sys(testCase)
            
    [A, B, C, R, Q, t] = getParams(testCase);

    n = length(t);
    u = sin(t/5);

    rng default
    w = sqrt(Q)*randn(n,1);
    v = sqrt(R)*randn(n,1);

    sys = ss(A,B,C,0,-1);
    [y, ~, x] = lsim(sys,u+w);
    x = x';
    yv = y + v;

end
        
function [A, B, C, R, Q, t] = getParams(testCase)
%%
    A = testCase.TestData.A;
    B = testCase.TestData.B;
    C = testCase.TestData.C; 
    R = testCase.TestData.R;  
    Q = testCase.TestData.Q; 
    t = testCase.TestData.t;
end

function dydt = vanderpoldemo(t, y, Mu)
%VANDERPOLDEMO Defines the van der Pol equation for ODEDEMO.

% Copyright 1984-2014 The MathWorks, Inc.
dydt = [y(2); ...
        Mu*(1-y(1)^2)*y(2) - y(1)];

end


















% function CONTROLS_SYS_KF(testCase)
%         
%     [A, B, C, R, Q, w, v, u, t] = getParams(testCase);
% 
%     Plant = ss(A, [B, B], C, 0, -1, ...
%         'inputname', {'u', 'w'}, ...
%         'outputname', 'y');
% 
%     [kalmf, L, P, M] = kalman(Plant, Q, R);
% 
% 
%     kalmf = kalmf(1, :);
% 
%     a = A;
%     b = [B B 0*B];
%     c = [C; C];
%     d = [0 0 0; 0 0 1];
%     P = ss(a, b, c, d, -1, ...
%         'inputname', {'u', 'w', 'v'}, ...
%         'outputname', {'y', 'yv'});
% 
%     sys = parallel(P, kalmf, 1, 1, [], []);
% 
%     SimModel = feedback(sys, 1, 4, 2, 1);
%     SimModel = SimModel([1 3], [1 2 3]);
% 
%     [out, x] = lsim(SimModel, [w, v, u]);
% 
%     y = out(:,1);   % true response
%     ye = out(:,2);  % filtered response
%     yv = y + v;     % measured response
% 
%     f = figure();
%     subplot(211), plot(t,y,'--',t,ye,'-'), 
%     xlabel('No. of samples'), ylabel('Output')
%     title('Kalman filter response')
%     legend('y true', 'y filtered')
%     subplot(212), plot(t,y-yv,'-.',t,y-ye,'-'),
%     xlabel('No. of samples'), ylabel('Error')
%     legend('measured err', 'estim err')
%         
% end

    





% (x_t_1 - xt)/dt = J*x
% x_t_1 = x_t + dt*J*x
% x_t_1 = J*x
% function dy = odefun(t,y)
% 
% 
% y = y(:);
% 
% dy = A*y + b;
% z = C*y;
% 
% dy = dy(:);
% 
% end