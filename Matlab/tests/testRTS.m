classdef testRTS < matlab.unittest.TestCase
                  
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
            
            t = 1:100;
            testCase.TestData.t = t;  
            
        end        
        
    end
    
    methods (Test)        
                
        function testKFDisc(testCase)

            
            [A, B, C, R, Q, t] = getParams(testCase);
            
            [y, yv, u, x, w] = sim_sys(testCase, @odefun, A, B, C, R, Q, t);
            
            obsfun = @(tt, x) C*x;
            odefun_ = @(tt, x) odefun(tt, x, A, B, t, u, w);

            P0 = B*Q*B';         % Initial error covariance
            x0 = zeros(3,1);     % Initial condition on the state
            
            sol = EKF(odefun_, obsfun, t, ...
                                yv, x0, P0, R, Q, 'disc-time');

%             [xks_rts, Pks_rts, zks_rts] = 
            [xks_rts, Pks_rts] = RTS_smoother(sol.xks_post, sol.Pks_post, ...
                         sol.xks_pri, sol.Pks_pri, ...
                         sol.Fks);

            
            plot_results=true;
            
            if plot_results
                 figure()
                 hold on 
                 plot(t, y, '-k', 'LineWidth', 1)
                 plot(t, yv, '--xb', 'LineWidth', 1)
                 plot(t, sol.zks_post, '-r', 'LineWidth', 1)
                 plot(t, obsfun(t, xks_rts), '-g', 'LineWidth', 1)
                 title('Time-varying Kalman filter response')
                 xlabel('No. of samples'), ylabel('Output')
                 s = {'truth','obs', 'kf filtered obs', 'rts smoothed obs'};
                 legend(s)
                     
            end
            
%             m1 = max(abs(yks_Y - yks_M));
%             m2 = max(abs(yks_C - yks_M));
%             % todo make different tests for the differt implementations
%             testCase.assertLessThan(m1, 0.22)
%             testCase.assertLessThan(m2, 0.31)
                        
        end
        
    end
    
end 

function [y, yv, u, x, w] = sim_sys(testCase, odefun, A, B, C, R, Q, t)
            
    n = length(t);
    u = sin(t/5);

    rng default
    w = sqrt(Q)*randn(n,1);
    v = sqrt(R)*randn(n,1);

    xs = zeros(3, length(t));
    x = xs(:, 1);
    for i=1:length(t)
        x = odefun(i, x, A, B, t, u, w);
        xs(:, i) = x;
    end
    
    y = C*xs;
    yv = (y' + v)';

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

function xt_1 = odefun(t, x, A, B, tu, u, w)
uk = interp1(tu, u, t, 'pchip');
wk = interp1(tu, w, t, 'pchip');
xt_1 = A*x + B*(uk + wk);
end




