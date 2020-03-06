classdef testInverter < matlab.unittest.TestCase
    
    properties
        TestData
    end
    
    methods (TestMethodSetup)
        function setupOnce(testCase)
            
            J = [(-2e-6 -2e-3) 2e-6; ...
                    2e-3      -2e-6];
                
            testCase.TestData.J = J;

            F = [0 ;0];
            
            testCase.TestData.F = F;
            
            M0 = [5e4; 5e7];
            
            testCase.TestData.M0 = M0;
            
        end
    end
    
    methods (Test)
        
        function [J, F, M0] = getParams(testCase)
            J = testCase.TestData.J;
            F = testCase.TestData.F;
            M0 = testCase.TestData.M0;
        end
                                    
        function testODEInverter(testCase) 
            
            % for a single large timestep inverting a stiff system 
            % using simple substitution (forward Euler) results in 
            % a horribly inaccurate delta. 
            % In contrast, the stable inversion results in a dy
            % much closer to the the symbolic result. 
            
            [J, F, M0] = getParams(testCase);
            odefunc = @(t, y) odefun(t, y, J, F);
            t = 0; dt = 1e6;
            dy_forward_euler = dt*odefunc(t, M0);
            dy_crank_nick = ode_inverter(odefunc, t, M0, dt);
            
            syms M1(t) M2(t)
            Y = [M1; M2];
            odes = (diff(Y) == J*Y + F);
            C = (Y(0) == M0);
            [M1(t), M2(t)] = dsolve(odes, C);
            
            symb_M1 = double(subs(M1(t), t, dt)) - double(subs(M1(t), t, 0));
            symb_M2 = double(subs(M2(t), t, dt)) - double(subs(M2(t), t, 0));
            
            dy_symb = [symb_M1; symb_M2];
            
            cond = ...
            (sum(abs(dy_symb - dy_forward_euler)) > 1.99e5) & ...
            (sum(abs(dy_symb - dy_crank_nick) < 100)) ;
            
            testCase.assertTrue(cond)
            
        end
        
    end
    
end 


function dy = odefun(t, y, J, F)

y = y(:);
    
dy = J*y + F;

dy = dy(:);

end



