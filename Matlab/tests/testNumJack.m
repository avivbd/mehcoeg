classdef testNumJack < matlab.unittest.TestCase
    
    methods (Test)
                                    
        function testNumJacSingle(testCase) 
            for method={'complex-step', 'first-order-fd', 'second-order-fd'}
                y = pi/4;
                Fn = @(t, x) exp(x)./((cos(x)).^3 + (sin(x)).^3);
                flnum = numjac(Fn, 1, y, 'method', method{1});

                syms x
                Fs = exp(x)/((cos(x))^3 + (sin(x))^3);
                Fps = diff(Fs);
                exact = subs(Fps, y);
                flexact = double(exact);

                testCase.assertEqual(flnum, flexact, ...
                    'RelTol', 1e-7, ...
                    'AbsTol', 1e-7)
            end
            
        end
        
        
        function testNumJacVec(testCase) 
            for method={'complex-step', 'first-order-fd', 'second-order-fd'}
                y = [1; 1]; t = 2;
                jnum = numjac(@odefun, t, y, 'method', method{1});

                syms f(y1, y2) 
                f(y1, y2) = odefun(t, [y1, y2]);
                jsym = jacobian(f(y1, y2));
                jexact = subs(jsym, {y1, y2}, [y(1), y(2)]);

                diff = double(sum(sum(jexact - jnum)));
                testCase.assertTrue(diff<1e-7)
            end
            
        end
        
        
        
    end
    
end 


function dy = odefun(t,y)

y = y(:);

dy(1) = -10*y(1)*t;
dy(2) = -12*y(2)*y(1);

dy = dy(:);

end

