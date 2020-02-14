function [sol] = homebrew_solver(t, y0, odefun)
% homebrew_solver accepts a vector t of timepoints, a vector of n initial 
% values and an ode function of the form dx = fun(x, t)   


%twiggle parameter (balance implicit and explicit)
theta = 0.5 ;

%inital time step size
% dt = 0.005;

% max iterations
% max_iter = 50000;

%%%%%%%%%%%%%%%%%%
%parse arguments

%%%%%%%%%%%%%%%
% intialize values for main loop
tstart = t(1);
tend = t(end);

tk = tstart;
k = 0;
y = y0(:);

% initialize variables to collect results
% Y = 
% T = 

not_done = true;

while not_done
    
    dy = main_loop_body(tk, y, theta, dt);
    
    %add to previous value of y
    y = y + dy;
    tkplus1 = tk + dt;
    tk = tkplus1;
    
    if tk >= tend
        break
    end
    
%     if k >= max_iter
%         print('Max iterations reached!')
%         break
%     end
    
end    


    %variable timestep
%     phi = max( abs((Y(:,n)-Y(:,n-1))./Y(:,n-1)) );

%         if phi < 0.01
%             dt = dt*2;

%         else
%             dt = dt/2;

%         end

%     n = n+1;


end


function [dy] = main_loop_body(tk, y, dt, theta)

    %evaluate derivative
    df = odefun(tk, y);

    %evaluate the jacobian
    J = num_jac(y, tk);

    %calculate change in function value (dy)
    dy = (eye(size(J))/dt-theta*J)\df;

    
end

function [ J ] = num_jac(y, tk)
%NUM_JAC Calculate numerical Jacobian of a function

%how sensitve is the jacobian
eps = 0.000001;

% options for method: 'complex-step', 'first-order-fd', 'second-order-fd'
method = 'complex-step';

J = zeros(length(y));

%calculate Jacobian
    for ii = 1:length(y) %perturb one value of y at a time
        h = y(ii)*eps; % get diff relative to y_i
        
        switch method
            case 'first-order-fd'
                df_y = odefun(tk, y); 
                y_plus_h = y;        
                y_plus_h(ii) = y(ii) + h;
                df_y_plus_h = odefun(tk, y_plus_h);
                J(:,ii)  = (df_y_plus_h - df_y)/h; 
                
            case 'second-order-fd'
                y_plus_h = y;        
                y_plus_h(ii) = y(ii) + h;
                df_y_plus_h = odefun(tk, y_plus_h);
                
                y_minus_h = y;
                y_minus_h(ii) = y(ii) - h;
                df_y_minus_h = odefun(tk, y_minus_h);
                J(:, ii) = (df_y_plus_h - df_y_minus_h)/(2*h);
                
            case 'complex-step'
                y_plus_ih = y;        
                y_plus_ih(ii) = y(ii) + 1i*h;
                df_y_plus_ih = odefun(tk, y_plus_ih);
                J(:, ii) = imag(df_y_plus_ih)/h;     
        end
    end
    
end


