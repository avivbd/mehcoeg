function  sol = EKF(odefun, obsfun, t, z, x0, P0, R, Q, ...
        state_eqn_type, varargin)
%EKF Implementation of the Kalman filter
%   This code implements the extended Kalman filter
%   with discrete time measurements
%   choose between continous time and discrete time state eqs
%   eg x_dot=f(x, t) *or* x(t+1) = f(x(t), t)


p = verify_params(odefun, obsfun, z, t, x0, P0, R, Q,...
                  state_eqn_type, varargin{:});
              
sol = main(p.Results.odefun, ...
                       p.Results.obsfun, ...
                       p.Results.z, ...
                       p.Results.t, ...
                       p.Results.x0, ...
                       p.Results.P0, ...
                       p.Results.R, ...
                       p.Results.Q, ...
                       p.Results.state_eqn_type,...
                       p.Results.implementation, ...
                       varargin{:});
                
                   
end 

function p = verify_params(odefun, obsfun, z, t, x0, P0, R, Q, ...
    state_eqn_type, varargin)
%% Parse Inputs
p = inputParser();
p.KeepUnmatched = true;
% required arguments
isFunctionHandle = @(x) isa(x, 'function_handle') & (nargin(x) == 2);
addRequired(p, 'odefun', isFunctionHandle)
addRequired(p, 'obsfun', isFunctionHandle)

addRequired(p, 'z', @isnumeric)
validationFcn = @(x) isnumeric(x);
addRequired(p, 't', validationFcn)
addRequired(p, 'x0', @isnumeric)
addRequired(p, 'P0', @isnumeric)
addRequired(p, 'R', @isnumeric)
addRequired(p, 'Q', @isnumeric)

valid_types = {'disc-time', 'cont-time'};
isValidType = @(state_eqn_type) any(strcmp(state_eqn_type, valid_types));
addRequired(p, 'state_eqn_type', isValidType)

valid_implementation = {'regular', 'sqrt'};
isValidImp = @(implementation) any(strcmp(implementation, valid_implementation));
default_validation = 'regular';
addOptional(p, 'implementation', default_validation, isValidImp)

parse(p, odefun, obsfun, z, t, x0, P0, R, Q, ...
      state_eqn_type, varargin{:})

end

function sol = main(odefun, obsfun, z, t, x0, P0, R, Q,...
    state_eqn_type, implementation, varargin)

len_t = length(t);
len_x = length(x0);
len_z = length(obsfun(t(1), x0));
P_shape = size(P0);

xks_post = zeros(len_x, len_t);
Pks_post = zeros(P_shape(1), P_shape(2), len_t);
zks_post = zeros(len_z, len_t);

xks_pri = xks_post;
Pks_pri  = Pks_post; 
zks_pri = zks_post;

Fks = Pks_pri;
I = eye(P_shape);

xk = x0(:);
Pk = P0;

    for i=1:length(t)

        tk = t(i);
        
        if i==1
            dt = t(2) - t(1);
        else
            dt = t(i) - t(i-1);
        end
        
        zk = z(i);
        
        if strcmp(state_eqn_type, 'cont-time')
            x_hat = ode_incrementor(odefun, tk, xk, dt, varargin{:});
            xks_pri(:, i) = x_hat;
            
            F = numjac(odefun, tk, x_hat);
            Fks(:, :, i) = F;
            
            P_dot = @(tk, Pk) mRiccati(tk, Pk, F, Q); 
            P_hat = ode_incrementor(P_dot, tk, Pk, dt, varargin{:});
            Pks_pri(:, :, i) = P_hat;

            if any(isnan(P_hat))
                error('Phat is a singular matrix!')
            end
            
        else
            x_hat = odefun(tk, xk);
            xks_pri(:, i) = x_hat;
            
            F = numjac(odefun, tk, xk);
            Fks(:, :, i) = F;
            
            P_hat = F*Pk*F' + Q;
            Pks_pri(:, :, i) = P_hat;
        end
        
        
        if strcmp(implementation, 'regular')
            H = numjac(obsfun, tk, x_hat); 
            S = (H*P_hat*H' + R);
            K = P_hat*H'*pinv(S);

            zhat = obsfun(tk, x_hat);
            zks_pri(:, i) = zhat;
            
            xk = x_hat + K*(zk - zhat);
            Pk = (I - K*H)*P_hat;
            zk = obsfun(tk, xk);
        else
            H = numjac(obsfun, tk, x_hat);
            P12 = P_hat*H';
            S = chol(H*P12+R);
            U = P12/S;
            
            zhat = obsfun(tk, x_hat);
            zks_pri(:, i) = zhat;
            
            xk = x_hat + U*(S'\(zk-zhat));
            Pk = P_hat - U*U';
            zk = obsfun(tk, xk);
        end
        
        xks_post(:, i) = xk;
        Pks_post(:, :, i) = Pk;
        zks_post(:, i) = zk;
        

    end

    
    sol = struct();
    v = {xks_post, Pks_post, xks_pri, ...
         Pks_pri, zks_post, zks_pri, ...
         Fks};
    k = {'xks_post', 'Pks_post', 'xks_pri', ...
         'Pks_pri', 'zks_post', 'zks_pri', ...
         'Fks'};
    
    for j=1:length(v)
        sol.(k{j}) = v{j};
    end
    
    
end

function dPdt = mRiccati(~, P, F, Q)
P = reshape(P, size(F)); 
dPdt = F*P + P*F' + Q; 
dPdt = dPdt(:); 
end

