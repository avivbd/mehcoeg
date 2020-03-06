function  [X0, Xs, X_hats, Ds, Ts] = ENKF(odefun, H, t, z, x0, P0, R, Q, ...
                                 state_eqn_type, cov_or_sd, varargin)
%ENKF Implementation of the Kalman filter
%   This code implements the ensemble Kalman filter
%   with discrete time measurements
%   choose between continous time and discrete time state eqs
%   eg x_dot=f(x, t) *or* x(t+1) = f(x(t), t)
% Inputs:

% odefun: function of the form f(t, x)
% H: observation matrix. Which of the states is observable?  
% t: time vector for observations
% z: observations
% x0: initial guess for state
% P0: initial sd matrix for state
% R: observation/measurement uncertainty
% Q: process uncertainty 
% state_eqn_type: continous (cont-time) or discrete (disc-time). 
% If continous expected form is dx/dt = f(t,x)
% If discrete expected form is x(t+1) = f(t, x(t))
% cov_or_sd: whether matrices for uncertainty ellipses are given as 
% covariance (in which case they must be symmetric and positive definite)
% or standard deviations. 

% returns: 
% X0: initial ensemble (particles) perturbed with noise drawn from P0
% Xs: state estimate (postirori, after data assimilation)
% X_hats: predicted states (priori, prior to data assimilation)
% Ds: ensemble of perturbed observations with noise drawn from R

p = verify_params(odefun, H, z, t, x0, P0, R, Q,...
                  state_eqn_type, cov_or_sd, varargin{:});
              
[X0, Xs, X_hats, Ds, Ts] = main(p.Results.odefun, ...
                       p.Results.H, ...
                       p.Results.z, ...
                       p.Results.t, ...
                       p.Results.x0, ...
                       p.Results.P0, ...
                       p.Results.R, ...
                       p.Results.Q, ...
                       p.Results.state_eqn_type,...
                       p.Results.implementation, ...
                       p.Results.ensemble_size, ...
                       p.Results.cov_or_sd, ...
                       varargin{:});

end 

function p = verify_params(odefun, H, z, t, x0, P0, R, Q, ...
    state_eqn_type, cov_or_sd, varargin)
%% Parse Inputs
p = inputParser;
p.KeepUnmatched = true;

% required arguments
isFunctionHandle = @(x) isa(x, 'function_handle') & (nargin(x) == 2);
addRequired(p, 'odefun', isFunctionHandle)
addRequired(p, 'H', @isnumeric)

addRequired(p, 'z', @isnumeric)
addRequired(p, 't', @isnumeric)
addRequired(p, 'x0', @isnumeric)
addRequired(p, 'P0', @isnumeric)
addRequired(p, 'R', @isnumeric)
addRequired(p, 'Q', @isnumeric)

valid_types = {'disc-time', 'cont-time'};
isValidType = @(state_eqn_type) any(strcmp(state_eqn_type, valid_types));
addParameter(p, 'state_eqn_type', isValidType)

valid_implementation = {'pert-obs', 'sqrt'};
isValidImp = @(implementation) any(strcmp(implementation, valid_implementation));
default_validation = 'pert-obs';
addParameter(p, 'implementation', default_validation, isValidImp)

default_ensemble_size = 10*numel(x0);
ensemble_validation_fcn = @(x) isscalar(x) && (x >= 0) ;
addParameter(p, 'ensemble_size', default_ensemble_size, ...
           ensemble_validation_fcn)
       
valid_cov_or_sd = {'cov', 'sd'};
isValidCov = @(cov_or_sd) any(strcmp(cov_or_sd, valid_cov_or_sd));
default_cov_or_sd = 'cov';
addParameter(p, 'cov_or_sd', default_cov_or_sd, isValidCov)
       
parse(p, odefun, H, z, t, x0, P0, R, Q, ...
      state_eqn_type, cov_or_sd, varargin{:})

end

function [X0, Xs, X_hats, Ds, Ts] = main(odefun, H, z, t, x0, P0, R, Q,...
    state_eqn_type, implementation, ensemble_size, cov_or_sd, varargin)

len_t = length(t);
len_x = length(x0);
len_z = length(H*x0);

Xs = zeros(len_x, ensemble_size, len_t);
X_hats = zeros(len_x, ensemble_size, len_t);
Ds = zeros(len_z, ensemble_size, len_t);


% generate initial ensemble
if strcmp(cov_or_sd, 'cov')
    eta = (randn(ensemble_size, len_x)*chol(P0))';
    
elseif strcmp(cov_or_sd, 'sd') 
    eta = (randn(ensemble_size, len_x)*(P0))';
end

X = bsxfun(@plus, x0(:), eta);
X0 = X;


% test for vectorization of odefun
ode_inc = @(xk) ode_incrementor(odefun, t(1), xk, t(2)-t(1), varargin{:});
try
    ode_inc(X);
    isvec=true;
catch ME
    if strcmp(ME.identifier, 'MATLAB:odearguments:SizeIC')
        cell2mat(arrayfun(@(i) ...
                 ode_inc(X(:, i)), 1:length(X), ...
                 'UniformOutput', false));    
        s = strcat(['Vectorized integration failed. ' ...
                   'Using arrayfun which is much slower. ' ...
                   'Consider modifying odefun so that it accepts a matrix of ' ...
                   'initial conditions (n states x n particles)' ...
                   'and using "myode" solver']);
        warning(s)
        isvec=false;
    else
        rethrow(ME)
    end
end
            
    

for i=1:length(t)

    tk = t(i);
    
    if i==1
        dt = t(2) - t(1);
    else
        dt = t(i) - t(i-1);
    end
    
    zk = z(i);
    
    % Add model noise from Q to the ensemble  
    if strcmp(cov_or_sd, 'cov')
        zeta = (randn(ensemble_size, len_x)*chol(Q))';
        
    elseif strcmp(cov_or_sd, 'sd') 
        zeta = (randn(ensemble_size, len_x)*(Q))';
    end

    X = bsxfun(@plus, X, zeta);

    
    if isvec
        if strcmp(state_eqn_type, 'cont-time')
            % propagate the particles through the system
            ode_inc = @(xk) ode_incrementor(odefun, tk, xk, dt, varargin{:});
            X_hat  = ode_inc(X);
            
        elseif strcmp(state_eqn_type, 'disc-time')
            ode_inc = @(xk) odefun(tk, xk);
            X_hat = ode_inc(X);
            
        end
    elseif isvec~=1 
        if strcmp(state_eqn_type, 'cont-time')
            ode_inc = @(xk) ode_incrementor(odefun, tk, xk, dt, varargin{:});
            X_hat = cell2mat(arrayfun(@(i) ...
                             ode_inc(X(:, i)), 1:length(X), ...
                             'UniformOutput', false)); 
        
        elseif strcmp(state_eqn_type, 'disc-time')                             
            ode_inc = @(xk) odefun(tk, xk);
            X_hat = cell2mat(arrayfun(@(i) ...
                        ode_inc(X(:, i)), 1:length(X), ...
                        'UniformOutput', false));
        end
    end
        
    if any(isnan(X_hat))
        error('Nan values in states')
    end
    
    X_hats(:, :, i) = X_hat;
    
    if strcmp(implementation, 'pert-obs')
        
        % make a matrix of perturbed observations with noise drawn from R
        zk_en = repmat(zk, 1, ensemble_size);  
        if strcmp(cov_or_sd, 'cov')
            epsilon = (randn(ensemble_size, len_z)*chol(R))'; 
        elseif strcmp(cov_or_sd, 'sd') 
            epsilon = (randn(ensemble_size, len_z)*(R))'; 
        end
        
        D = zk_en + epsilon;
        Ds(:, :, i) = D;
        
        % calculate the ensemble sample mean and covariance
        x_bar = mean(X_hat, 2);
        X_bar = bsxfun(@minus, X_hat, x_bar); % aka A
        C = cov(X_bar'); 

        % Calculate the gain matrix
        S = H*C*H'+ R;
        K = C*H'/S;
        correction = K*(D  - H*X_hat);
        X = X_hat + correction;
        Xs(:, :, i) = X; 
        
    end
      
    

end

Ts = t;

end

