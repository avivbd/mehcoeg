function J = numjac(odefun, t, y, varargin )
%NUMJAC Calculate the Jacobian of a function numerically
%   Calculate the matrix of derivatives of a set of differential
%   equations with respect to the state variables.
%   The equations should be held in a function of the form dy/dt = f(t, y). 
%   Have the option of three different methods: 'complex-step', 
%   'first-order-fd', 'second-order-fd'. Also have a choice for the error
%   tolerance on the finite difference  (fd) calculation.

%%
p = parse_inputs(odefun, t, y, varargin{:});
J = main(p.Results.odefun, ...
         p.Results.t, ...
         p.Results.y, ...
         p.Results.eps, ...
         p.Results.method);

end

function p = parse_inputs(odefun, t, y, varargin)
%% Parse Inputs
p = inputParser;
p.KeepUnmatched = true;

% required arguments
isFunctionHandle = @(x) isa(x, 'function_handle');
addRequired(p, 'odefun', isFunctionHandle)

addRequired(p, 't', @isnumeric)
addRequired(p, 'y', @isnumeric)

% options for fd method
validOptions = {'complex-step', 'first-order-fd', 'second-order-fd'};
checkMethod = @(x) any(validatestring(x, validOptions));
defultMethod = 'complex-step';
addParameter(p, 'method', defultMethod, checkMethod)

%how sensitve is the jacobian
defaultEps = 1e-8;
validationFcn = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addParameter(p, 'eps', defaultEps, validationFcn)


parse(p, odefun, t, y, varargin{:})

end

function J = main(odefun, t, y, eps, method)
%% Calculate Jacobian by purturbing the inputs one at a time

num_y = size(y, 2);

if num_y>1
    single_y = y(:, 1);
    m = numel(odefun(t, single_y ));
    n = numel(single_y );
    J = zeros(m, n, num_y );
else
    m = numel(odefun(t, y));
    n = numel(y);
    J = zeros(m,n);
end

for ii = 1:size(y, 1) %perturb one value of y at a time
    if num_y>1
        h = y(ii, :)*eps;
    else
        h = y(ii)*eps; % get diff relative to y_i
    end
    
    h(h==0) = eps;

    switch method
        case 'first-order-fd'
            df_y = odefun(t, y); 
            y_plus_h = y;        
            y_plus_h(ii) = y(ii) + h;
            df_y_plus_h = odefun(t, y_plus_h);
            J(:,ii)  = (df_y_plus_h - df_y)/h; 

        case 'second-order-fd'
            y_plus_h = y;        
            y_plus_h(ii) = y(ii) + h;
            df_y_plus_h = odefun(t, y_plus_h);

            y_minus_h = y;
            y_minus_h(ii) = y(ii) - h;
            df_y_minus_h = odefun(t, y_minus_h);
            J(:, ii) = (df_y_plus_h - df_y_minus_h)/(2*h);

        case 'complex-step'
            if num_y>1
                y_plus_ih = y;        
                y_plus_ih(ii, :) = y(ii, :) + 1i*h;
                df_y_plus_ih = odefun(t, y_plus_ih);
                imag_df_y_plus_ih = imag(df_y_plus_ih);
                J(:, ii, :) = bsxfun(@rdivide, imag_df_y_plus_ih, h); 
            else
                y_plus_ih = y;        
                y_plus_ih(ii) = y(ii) + 1i*h;
                df_y_plus_ih = odefun(t, y_plus_ih);
                imag_df_y_plus_ih = imag(df_y_plus_ih);
                J(:, ii) = imag_df_y_plus_ih/h;     
            end
            
            
    end
end
    
end