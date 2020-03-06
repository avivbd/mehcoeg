function [ T, Y ] = ode_solver(odefun, t, y0, varargin )
%% ODE_SOLVER Function implementing a simple ODE solver.  
%

p = parse_inputs(odefun, t, y0, varargin{:});
[T, Y] = main(p.Results.odefun, p.Results.t,...
              p.Results.y0, p.Results.max_attempts,...
              p.Results.tol, varargin{:});


end

function p = parse_inputs(odefun, t, y0, varargin)
%% Parse Inputs
p = inputParser;
p.KeepUnmatched = true;

% required arguments
isFunctionHandle = @(x) isa(x, 'function_handle') & (nargin(x) == 2);
addRequired(p, 'odefun', isFunctionHandle)

is_valid_t = @(x) isnumeric(x) & isvector(x);
addRequired(p, 't', is_valid_t)
addRequired(p, 'y0', @isnumeric)

% optional parameters
default_max_attempts = 100;
validationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
addOptional(p, 'max_attempts', default_max_attempts, validationFcn)

default_tol = 1e-3;
addOptional(p, 'tol', default_tol, validationFcn)

parse(p, odefun, t, y0, varargin{:})


end

function [t, Y] = main(odefun, t, y0, max_attempts, tol, varargin)
%% Call the iterator

tstart = t(1);
tend = t(end);

% make initial dt guess according to 100 equally spaced data points
if length(t) == 2
    t = linspace(t(1), t(end), 100);
end
dt = t(2) - t(1);
assert(dt~=0, 'dt cannot be zero')

% over allocate arrays
Y = zeros(length(y0), 100000);
T = zeros(1, 100000);

T(:, 1) = t(1);
Y(:, 1) = y0;

yk = y0;
tk = tstart;
k = 1;

while true
    
    n_attempts = 0;

    while n_attempts <= max_attempts

        yk_plus_1 = ode_incrementor(odefun, tk, yk, dt, varargin{:});
    
        %variable timestep
        phi = abs(yk_plus_1 - yk)./abs(yk);
    
        if min(phi) > tol
            dt = dt/2;
            n_attempts = 1 + n_attempts;
            continue
        else
            tk = tk + dt;
            dt = dt*2;
            yk = yk_plus_1;
            k = k + 1;
            break
        end
        
    end    
    
if (n_attempts == max_attempts)
    error('max iterations reached without solution!')
end

Y(:, k) = yk;
T(1, k) = tk;

if dt>0
    if tend<tk 
        break
    end
    
elseif dt<0
    if tend>tk
        break
    end
end


end    

Y(:, k:end) = [];
T(:, k:end) = [];

if length(t)>2 
    YY = interp1(T', Y', t(:)', 'pchip');
    if size(t, 2)>1
        Y = YY';
    else
        Y = YY;
    end
end

end

