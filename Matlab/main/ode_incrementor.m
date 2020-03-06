function yk_plus_1 = ode_incrementor(odefun, tk, yk, dt, varargin)
% ODE_ITERATOR Increments the solution of an ode dt/dt = f(x, t)
% forward in time one step.

p = parse_args(odefun, tk, yk, dt, varargin{:});
yk_plus_1 = main(p.Results.odefun, p.Results.tk, ...
                 p.Results.yk, p.Results.dt, ...
                 p.Results.solver, varargin{:});

end

function p = parse_args(odefun, tk, yk, dt, varargin)
%% Parse Inputs
p = inputParser;
p.KeepUnmatched = true;

% required arguments
isFunctionHandle = @(x) isa(x, 'function_handle') & (nargin(x) == 2);
addRequired(p, 'odefun', isFunctionHandle)

addRequired(p, 'tk', @isnumeric)
addRequired(p, 'yk', @isnumeric)

validationFcn = @(x) isnumeric(x) && isscalar(x);
addRequired(p, 'dt', validationFcn)

valid_solvers = {'ode15s', 'myode'};
isValidSolver = @(solver) any(strcmp(solver, valid_solvers));

defaultSolver = 'myode';
addParameter(p, 'solver', defaultSolver, isValidSolver)

parse(p, odefun, tk, yk, dt, varargin{:})

end

function yk_plus_1 = main(odefun, tk, yk, dt, solver, varargin)
%% return new value of state
    
    % get dy in a stable way
    if strcmp(solver, 'ode15s')
        sol = ode15s(odefun, [tk, tk+dt], yk);
        yk_plus_1 = sol.y(:, end);
        
    elseif strcmp(solver, 'myode')
        dy = ode_inverter(odefun, tk, yk, dt, varargin{:});
        yk_plus_1 = yk + dy; %add to previous value of y
    end
    
end





