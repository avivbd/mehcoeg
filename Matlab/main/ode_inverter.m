function dy = ode_inverter(odefun, t, y, dt, varargin )
%   ODE_INVERTER Given a handle to an odefun of the form dy/dt = f(t, y)
%   calculates the dy in a way that is stable for stiff equations
%   The balance between implicit and explicit is set by Theta.
%   A theta=1 is fully implicit, so stable. 
%   A theta=0 is fully explicit, so more accurate. 

p = parse_args(odefun, t, y, dt, varargin{:});

dy = main(p.Results.odefun, p.Results.t, ...
         p.Results.y, p.Results.dt, ...
         p.Results.theta, varargin{:});

end

function p = parse_args(odefun, t, y, dt, varargin)
%% Parse Inputs
p = inputParser;
p.KeepUnmatched = true;

% required arguments
isFunctionHandle = @(x) isa(x, 'function_handle') & (nargin(x) == 2);

addRequired(p, 'odefun', isFunctionHandle)
addRequired(p, 't', @isnumeric)
addRequired(p, 'y', @isnumeric)
validationFcn = @(x) isnumeric(x) && isscalar(x);
addRequired(p, 'dt', validationFcn)

%balance between implicit and explicit solutions
defaultTheta = 0.5;
addOptional(p, 'theta', defaultTheta, validationFcn)
parse(p, odefun, t, y, dt, varargin{:})

end

function dy = main(odefun, t, y, dt, theta, varargin)
%% Calculate dy

% evaluate derivative
df = odefun(t, y);

%evaluate the jacobian
J = numjac(odefun, t, y, varargin{:});

%calculate change in function value (dy)
% see Equation 3.49 in Slingerland and Kump

if size(y, 2)>1
    [m,n,p] = size(J);
    shape_df = size(df);
    block_J = kron(speye(p), ones(m,n));
    block_J(logical(block_J)) = J(:);
    I  = eye(size(block_J));
    dy = (I/abs(dt) - theta*block_J)\df(:);
    dy = reshape(dy, shape_df);
    
else 
    I = eye(size(J));
    dy = (I/abs(dt) - theta*J)\df; 
    
end


end

