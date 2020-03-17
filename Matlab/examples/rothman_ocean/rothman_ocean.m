function rothman_ocean()
%% Stiff continous linear system example
%  based on the Rothman Ocean example in Slingerland and Kump
%   
%   DIC            DOC
%  Reservoir1   Reservoir2 
% 
%   0.1
%    |
%    v     100
%  -----   --->   -----
% | 5e4 |        | 5e7 |
% |_____|  <---  |_____|
%    |      100
%    v
%   0.1
% 
% 

clc
clear
close all

p = params();

z = simulate_system_no_noise(p);
z_w = simulate_system_process_noise(p);
z_v = simulate_system_measurement_noise(p);
plot_results(z, z_w, z_v, p)

end

function  p = params()
%% get system parameters
p = struct();

p.A = [(-2e-6 -2e-3)  2e-6;...
         2e-3        -2e-6];
  
p.B = [0.2; 0];

p.C = eye(2);

p.x0 = [5e4;...
        5e7];

p.Q = [0.1 0;...
        0  0];
   
p.R = [0.1 0;...
       0   0];
   
p.t = linspace(0, 2e8, 100);

rng default
p.v = randn(0, sqrt(p.R));
p.w = sqrt(p.Q)*randn(p.n,1);

end

function z = simulate_system_no_noise(p)
    pp = p;
    pp.w = 0;
    pp.v = 0;
    z = simulate_system(pp);
end


function z_w = simulate_system_process_noise(p)
    pp = p;
    pp.v = 0;
    z_w = simulate_system(p);
end

function z_v = simulate_system_measurement_noise(p)
    z_v = simulate_system(p);
end

function z = simulate_system(p)

sol = ode15s(@(t, x) odefun(t, x, p),...
            [p.t(1) p.t(end)], p.x0);
x = deval(sol, p.t);
z = (p.C*x)';

end

function dx = odefun(t, x, p)
% u_t = interp1(p.t, p.u, t);
w_t = interp1(p.t, p.w, t);
dx = p.A*x + p.B + w_t; %*(u_t + w_t);
end

function plot_results(z, p)

subplot(211)
plot(p.t, z(:, 1),'--')
title('DIC reservoir')
xlabel('t (yr)'), ylabel('GtC')

subplot(212)
plot(p.t, z(:, 2),'--')
title('DOC reservoir')
xlabel('t (yr)'), ylabel('GtC')

end