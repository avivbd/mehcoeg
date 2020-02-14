function [ output_args ] = forward_model( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = params();
p.tspan = linspace(-253, -246, 100);

A1 = 0.75;
A2 = 0.25;
omega1 = 0.5;
omega1 = 1.3;
phi1 = pi/2;
phi2 = pi/3;
p.F_C_volc_t = p.F_C_volc_0*( 1 + A1*sin(omega1*p.tspan + phi1) + ...
                                  A2*sin(omega1*p.tspan + phi2)); 

M0 = [p.M_C , p.delC, p.M_ALK, p.M_Sr, p.R_Sr];
[T, Y] = ode45( @(t, x) ode_eqs(x, t, p), p.tspan, M0);

plot_ode(T,Y, p)


end

