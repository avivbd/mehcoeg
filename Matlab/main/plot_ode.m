function plot_ode(T,Y, p)

ylab = {'MC', 'delC', 'MALK','M_Sr', 'RSr'};

for i = 1:length(ylab)
    subplot(3,2,i)
    plot(T, Y(:, i))
    ylabel(ylab(i))
    grid on
    
end

subplot(3, 2, i+1)
plot(p.tspan, p.F_C_volc_t)
ylabel('Fvolc')

xlabel('time [Ma]')
grid on


end