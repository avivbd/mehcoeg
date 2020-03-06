function plot_enkf(t, z, X0, Xs, X_hats, Ds, state_names, H, z_marker)

plot_results(t, z, X0, X_hats, Ds, Xs, state_names, H, z_marker)


end


function plot_results(t, z, X0, X_hats, Ds, Xs, state_names, H, z_marker)
%%

figure();
n_ensemble = size(X0, 2);
n_states = size(X0, 1);

for k=1:n_states
    
    subplot(n_states, 1, k);
    hold on

    ensemble = reshape(squeeze(Xs(k, :, :)), 1, []);
    t_jigger = repmat(t, n_ensemble, 1) + ...
                    mean(diff(t))/5*randn(n_ensemble, length(t));
    
    t_jigger = t_jigger(:);
    scatter(t_jigger, ensemble', 5, ...
           'Marker', 'o', ...
           'MarkerFaceAlpha', 0.25,...
           'MarkerEdgeAlpha', 0.25)
     grid on  
    
     
       
   if H(k)~=0
       plot(t, z, z_marker)
   end

   ylabel(state_names{k})
   
end




end


function interactive_filter_progression()

%%
figure();

jigger = 0.001*randn(size(X0, 2), 1);
t0 = p.t(1)-p.dt;

ax1 = subplot(2, 1, 1);
ax2 = subplot(2, 1, 2);
% ax3 = subplot(4, 1, 3);
% ax4 = subplot(4, 1, 4);
hold(ax1, 'on')
hold(ax2, 'on')
% hold(ax3, 'on')
% hold(ax4, 'on')
plot(t0+jigger, X0(1, :), '.k', 'Parent', ax1);
plot(t0+jigger, X0(2, :), '.k', 'Parent', ax2);
% plot(t0+jigger, X0(3, :), '.k', 'Parent', ax3);
% plot(t0+jigger, X0(4, :), '.k', 'Parent', ax4);


for i=1:length(p.t)
    
    hold on
    plot(p.t(i)+jigger, X_hats(1, :, i), '.b', 'Parent', ax1, 'MarkerSize', 15)
    plot(p.t(i), mean(X_hats(1, :, i)), 'ob', 'Parent', ax1, 'MarkerSize', 15)
    
    plot(p.t(i)+jigger, X_hats(2, :, i), '.b', 'Parent', ax2, 'MarkerSize', 15)
    plot(p.t(i), mean(X_hats(2, :, i)), 'ob', 'Parent', ax2, 'MarkerSize', 15)
    
%     plot(p.t(i)+jigger, X_hats(3, :, i), '.b', 'Parent', ax3, 'MarkerSize', 15)
%     plot(p.t(i), mean(X_hats(3, :, i)), 'ob', 'Parent', ax3, 'MarkerSize', 15)
    
%     plot(p.t(i)+jigger, X_hats(4, :, i), '.b', 'Parent', ax4, 'MarkerSize', 15)
%     plot(p.t(i), mean(X_hats(4, :, i)), 'ob', 'Parent', ax4, 'MarkerSize', 15)
    
    waitforbuttonpress;
    plot(p.t(i)+jigger, Ds(1, :, i), 'or', 'Parent', ax2)
    plot(p.t(i), mean(Ds(1, :, i)), 'or', 'Parent', ax2, 'MarkerSize', 15)
    waitforbuttonpress;
    
    plot(p.t(i)+jigger, Xs(1, :, i), '.k', 'Parent', ax1, 'MarkerSize', 15);
    plot(p.t(i), mean(Xs(1, :, i)), 'ok', 'Parent', ax1, 'MarkerSize', 15);
    
    plot(p.t(i)+jigger, Xs(2, :, i), '.k', 'Parent', ax2, 'MarkerSize', 15);
    plot(p.t(i), mean(Xs(2, :, i)), 'ok', 'Parent', ax2, 'MarkerSize', 15);
    
%     plot(p.t(i)+jigger, Xs(3, :, i), '.k', 'Parent', ax3, 'MarkerSize', 15);
%     plot(p.t(i), mean(Xs(3, :, i)), 'ok', 'Parent', ax3, 'MarkerSize', 15);
    
%     plot(p.t(i)+jigger, Xs(4, :, i), '.k', 'Parent', ax4, 'MarkerSize', 15);
%     plot(p.t(i), mean(Xs(4, :, i)), 'ok', 'Parent', ax4, 'MarkerSize', 15);
    
    
    waitforbuttonpress;
%     set(h_line, 'XData', [], 'YData', []);
    
end
ylabel(p.state_names{k})



end
