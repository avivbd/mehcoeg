function [] = vis_and_get_data()
%% Gets data for model.

%PT_INV Top level function for calling routines that carry out inverse
%modeling of the PT carbon cycle perturbation and early Triassic warmth

Data = load_data();

f = plot_raw_data(Data);

[TS, TI] = smooth_and_interp_data(Data);

plot_smoothed_data(TS, f)


end

function [TS, TI] = smooth_and_interp_data(Data)
%% split into arrays
Age_d13C = Data.d13_data(:, 1);
d13C = Data.d13_data(:, 2);

Age_RSr = Data.sr_data(:, 1);
sr = Data.sr_data(:, 2);

Age_Temp = Data.temp_data(:, 1);
temp = Data.temp_data(:, 2);

%% make smoothed curves
d13C_fit = smooth_d13C(Age_d13C, d13C);
temp_fit = smooth_temp(Age_Temp, temp);
sr_fit = smoooth_sr(Age_RSr, sr);

t_min = max([min(Age_d13C), min(Age_RSr), min(Age_Temp)]);
t_max = min([max(Age_d13C), max(Age_RSr), max(Age_Temp)]);

t_interp = linspace(t_min, t_max, 150);
t_interp = reshape(t_interp, numel(t_interp), 1);

d13C = feval(d13C_fit, t_interp);
temp = feval(temp_fit, t_interp);
sr = feval(sr_fit, t_interp);

TS = table(t_interp, d13C, sr, temp);

TI = get_deduped_data(Data, t_interp);

end

function TI = get_deduped_data(Data, t_interp)
%% make deduped data
d13C_deduped = avg_dups(Data.d13_data);
sr_deduped = avg_dups(Data.sr_data);
temp_deduped = avg_dups(Data.temp_data);

d13C = interp1(d13C_deduped(:, 1), d13C_deduped(:, 2), t_interp, 'pchip');
sr = interp1(sr_deduped(:, 1), sr_deduped(:, 2), t_interp, 'pchip');
temp = interp1(temp_deduped(:, 1), temp_deduped(:, 2), t_interp, 'pchip');

TI = table(t_interp, d13C, sr, temp);

end

function fitresult = smooth_d13C(Age_d13C, d13C)

[xData, yData] = prepareCurveData( Age_d13C, d13C );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
excludedPoints = excludedata( xData, yData, 'Indices', [142 348 357 363 372 380 422] );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.975127070508847;
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData, excludedPoints, 'predobs', 0.9 );
% legend( h, 'd13C vs. Age_d13C', 'Excluded d13C vs. Age_d13C', 'untitled fit 1', 'Lower bounds (untitled fit 1)', 'Upper bounds (untitled fit 1)', 'Location', 'NorthEast' );
% % Label axes
% xlabel Age_d13C
% ylabel d13C
% grid on

end

function fitresult = smooth_temp(Age_Temp, temp)

[xData, yData] = prepareCurveData( Age_Temp, temp );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
excludedPoints = excludedata( xData, yData, 'Indices', [51 76 80 110] );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.990703596121366;
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData, excludedPoints, 'predobs', 0.9 );
% legend( h, 'temp vs. Age_Temp', 'Excluded temp vs. Age_Temp', 'untitled fit 1', 'Lower bounds (untitled fit 1)', 'Upper bounds (untitled fit 1)', 'Location', 'NorthEast' );
% % Label axes
% xlabel Age_Temp
% ylabel temp
% grid on

end

function fitresult = smoooth_sr(Age_RSr, sr )

[xData, yData] = prepareCurveData( Age_RSr, sr );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
excludedPoints = excludedata( xData, yData, 'Indices', [51 72 74 81 96 99 109 114] );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.762925825619236;
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData, excludedPoints, 'predobs', 0.9 );
% legend( h, 'sr vs. Age_RSr', 'Excluded sr vs. Age_RSr', 'untitled fit 1', 'Lower bounds (untitled fit 1)', 'Upper bounds (untitled fit 1)', 'Location', 'NorthEast' );
% % Label axes
% xlabel Age_RSr
% ylabel sr
% grid on

end


function f = plot_raw_data(Data)

ax1 = subplot(3,1,1);
plot(Data.d13_data(:, 1), Data.d13_data(:, 2), 'o')
ax1.XLim = [-254 -244];

ax2 = subplot(3,1,2);
plot(Data.sr_data(:, 1), Data.sr_data(:, 2), 'o')
ax2.XLim = [-254 -244];

ax3 = subplot(3,1,3);
plot(Data.temp_data(:, 1), Data.temp_data(:, 2), 'o')
ax3.XLim = [-254 -244];

f = gcf();

end

function plot_smoothed_data(TS, f)

ax3 = f.Children(3);
hold( ax3, 'on')
plot(TS.t_interp, TS.d13C, 'color', 'k', 'linewidth', 2, 'Parent', ax3)


ax2 = f.Children(2);
hold( ax2, 'on')
plot(TS.t_interp, TS.sr, 'color', 'k', 'linewidth', 2, 'Parent', ax2)


ax1 = f.Children(1);
hold( ax1, 'on')
plot(TS.t_interp, TS.temp, 'color', 'k', 'linewidth', 2, 'Parent', ax1)


end

function deduped_mat = avg_dups(array)

[C,ia,idx] = unique(array(:,1),'stable');
val = accumarray(idx,array(:,2),[],@mean); 
deduped_mat = [C val];

end




