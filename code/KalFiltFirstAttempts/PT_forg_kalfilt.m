function [  ] = PT_forg_kalfilt(  )
%PT_fborg_kalfilt Applies a kalman filter to forg estimates derived from
%the Early Triassic record

clear
clc
close all

plotdata = 'on'; %on/off


[t, delC, h] = plotfun(plotdata);

numSamples = length(t);
Q = 1e-5;
R = 0.0002;
hSig = dsp.SignalSource;
hSig.Signal = delC;
      
hTScope = dsp.TimeScope('NumInputPorts', 2, 'TimeSpan', numSamples, ...
          'TimeUnits', 'Seconds', 'YLimits',[-5 5], ...
          'ShowLegend', true); % Create the Time Scope
hKalman = dsp.KalmanFilter('ProcessNoiseCovariance', Q,...
          'MeasurementNoiseCovariance', R,...
          'InitialStateEstimate', 5,...
          'InitialErrorCovarianceEstimate', 1,...
          'ControlInputPort',false); %Create Kalman filter
      
i = 1;                      
while(~isDone(hSig))
%    trueVal = step(hSig);
%    noisyVal = trueVal + sqrt(R)*randn;
%    estVal = step(hKalman, noisyVal);
%    step(hTScope,noisyVal,trueVal,estVal);
   noisyVal = step(hSig);
%    noisyVal = trueVal + sqrt(R)*randn;
   estVal = step(hKalman, noisyVal);
   step(hTScope,noisyVal,estVal);
   estVals(i) =  estVal;
   i = i+1;
end         

plot(h.ax,t,estVals,'Color','r')

% %calculate forg from data
% forg = (Y-(-5))/25;
% 
% %calculate first differences
% difforg = diff(forg);
% %turn it into a "dynamical" system
% forg_test = zeros(size(forg));
% forg_test(1) = forg(1);
% 
% for i=1:length(difforg)
%     forg_test(i+1) = forg(i) + difforg(i);
% end

end

function [t_PT, delC_PT,h] = plotfun(plotdata)
%%load P-T data and prep it for regression

load('/Users/avivbachan/My_Documents_Google_Drive/Research/PennStatePostdoc/Phanerozoic_d13C/Data/Jon_Data_Including_Neoprot/delC_full.mat')
load('/Users/avivbachan/My_Documents_Google_Drive/Research/PennStatePostdoc/Phanerozoic_d13C/Data/Jon_Data_Including_Neoprot/t_full.mat')

%remove NaN datapoints
nanind = isnan(delC_full);
delC_full(nanind) = [];
t_full(nanind) = [];

%replace duplicates with their average
[t_unique,~,idx] = unique(t_full);
delC_unique = accumarray(idx,delC_full,[],@mean);

%remove outliers (point to point difference larger than plus minus 2 permil)
ind_diff_delC_outliers = find(or(diff(delC_unique)>2,diff(delC_unique)<-2));
delC_no_outliers = delC_unique;
t_no_outliers = t_unique;
delC_no_outliers(ind_diff_delC_outliers) = [];
t_no_outliers(ind_diff_delC_outliers) = [];

%restrict data to P-T interval
idx_PT = find(and(t_no_outliers>245,t_no_outliers<255));
t_PT = t_no_outliers(idx_PT);
delC_PT = delC_no_outliers(idx_PT);

t_PT = -1e6*t_PT;

%plot data
switch plotdata
    case 'on'
        scrsz = get(0,'ScreenSize');

        h.fig = figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2]);

        h.lines = plot(t_PT,delC_PT,'o');

        hold on

        h.ax = gca;

        xlabel('Time [Ma]')

        ylabel('\delta^{13}C  [permil]')

        legend off

        box on

        grid on

        ylim([-2 +8])

        set(h.fig,'PaperPositionMode','auto');

        %calculate forg from delbcarb = delin + f*eps
        pos = get(h.ax,'Position');
        Yticklabs = get(h.ax,'YTickLabel');
        h.ax2 = axes('Position',pos,'Parent',h.fig,'XTick',[],'YAxisLocation','right',...
            'Color','none', 'YTickLabel',num2str((str2num(Yticklabs)+5)/25));
        ylabel('f_{org}')
        uistack(h.ax2,'down')
        linkaxes([h.ax, h.ax2],'xy')
%         title('\delta^{13}C data of the past 1 Ga')
    otherwise
        h = [];
end

end