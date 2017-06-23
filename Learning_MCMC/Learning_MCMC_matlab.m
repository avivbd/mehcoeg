function Learning_MCMC()

    clc
    clear
    close all
    % rng('default');  % For reproducibility
    
    
    
[t,Data,pl] = getdata();
   
% [] = getcostfun

MCMC(t,Data,pl,a,bounds_flag,cost_fun,fun,theta0,theta_bounds)

end





function [t,Data,pl] = getdata(flag)
%% Make synthetic data

switch flag
    case 1 %well shaped 1D targed
        
        
    
    case 2 %double well shaped 1D target
        
    case 3 %2D surface from 2 parameter timeseries
                
    case 4 %3D surface from 3 parameter timeseries
        
        n = 1000; %make time series
        t = linspace(0,4000,n); %time span of 0 to 1
        r = randn(1,length(t)); %gaussian random noise
        omega(1) = 1.5; %parameter value
        omega(2) = 1;
        omega(3) = 0.8;
        fun = @(omega) omega(3)*exp(-0.5*omega(2)*t*(2*pi/4000)).*(sin(2*pi*1*omega(1)*t*(2*pi/4000))+1)/400; %model
        xx = fun(omega); %process
        sigma = 0.02;
        Data = (xx + sigma^2*r); % process + noise
        
        pl.hfig1 = figure('Position',[916   502   533   303]);
        plot(t,Data)
        pl.hax1 = gca;
        
        a  = [1 1 1]; %sdt for guesses     
        theta0 = [0.1; 0.1; 0.1]; %initial guess for theta
        theta_bounds = [0.4 3; 0 3; -1 2]; %bounds on theta (only used if bounds flag = 'yes')
        sigma_d = .00025; %sigma of SSE
        bounds_flag = 'no';
    
        cost_fun = @(theta) (1/length(Data))*sum(((Data - fun(theta)).^2))/(sigma_d^2); %cost function
    
    
   
end
        

end





function MCMC(t,D,pl,a,bounds_flag,cost_fun,model,theta1,theta_bounds)


burntime = 5000;%toss out first however many steps
nsteps = 20000;%markov chain with n steps

ntheta = length(theta1);
nt = length(t);

chi1 = cost_fun(theta1); %get result for our initial parameters guess 

%preallocate
thetasave = zeros(ntheta,nsteps - burntime);%preallocate mem
chisave = zeros(ntheta,nsteps - burntime);

for i = 1:nsteps
	theta2  = theta1 + (a.*randn(1,length(theta1)))'; %make a second guess based on the initial one    	
	chi2 = cost_fun(theta2); %calculate the goodness of fit of the model output under the second guess
	Lratio = exp(-chi2 + chi1); %is the second guess better than the first one (Lratio>1) or worse (L<1)
    
    %compare to uniform random number on [0 1]
    switch bounds_flag
    case 'yes'
        c = sum(and(theta_bounds(:,1)<theta2 , theta2<theta_bounds(:,2)) ) == ntheta;  %check if proposed theta is out of bounds
        if  c*Lratio > rand(1)
            theta1 = theta2;
            chi1 = chi2;
        end
        
    case 'no'
        if  Lratio > rand(1)
            theta1 = theta2;
            chi1 = chi2;
        end
    end
    
    %save theta values
    thetasave(:,i)=theta1;
    chisave(:,i)=chi1;    

end

%remove burntime samples
ThetaNoBurn = thetasave(:,burntime:end);
ChiNoBurn = chisave(:,burntime:end);

hfig2 = figure('Position',[511    13   515   792]);
for i = 1:ntheta
    subplot(3,1,i)
    [P(i,:),xi(i,:)] = ksdensity(ThetaNoBurn(i,:));
    plot(xi(i,:),P(i,:))
%     xlim()
    xlabel(['theta ' num2str(i)])
    ylabel('P')
end

figure(pl.hfig1)
hold on
for i = 1:ntheta
[c,I] = max(P(i,:));
most_probable_theta(i) = xi(i,I);
end
plot(t,model(most_probable_theta),'-r','Parent',pl.hax1) %plot model
xlabel('time')
ylabel('Data')
legend('noisy data','model with best param value')

%test for normalacy of the residuals
residuals = model(most_probable_theta) - D;
figure
qqplot(residuals)


pl.hfig3 = figure('Position',[ 0     5   513   798]);
for i = 1:ntheta
    subplot(3,1,i)
    plot(1:nsteps,thetasave(i,:),'.')
    switch bounds_flag
    case 'yes'
        ylim(theta_bounds(i,:))
    end
%     plot([burntime burntime],[0 1],'r');  
    xlabel('MCMC iterations')
    ylabel('theta value')
%     legend('parameter choice','burntime')
end


end
	