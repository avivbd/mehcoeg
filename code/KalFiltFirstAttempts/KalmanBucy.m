function KalmanBucy
%% Kalfilt 04_OCT_2015. Applies the Kalman filter to a simple state space
%model of the carbon cycle

clc
clear
close all
rng(1)


[t_interp,delC_interp] = get_data();

%run params
t_run_start = min(-t_interp)*1e6;
t_run_end = max(-t_interp)*1e6;

n_steps = length(t_interp);
Dt = (t_run_end - t_run_start)/n_steps;
t = linspace(t_run_start,t_run_end,n_steps);


u = delC_interp/10;


%initial mass
M0.P = 2e15;
M0.C = 3.8e18;


%isotpic parameters
p.delC_in = -5;
p.epsilon = 25;

%weathering fluxes
F0.in = 50e12;
p.forg = 0.2;
F0.volc = 5e12;
F0.volc_red = p.forg*F0.volc;
F0.volc_ox = (1 - p.forg)*F0.volc;
F0.worg = p.forg*(F0.in - F0.volc);
F0.wcarb = (1- p.forg)*(F0.in - F0.volc);
F0.wsil  = F0.volc_ox;
F0.wp = 36e9;
F0.bp = 36e9;

%isotopic compositions of weathering fluxes
p.delwcarb = p.delC_in + (F0.worg + F0.volc_red)/F0.in*p.epsilon;
p.delworg = p.delwcarb - p.epsilon;

%burial fluxes
F0.bcarb = (1- p.forg)*F0.in;
F0.borg = p.forg*F0.in;

%isotopic composition of carbonate burial flux
p.delC_0 = p.delC_in + p.forg*p.epsilon;
%     Fworg_0 = 14e12;

%flux coefficients
% k0.wo = F0.worg/M0.O;
k0.bc = F0.bcarb/M0.C;
k0.ws = F0.wsil/M0.C;
k0.bo = F0.borg/M0.P;
k0.wp = F0.wp/M0.C; 
k0.bp = F0.bp/M0.P;

%         kbo  kbc   kwp   kbp   
factor = [1     1     1     1   ];


%modify initial values by certain amounts set in factor above
k.bo = factor(1)*k0.bo;
k.bc = factor(2)*k0.bc;
k.ws = factor(2)*k0.ws;
k.wp = factor(3)*k0.wp;
k.bp = factor(4)*k0.bp;

%Coefficient matrix
%state matrix
A = [ -k.bp  k.wp; -k.bo -k.bc];

%control matrix
B = [0; F0.in];

%carbon isotopes observability matrix 
%linearization of delC = delC_in + f*epsilon
%f = Fbo/(Fbo + Fbc) = kbo*M.P/(kbo*M.P + kbc*M.C)



% ddelta_dMp = @(MP,MC) (p.epsilon*k.bo)/(MC*k.bc + MP*k.bo) ...
%     - (MP*p.epsilon*k.bo^2)/(MC*k.bc + MP*k.bo)^2;

% ddelta_dMc = @(MP,MC) -(MP*p.epsilon*k.bc*k.bo)/(MC*k.bc + MP*k.bo)^2;

%testing
% ddelta = [ddelta_dMp ddelta_dMc]*[M0.P;1.1*M0.C] + p.delC_0;

% C = [ddelta_dMp(M0.P,M0.C) ddelta_dMc(M0.P,M0.C)];

C = @(MP,MC) [(p.epsilon*k.bo)/(MC*k.bc + MP*k.bo) ...
    - (MP*p.epsilon*k.bo^2)/(MC*k.bc + MP*k.bo)^2, ...
       -(MP*p.epsilon*k.bc*k.bo)/(MC*k.bc + MP*k.bo)^2 ];

%with feedforward
%observability and feed-forward matrices. 
% C = eye(size(A));
D = 0 ;

%Covariances of the process noise 
mu_Q = [0 0];
% Sigma_Q = [   1e-2*(F0.wp)^2      0 ; 
%             0              1*(F0.in)^2];
Sigma_Q = [   1e-8*(F0.wp)^2      0 ; 
            0              1e-1*(F0.in)^2];


Q = chol(Sigma_Q);
%what percent of steady state flux does this variance represent?
percent_in_Q = [Q(1,1)/F0.wp*100,Q(2,2)/F0.in*100];

%Cov of measurement noise
% mu_R = [0 0];
% Sigma_R = [   1*(F0.wp)^2      0 ; 
%             0              1e2*(F0.in)^2];
mu_R = 0;
Sigma_R = (1e-4)^2; %std of plus minus 2 permil 
R = chol(Sigma_R);
%what percent of steady state flux does this variance represent?
% percent_in_R = [R(1,1)/F0.wp*100,R(2,2)/F0.in*100];


%Discretize the state matrix
A_disc = expm(A*Dt);
A_inv = eye(size(A))/A;

%Discretize input matrix
B_disc = (A_disc- eye(size(A_disc)))*A_inv*B;

%Discritize Q
QQ = [-A Q ; zeros(size(A)) A']*Dt;
QQQ = expm(QQ);
Q_disc = (QQQ(3:4,3:4))'*(QQQ(1:2,3:4));

%Discritize R
R_disc = R/Dt;

%make w: a random variable with an average of mu and sigma = Q_disc
% w_disc = repmat(mu_Q',1,n_steps) + Q_disc*randn(2,n_steps);

%make v: a random variable with an average of mu and sigma = R_disc
% v_disc = repmat(mu_R',1,n_steps) + R_disc*randn(length(R),n_steps);

%preallocate mem
% x_clean = zeros(2,n_steps);
% x = zeros(2,n_steps);
% y_clean = zeros(2,n_steps);
% y_proc = zeros(1,n_steps);
% y = zeros(1,n_steps);

%Initial conditions
%all zero


% for i = 1:n_steps-1
%             
%          %solution with no measurment or process noise
%          x_clean(:,i+1) = A_disc*x_clean(:,i) + B_disc*u(i);
%          y_clean(:,i+1) = C*x_clean(:,1);
%          
%          %solution with process noise
%          x(:,i+1) = A_disc*x(:,i) + B_disc*u(i) + w_disc(:,i);
%          y_proc(:,i+1) = C*x(:,i);
%          
%          %Add measurement noise to observables
%          y(:,i+1) = C*x(:,i) + v_disc(:,i);
%          
%          
%    
% end

%det PT data


%Kalman filter loop: operates on noise measurements z = y
z = delC_interp;

% Initial estimate for P, the uncertainty in the initial conditions
% where P is the a-posteriori estimate covariance matrix
% set to be the same as Q 
P = Q_disc;

%inital values estimates are set at zero
x_hat = zeros(2,1);

%preallocate
x_hat_kf = zeros(2,n_steps);
P_kf = zeros(2,2,n_steps);


for i = 1:n_steps    

    %predict
    [x_hat,P,y_hat] = KF_predict(x_hat,P,A_disc,Q_disc,B_disc,C,D,u(i));
    
    %update    
%     [x_hat,P] = KF_update(x_hat,P,z(:,i),C,R_disc);
    
    %test against skarka code
%     [x_hat,P] = kf_predict_skark(x_hat,P,A_disc,Q_disc,B_disc,u(i));
    
%     [x_hat,P,~,~,~,~] = kf_update_skark(x_hat,P,z(:,i),C,R_disc);

   %save the measurments
   x_hat_kf(:,i) = x_hat;      
   P_kf(:,:,i) = P;
   y_hat_kf(:,i) = y_hat;
   
end


%now run an RTS smoother on the estimates.  

% [x_hat_rts,~] = rts_smooth(x_hat_kf,P_kf,A_disc,Q_disc);

%compare to 
% [x_hat_rts,~,~] = rts_smooth_skark(x_hat_kf,P_kf,A_disc,Q_disc);

% [x_hat_rts,~] = urts_smooth1(x_hat_kf,...
%     reshape(cell2mat(P_kf),[2 2 100]),...
%     A_disc,Q_disc);


% Reconstruct solution
% M.P_clean = x_clean(1,:) + M0.P;
% M.C_clean = x_clean(2,:) + M0.C;
% M.P = x(1,:) + M0.P ;
% M.C = x(2,:) + M0.C;
% y_proc = y_proc + p.delC_0;
% y = y + p.delC_0;


M.P_hat_kf = x_hat_kf(1,:) + M0.P;
M.C_hat_kf = x_hat_kf(2,:) + M0.C;
delC_hat = y_hat_kf + p.delC_0;


% M.P_hat_rts = x_hat_rts(1,:) + M0.P;
% M.C_hat_rts = x_hat_rts(2,:) + M0.C;


m_plots = 2;
n_plots = 2;

hfig = figure;
set(hfig,'Position',[453     5   560   800])


%plot of forcing
subplot(m_plots,n_plots,1)

plot(t,u)
ylabel('Forcing')


%plot of model responses with no noise
% subplot(m_plots,n_plots,2)
% [hAx,hl1 hl2] = plotyy(t,M.P_clean,t,M.P);

% set(hl1,'LineStyle','-','Marker','none')
% set(hl2,'LineStyle','-','Marker','none')

% ylabel(hAx(1),'M_P no noise') % left y-axis
% ylabel(hAx(2),'M_P') % right y-axis

% set(hAx(1),'YLim',get(hAx(2),'YLim'),'YTickMode','auto')
% set(hAx(2),'YLim',[0 8e18],'YTickMode','auto')

%plot of model responses P with process noise (line)
%and measurement noise (dots)
% subplot(m_plots,n_plots,3)
% [hAx,hl1 hl2] = plotyy(t,M.C_clean,t,M.C);
% set(hAx(1),'YLim',get(hAx(2),'YLim'),'YTickMode','auto')
% ylabel(hAx(1),'M_C no noise') % left y-axis
% ylabel(hAx(2),'M_C') % right y-axis

% set(h(1),'LineStyle','-','Marker','none','Color','b')
% set(h(2),'LineStyle','none','Marker','o','MarkerSize',3,'Color','b')


%plot of model responses M with process noise (line)
%and measurement noise (dots)
% subplot(m_plots,n_plots,4)
% h = plot(t,y_proc,t,y);
% ylabel('M_C noisy') % left y-axis
% set(h(1),'LineStyle','-','Marker','none','Color',[0 0.5 0])
% set(h(2),'LineStyle','none','Marker','o','MarkerSize',3,'Color','k')


subplot(m_plots,n_plots,3)
h = plot(t,delC_hat);
ylabel('delC')


%plot kalman filter and RTS smoother results
subplot(m_plots,n_plots,2)
h = plot(t,M.P_hat_kf);
% set(h(1),'LineStyle','none','Marker','o','MarkerSize',3,'Color','b')
% set(h(2),'LineStyle','--','Marker','none','Color','b')
% set(h(3),'LineStyle','-','Marker','none','Color','b')
% ylabel('M_P filtered') % left y-axis
% ylim([0 4e15])




subplot(m_plots,n_plots,4)
h = plot(t,M.C_hat_kf);
% set(h(1),'LineStyle','none','Marker','o','MarkerSize',3,'Color',[0 0.5 0])
% set(h(2),'LineStyle','--','Marker','none','Color',[0 0.5 0])
% set(h(3),'LineStyle','-','Marker','none','Color',[0 0.5 0])
% ylabel('M_P filtered') % left y-axis





%%
%%Test against packaged kalman filters and rts smoothers



% line(t,M.P_hat,'Parent',hAx(1))
% line(t,M.C_hat,'Parent',hAx(2))

% assignin('base','M',M)
% assignin('base','y',y)

end

function [x,P,y] = KF_predict(x,P,A,Q,B,C,D,u)

    %prediction 
     x = A*x + B*u;
     %propegat error in prediction
     P = A * P * A' + Q;

     y  = C(x(1),x(2))*x + D*u;
%    y  = C*x + D*u;
end

function [x,P] = KF_update(x,P,z,C,R)

%update 
%    K = P*C'/(C*P*C' + R); %compute Kalman gain
%    x = x + K*(z - C*x); %correct values
%    P = P- K*C*P; %update cov
    K = P*C(x(1),x(2))'/(C(x(1),x(2))*P*C(x(1),x(2))' + R); %compute Kalman gain
   x = x + K*(z - C(x(1),x(2))*x); %correct values
   P = P- K*C(x(1),x(2))*P; %update cov



end

function [x,P] = rts_smooth(x,P,A,Q)

for j = length(x)-1:-1:1
    
    P_pred   = A * P(:,:,j) * A' + Q;
    
    D = P(:,:,j) * A' / P_pred;
    
    x(:,j) = x(:,j) + D*( x(:,j+1) - A*x(:,j));
    
    P(:,:,j) = P(:,:,j) + D*(P(:,:,j+1) - P_pred )*D';

    
    
%     D = p(:,:,j)*A'/p(:,:,j+1);
       
%       X(:,j) = x(:,j)     + D*(X(:,j+1) - A*x(:,j));
      
%       P(:,:,j) = p(:,:,j) + D*(P(:,:,j+1) - p(:,:,j+1))*D';

    
end

end


function [x,P] = kf_predict_skark(x,P,A,Q,B,u)

  %
  % Check arguments
  %
  if nargin < 3
    A = [];
  end
  if nargin < 4
    Q = [];
  end
  if nargin < 5
    B = [];
  end
  if nargin < 6
    u = [];
  end
  
  %
  % Apply defaults
  %
  if isempty(A)
    A = eye(size(x,1));
  end
  if isempty(Q)
    Q = zeros(size(x,1));
  end
  if isempty(B) & ~isempty(u)
    B = eye(size(x,1),size(u,1));
  end

  %
  % Perform prediction
  %
  if isempty(u)
    x = A * x;
    P = A * P * A' + Q;
  else
    x = A * x + B * u;
    P = A * P * A' + Q;
  end
  
end

function [X,P,K,IM,IS,LH] = kf_update_skark(X,P,y,H,R)

  %
  % Check which arguments are there
  %
  if nargin < 5
    error('Too few arguments');
  end

  %
  % update step
  %
  IM = H*X;
  IS = (R + H*P*H');
  K = P*H'/IS;
  X = X + K * (y-IM);
  P = P - K*IS*K';
  if nargout > 5
    LH = gauss_pdf(y,IM,IS);
  end
  
end

function [M,P,D] = rts_smooth_skark(M,P,A,Q)

  %
  % Check which arguments are there
  %
  if nargin < 4
    error('Too few arguments');
  end

  %
  % Extend A and Q if they are NxN matrices
  %
  if size(A,3)==1
    A = repmat(A,[1 1 size(M,2)]);
  end
  if size(Q,3)==1
    Q = repmat(Q,[1 1 size(M,2)]);
  end

  %
  % Run the smoother
  %
  D = zeros(size(M,1),size(M,1),size(M,2));
  for k=(size(M,2)-1):-1:1
    P_pred   = A(:,:,k) * P(:,:,k) * A(:,:,k)' + Q(:,:,k);
    D(:,:,k) = P(:,:,k) * A(:,:,k)' / P_pred;
    M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - A(:,:,k) * M(:,k));
    P(:,:,k) = P(:,:,k) + D(:,:,k) * (P(:,:,k+1) - P_pred) * D(:,:,k)';
  end
  
end


function [t_interp,delC_interp] = get_data()
%% Load data and make period groups and then Make a figure with all the data

%load previously saved data - put in same directory as this file 
load('/Users/avivbachan/My_Documents_Google_Drive/Research/PennStatePostdoc/Phanerozoic_d13C/Data/PhanDataSaltzmanThomas.mat')

% delC_full = delC_full;
% t_full = t_full;

%remove NaN datapoints
nanind = isnan(delC_full);
delC_full(nanind) = [];
t_full(nanind) = [];

%replace duplicates with their average
[t_unique,~,idx] = unique(t_full);
delC_unique = accumarray(idx,delC_full,[],@mean);


%assign groups for plot colors. Based on GTS 2012.
Quaternary      = [2.587,   0];
Neogene         = [23.03,   2.587];
Paleogene       = [66.14,   23.03];
Cretaceous      = [146.39,  66.14];
Jurassic		= [201.30,  146.39];
Triassic        = [252.16,  201.30]; 
Permian         = [298.88,  252.16];
Carboniferous   = [358.94,  298.88];
Devonian        = [419.20,  358.94];
Silurian        = [443.83,  419.20];
Ordovician      = [485.37,  443.83];
Cambrian        = [541,     485.37];
Ediacaran       = [635,     541];
Cryogenian      = [850,     635];
Tonian          = [1000,    850];

gr = zeros(size(t_full));

names = {
'Quaternary'
'Neogene'
'Paleogene'
'Cretaceous'
'Jurassic'
'Triassic'
'Permian'
'Carboniferous'
'Devonian'
'Silurian'
'Ordovician'
'Cambrian'
'Ediacaran'
'Cryogenian'
'Tonian'
};


for i = 1:length(t_full)
    
    if and(Quaternary(1)>=t_full(i),t_full(i)>Quaternary(2))
        gr(i) = 1;
    end
    
    if and(Neogene(1)>=t_full(i),t_full(i)>Neogene(2))
        gr(i) = 2;
    end
    
    if and(Paleogene(1)>=t_full(i),t_full(i)>Paleogene(2))
        gr(i) = 3;
    end 
    
    if and(Cretaceous(1)>=t_full(i),t_full(i)>Cretaceous(2))
        gr(i) = 4;
    end 
    
    if and(Jurassic(1)>=t_full(i),t_full(i)>Jurassic(2))
        gr(i) = 5;
    end 
    
    if and(Triassic(1)>=t_full(i),t_full(i)>Triassic(2))
        gr(i) = 6;
    end
    
    if and(Permian(1)>=t_full(i),t_full(i)>Permian(2))
        gr(i) = 7;
    end
    
    if and(Carboniferous(1)>=t_full(i),t_full(i)>Carboniferous(2))
        gr(i) = 8;
    end
    
    if and(Devonian(1)>=t_full(i),t_full(i)>Devonian(2))
        gr(i) = 9;
    end
    
    if and(Silurian(1)>=t_full(i),t_full(i)>Silurian(2))
        gr(i) = 10;
    end
    
    if and(Ordovician(1)>=t_full(i),t_full(i)>Ordovician(2))
        gr(i) = 11;
    end
    
    if and(Cambrian(1)>=t_full(i),t_full(i)>Cambrian(2))
        gr(i) = 12;
    end
    
    if and(Ediacaran(1)>=t_full(i),t_full(i)>Ediacaran(2))
        gr(i) = 13;
    end
    
    if and(Cryogenian(1)>=t_full(i),t_full(i)>Cryogenian(2))
        gr(i) = 14;
    end
    
    if and(Tonian(1)>=t_full(i),t_full(i)>Tonian(2))
        gr(i) = 15;
    end
    
end



%remove outliers
n_std = 2;%anypoint that is more then 2 permil from its neighbour
  
diff_delC = diff(delC_unique);

%get locations of outliers
ind_diff_delC_outliers = find(or(diff_delC>n_std,diff_delC<-n_std));

%make a datase that does not contain them
delC_no_outliers = delC_unique;
t_no_outliers = t_unique;
delC_no_outliers(ind_diff_delC_outliers) = [];
t_no_outliers(ind_diff_delC_outliers) = [];

%put a first order interpolation in the signal
t_interp = linspace(min(t_no_outliers),max(t_no_outliers),length(t_no_outliers));
delC_interp = interp1(t_no_outliers,delC_no_outliers,t_interp,'pchip');


end

function y = movingAverage(x, w)
   k = ones(1, w) / w;
   y = conv(x, k, 'same');
end

