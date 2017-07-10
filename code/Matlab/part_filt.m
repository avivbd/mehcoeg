function [S, E, p] = part_filt(y, p)

muM = p.M; %augemented state vector
Z_cov = diag(p.z.^2); %create measurement noise covariance matrix
wts = 1/p.n_ensemble*ones(1, p.n_ensemble);

mus = zeros(length(muM), length(y));
ensembles =  zeros(length(muM), p.n_ensemble, length(y));

for k=1:(length(y))

    W = repmat(p.w, 1, p.n_ensemble) .* randn(length(muM), p.n_ensemble);            %create process noise
    V = repmat(p.z, 1, p.n_ensemble) .* randn(size(y,1), p.n_ensemble);     %create measurement noise
    M = repmat(muM, 1, p.n_ensemble);     %state ensemble
    
    M = M + p.dt*p.odeFun(p.dt*k, M, p) + W;     %Update states. 
    Y = p.obsFun(M, p) + V; %make ensemble of predictions. 
    
    y_obs = repmat(y(:,k),1, p.n_ensemble) + V;  %make ensemble of observations
    dif = y_obs - Y;
    probs = mvnpdf(dif', 0, Z_cov);
    
    wts = probs'.*wts;
    wts = wts/sum(wts);
    muM = sum(repmat(wts, length(muM),1).*M, 2);
    M = (datasample(M', length(M), 'Weights', wts'))';
    wts = 1/p.n_ensemble*ones(1, p.n_ensemble);
    sigma = cov(M');
    M = (mvnrnd(muM, sigma, p.n_ensemble))';
    M = M+M*eps;
        
    mus(:,k) = muM; 
    ensembles(:,:,k) = M;

end

S.Mp_filt = mus(1,:)';
S.Mc_filt = mus(2,:)';
S.delC_filt = mus(3,:)';
S.F_forcing_filt = mus(4,:)';


E.Mp =  squeeze(ensembles(1,:,:))';
E.Mc = squeeze(ensembles(2,:,:))';
E.delC = squeeze(ensembles(3,:,:))';
E.F_forcing = squeeze(ensembles(4,:,:))';



end