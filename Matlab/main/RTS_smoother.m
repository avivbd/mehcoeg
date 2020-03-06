function [X, P] = RTS_smoother(xks_post, Pks_post, xks_pri, Pks_pri, Fks)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

lsim = size(xks_pri, 2);
nx = size(xks_pri, 1);

X  = zeros(nx, lsim);
P  = zeros(nx,nx,lsim);

X(:,lsim)   = xks_post(:, lsim);
P(:,:,lsim) = Pks_post(:, :, lsim);


for it = lsim-1:-1:1
    C         = Pks_post(:,:,it) * Fks(:, :, it+1)' / Pks_pri(:,:,it+1);
    X(:,it)   = xks_post(:,it)   + C * ( X(:,it+1)   - xks_pri(:,it+1)  );
    P(:,:,it) = Pks_post(:,:,it) + C * ( P(:,:,it+1) - Pks_pri(:,:,it+1)) * C';
    
end

end

