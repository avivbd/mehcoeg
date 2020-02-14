function [pCO2, CO3] = c_sys(M_ALK, M_C, p)

DIC = M_C*1e21/p.Voc;
ALK = M_ALK*1e21/p.Voc;
alpha = 3.4e-5; % mM/ppm
k1 = 1e-6;
k2 = 1e-9;
k = k1/k2;
bar = (4*ALK + DIC*k - ALK*k)^2 + 4*(k-4)*ALK^2;
foo = ALK*k - DIC*k - 4*ALK + sqrt(bar);
H2CO3 = ((DIC - ALK) + foo/(2*(k-4)));
pCO2 = H2CO3/alpha;
HCO3 = H2CO3*k1;
CO3 = HCO3*k2;
end