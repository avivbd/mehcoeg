%%

syms p_Fv p_d13C_volc p_Fwo p_delwo p_Fwcarb p_delC delC p_kbo Mp p_eps Mc p_kwp p_kbp p_kws

dMp = p_kwp*Mc - p_kbp*Mp;
dMc = p_Fv + p_Fwo - p_kws*Mc - p_kbo*Mp;
ddelta = ( p_Fv*(p_d13C_volc - delC) + p_Fwo*(p_delwo - delC) + ...
               p_Fwcarb*(p_delC - delC) + p_kbo*Mp*p_eps )/Mc;
           
           
s = [dMp, dMc, ddelta];

J = jacobian(s, [Mp; Mc; delC ])