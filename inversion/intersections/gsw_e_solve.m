function [SAns,CTns,pns,iter] = gsw_e_solve(SA,CT,p,e,k,SA0,CT0,p0)

% gsw_e_solve                                 find the zero of the function
%==========================================================================
% 
% USAGE:         
%  [SAns,CTns,pns,iter] = gsw_e_solve(SA,CT,p,e,k,SA0,CT0,p0)
%
% DESCRIPTION:
%  Find the zero of the e function using a bisection method, designed for
%  neutral tangent plane
%   
% INPUT:
%  SA	=	array of cast Absolute Salinities                      [ g/kg ]
%  CT	=	array of cast Conservative Temperatures               [ deg C ]
%  p	=	array of cast pressures                                  [ db ]
%  e	=	array of cast e values
%  k	=	interval (k-1, k) contains the zero  
%  SA0  =	the bottle Absolute Salinity                           [ g/kg ]
%  CT0  =	the bottle Conservative Temperature                   [ deg C ]
%  p0   =	the bottle pressure                                      [ db ]
%
% OUTPUT:	
%  SAns  =  Absolute Salinity of e zero                            [ g/kg ]
%  CTns  =	Conservative Temperature of e zero                    [ deg C ]
%  pns	 =  pressure of e zero                                       [ db ]       
%
% AUTHOR: 
%  DRJ on 17/06/03
%  Revised by Guillaume Sérazin 
%  
%==========================================================================

kl = k(1);
ku = k(2);
pl = p(kl);
el = e(1);
pu = p(ku);
eu = e(2);

iter = 0;
success = 0;

while success == 0
    iter = iter + 1;
    pm = 0.5*(pl + pu);
    [SAm, CTm] = gsw_stp_interp([SA(kl),SA(ku)],[CT(kl),CT(ku)],[p(kl),p(ku)],pm);
%--------------------------------------------------------------------------        
% gsw_stp_interp
%--------------------------------------------------------------------------        
% n = length(p);
% if p(1) < p0 && p0 < p(n)
%     I_p = find(p >= p0);
%     k = I_p(1) - 1;
% elseif p0 <= p(1)
%     k = 1;
% elseif p0 >= p(n)
%     k = n - 1;
% else
%     SAm = NaN;
%     CTm = NaN;
%     return
% end
% r = (p0 - p(k))/(p(k+1) - p(k));
% SAm = SA(k) + r*(SA(k+1) - SA(k));
% CTm = CT(k) + r*(CT(k+1) - CT(k));
%--------------------------------------------------------------------------    

%     [sigl, sigu] = gsw_sig_vals(SA0,CT0,p0,SAm,CTm,pm);
    pmid = 0.5*(p0 + pm);
%     sigl = gsw_rho(SA0,CT0,pmid) - 1000;
%     sigu = gsw_rho(SAm,CTm,pmid) - 1000;        
%     em = sigu - sigl;
    em = gsw_rho(SAm,CTm,pmid) - gsw_rho(SA0,CT0,pmid);
    if el*em < 0
        pu = pm;
        eu = em;
    elseif em*eu < 0
        pl = pm;
        el = em;
    elseif em == 0
        SAns = SAm;
        CTns = CTm;
        pns = pm;
        success = 1;
    end
    
    if success == 0
        if (abs(em) <= 5e-5) && (abs(pu-pl) <= 5e-3 && pm > 0)
            SAns = SAm;
            CTns = CTm;
            pns = pm;
            success = 1;
        elseif iter <= 10
            success = 0;
        else
%             disp('WARNING in e-solve')
%             disp(['iter: ', int2str(iter), '  em: ', num2str(em), '  dp: ', num2str(abs(pu-pl))])
            SAns = NaN;
            CTns = NaN;
            pns = NaN;
            success = 0;
            return
        end
    end
end

end