function [SAns,CTns,pns,dSAns,dCTns,dpns] = gsw_neutral_surfaces0(SA,CT,p,gamma,glevels)

%
% [SAns,CTns,pns,dSAns,dCTns,dpns] = gsw_neutral_surfaces0(SA,CT,p,gamma,glevels)
%
% Description
%  For a cast of hydrographic data (SA,CT,p) that has been pre-labelled with gamma,
%  find the Absolute Salinities, Conservative Temperatures and pressures of neutral
%  density surfaces specified by glevels.
%
% Input        
%  SA          Absolute Salinity of bottles
%  CT          Conservative Temperature of bottles
%  p           pressure of bottles
%  gamma       gamma values of bottles
%  glevels     neutral density surfaces to be outputed
%
% Output 
%  SAns        Absolute Salinity on the neutral density surfaces
%  CTns        Conservative Temperature on the neutral density surfaces
%  pns         pressure of the neutral density surfaces
%  dSAns       Absolute Salinity errors
%  dCTns       Conservative Temperature errors
%  dpns        pressure errors
%
% Note that NaN's indicate under or outcropping.
%
% DRJ on 16/06/03


%n2 = 2;
ptol = 1e-3;
n = length(SA);
ng = length(glevels);

SAns = nan(size(glevels));
CTns = SAns;
pns = SAns;
dSAns = SAns;
dCTns = SAns;
dpns = SAns;

% detect error condition & adjust cast to avoid an exact crossing
in_error = 0;
for k = 1:n
    if gamma(k) < -10
        in_error = 1;
    end
    for ig = 1:ng
        if gamma(k) == glevels(ig) 
            gamma(k) = gamma(k) + 1e-10;
        end
    end
end

if in_error == 1
    %gamma,glevels
    %error('neutral_surfaces0: gamma value is too small');
    return
end

% loop over the surfaces
ierr = 0;

for ig = 1:ng
    % find the intervals of intersection
    nint = 0;
    intz = nan(n,1);
    for k = 1:(n-1)
        gmin = min(gamma(k),gamma(k+1));
        gmax = max(gamma(k),gamma(k+1));
        if gmin <= glevels(ig) & glevels(ig) <= gmax
            nint = nint + 1;
            intz(nint) = k;
        end
    end
    intz(nint+1:n,1) = []; 
    
    % find point(s) of intersection
    if nint == 0
        SAns(ig) = nan; 
        CTns(ig) = nan; 
        pns(ig) = nan;
        dSAns(ig) = 0; 
        dCTns(ig) = 0; 
        dpns(ig) = 0;
    else  
        % choose the central interval
        if mod(nint,2) == 0 & intz(1) >= 0.5*n
            int_middle = 0.5*(nint + 2);
        else
            int_middle = floor(0.5*(nint + 1));
        end
        
        % loop over all intersections
        for i_int = 1:nint
            k = intz(i_int);
            
            % coefficients of a quadratic for gamma
            [rho_l, alpha_l, beta_l] = gsw_rho_alpha_beta(SA(k),CT(k),p(k));
            
            [rho_u, alpha_u, beta_u] = gsw_rho_alpha_beta(SA(k+1),CT(k+1),p(k+1));

            [rho_mid, alpha_mid, beta_mid] = gsw_rho_alpha_beta(0.5*(SA(k) + SA(k+1)),0.5*(CT(k) + CT(k+1)),0.5*(p(k) + p(k+1)));

            delSA = SA(k+1) - SA(k); 
            delCT = CT(k+1) - CT(k);          
            delp = p(k+1) - p(k); 
            delp2 = delp*delp;
            
            bden = rho_mid*(beta_mid*delSA - alpha_mid*delCT);          
            if abs(bden) <= 1e-6
                bden = 1e-6;
            end
            
            b_mid = (gamma(k+1) - gamma(k))/bden;
            
            % coefficients           
            a_part = delSA*(beta_u - beta_l) - delCT*(alpha_u - alpha_l);
            a = (a_part*b_mid*rho_mid)/(2*delp2);
            
            b_part = delSA*(p(k+1)*beta_l - p(k)*beta_u) - delCT*(p(k+1)*alpha_l - p(k)*alpha_u);
            b = (b_part*b_mid*rho_mid)/delp2;
            
            c_part1 = delSA*(beta_l*(p(k) - 2*p(k+1)) + beta_u*p(k)) -  ...
                       delCT*(alpha_l*(p(k) - 2*p(k+1)) + alpha_u*p(k));
            c_part2 = gamma(k) + (b_mid*rho_mid*p(k)*c_part1)/(2*delp2);
            c = c_part2 - glevels(ig);
            
            delta =  b*b - 4*a*c;
            
            % solve the quadratic
            if (a ~= 0) & (bden ~= 1e-6) & (delta >= 0) % | isnan(a)
                q = -0.5*(b + sign(b)*sqrt(delta));
                pns1 = q/a; 
                pns2 = c/q;
                if (pns1 >= p(k) - ptol) & (pns1 <= p(k+1) + ptol)
                    pns(ig) = min(p(k+1),max(pns1,p(k)));
                elseif (pns2 >= p(k) - ptol) & (pns2 <= p(k+1) + ptol)
                    pns(ig) = min(p(k+1),max(pns2,p(k)));
                else
                    rg = (glevels(ig) - gamma(k))/(gamma(k+1) - gamma(k));
                    pns(ig) = p(k) + rg*(p(k+1) - p(k));
                end
            else
                rg = (glevels(ig) - gamma(k))/(gamma(k+1) - gamma(k));
                pns(ig) = p(k) + rg*(p(k+1) - p(k));
            end
            
            [SAns(ig),CTns(ig)] = gsw_stp_interp(SA,CT,p,pns(ig));
            
            % write multiple values to file
            if nint > 1
                % find median values and errors
                if i_int == 1
                    SAns_top = SAns(ig); 
                    CTns_top = CTns(ig); 
                    pns_top = pns(ig);
                end
                
                if i_int == int_middle
                    SAns_middle = SAns(ig);
                    CTns_middle = CTns(ig); 
                    pns_middle = pns(ig);
                end
                
                if i_int == nint
                    SAns(ig) = SAns_middle;
                    CTns(ig) = CTns_middle;
                    pns(ig) = pns_middle;
                    if nargout > 3
                        if (pns_middle - pns_top) > (pns(ig) - pns_middle)
                            dSAns(ig) = SAns_middle - SAns_top;
                            dCTns(ig) = CTns_middle - CTns_top;
                            dpns(ig) = pns_middle - pns_top;
                        else
                            dSAns(ig) = SAns(ig) - SAns_middle;
                            dCTns(ig) = CTns(ig) - CTns_middle;
                            dpns(ig) = pns(ig) - pns_middle;
                        end
                    end
                end
                
            else
                if nargout > 3
                    dSAns(ig) = 0;
                    dCTns(ig) = 0;
                    dpns(ig) = 0;
                end
            end
        end
    end
end

SAns = SAns(:); 
CTns = CTns(:); 
pns = pns(:);
dSAns = dSAns(:); 
dCTns = dCTns(:); 
dpns = dpns(:);

end