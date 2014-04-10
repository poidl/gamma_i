function [SPns,ptns,pns] = depth_ntp_pt25(SP0,pt0,p0,SP,pt,p)

%  Find the position where the neutral tangent plane passing through a bottle
%  intersects a neighbouring cast
%
%  Usage :        [SPns,ptns,pns] = depth_ntp(SP0,pt0,p0,SP,pt,p)
%
%  Input :        SP0    the bottle salinity
%                 pt0    the bottle Conservative Temperature
%                 p0    the bottle pressure
%                 SP     vector of cast salinities
%                 pt     vector of cast Conservative Temperatures
%                 p     vector of cast pressures
%
%  Output :       SPns   salinity of the ntp intersection
%                 ptns   Conservative Temperature of the intersection
%                 pns   pressure of the intersection
%
%  Units :        SP	  g/kg (TEOS-10)
%                 Conservative Temperatures   degrees C (IPS-90)
%                 pressures      db
%  DRJ on 17/06/03


n = length(SP);
e = zeros(n,1);

%		find the bottle pairs containing a crossing

ncr = 0;
for k = 1:n
    [sigl,sigu] = sig_vals_pt25(SP0,pt0,p0,SP(k),pt(k),p(k));
    e(k) = sigu - sigl;
    if k > 1
        if e(k-1) == 0                     %  an exact crossing at the k-1 bottle
            ncr = ncr + 1;
            SPns = SP(k-1);
            ptns = pt(k-1);
            pns = p(k-1);
        elseif e(k)*e(k-1) < 0            %  a crossing between k-1 and k bottles
            ncr = ncr+1;                        %  some Newton-Raphson iterations
            pc0 = p(k-1) - e(k-1)*(p(k) - p(k-1))/(e(k) - e(k-1));
            iter = 0;
            success = 0;
            while success == 0
                iter = iter + 1;
                [SPc0,ptc0] = stp_interp([SP(k-1),SP(k)],[pt(k-1),pt(k)],[p(k-1),p(k)],pc0);
                [sigl,sigu] = sig_vals_pt25(SP0,pt0,p0,SPc0,ptc0,pc0);
%                  pmid = 0.5*(p0 + pc0);
%                  sig1 = gsw_rho_pt(SP0,pt0,pmid) - 1000;
%                  sig2 = gsw_rho_pt(SPc0,ptc0,pmid) - 1000;
                ec0 = sigu - sigl;
                p1 = 0.5*(p(k-1) + pc0);
                ez1 = (e(k-1) - ec0)/(pc0 - p(k-1));
                p2 = 0.5*(pc0 + p(k));
                ez2 = (ec0 - e(k))/(p(k) - pc0);
                r = (pc0 - p1)/(p2 - p1);
                ecz_0 = ez1 + r*(ez2 - ez1);
                if iter == 1
                    ecz0 = ecz_0;
                else
                    ecz0 = -(ec0 - ec_0)/(pc0 - pc_0);
                    if ecz0 == 0 
                        ecz0 = ecz_0;
                    end
                end
                pc1 = pc0 + ec0/ecz0;
                %  strategy when iteration jumps out of inteval
                if (pc1 < p(k-1)) | (pc1 > p(k))
                    SP = SP(:);
                    pt = pt(:);
                    p = p(:);
                    data = SP.*pt.*p;
                    indsp = find(isfinite(data));
                    [SPns,ptns,pns,niter] = e_solve(SP(indsp),pt(indsp),p(indsp),e(indsp),k,SP0,pt0,p0);
                    if pns < p(k-1) | pns > p(k)
                        disp('****ERROR**** in depth_ntp')
                        return
                    else
                        success = 1;
                    end
                else
                    %  test accuracy of the iterate
                    eps = abs(pc1-pc0);
                    if abs(ec0) <= 5e-5 & eps <= 5e-3
                        SPns = SPc0; 
                        ptns = ptc0; 
                        pns = pc0;
                        success = 1; 
                        niter = iter;
                    elseif iter > 10
                        [SPns,ptns,pns,niter] = e_solve(SP,pt,p,e,k,SP0,pt0,p0);
                        success = 1;
                    else
                        pc_0 = pc0; 
                        ec_0 = ec0; 
                        pc0 = pc1;
                        success = 0;
                    end
                end
            end
        end
    end
    if k == n & e(k) == 0                  %  the last bottle
        ncr = ncr + 1;
        SPns = SP(k); 
        ptns = pt(k); 
        pns = p(k);
    end
end

%  multiple and no crossings

if ncr == 0
    if e(1) > 0                                 %  outcropping
        SPns = -99.1;
        ptns = -99.1;
        pns = -99.1;
    else                                      %  undercropping
        SPns = -99.2;
        ptns = -99.2;
        pns = -99.2;
    end
elseif ncr >= 2                               %  multiple crossings
    SPns = -99.3;
    ptns = -99.3;
    pns = -99.3;
end


return