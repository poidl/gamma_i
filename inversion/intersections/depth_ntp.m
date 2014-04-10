function [SAns,CTns,pns] = depth_ntp(SA0,CT0,p0,SA,CT,p)

%  Find the position where the neutral tangent plane passing through a bottle
%  intersects a neighbouring cast
%
%  Usage :        [SAns,CTns,pns] = depth_ntp(SA0,CT0,p0,SA,CT,p)
%
%  Input :        SA0    the bottle salinity
%                 CT0    the bottle Conservative Temperature
%                 p0    the bottle pressure
%                 SA     vector of cast salinities
%                 CT     vector of cast Conservative Temperatures
%                 p     vector of cast pressures
%
%  Output :       SAns   salinity of the ntp intersection
%                 CTns   Conservative Temperature of the intersection
%                 pns   pressure of the intersection
%
%  Units :        SA	  g/kg (TEOS-10)
%                 Conservative Temperatures   degrees C (IPS-90)
%                 pressures      db
%  DRJ on 17/06/03


n = length(SA);
e = zeros(n,1);

%		find the bottle pairs containing a crossing

ncr = 0;
for k = 1:n
    [sigl,sigu] = sig_vals(SA0,CT0,p0,SA(k),CT(k),p(k));
    e(k) = sigu - sigl;
    if k > 1
        if e(k-1) == 0                     %  an exact crossing at the k-1 bottle
            ncr = ncr + 1;
            SAns = SA(k-1);
            CTns = CT(k-1);
            pns = p(k-1);
        elseif e(k)*e(k-1) < 0            %  a crossing between k-1 and k bottles
            ncr = ncr+1;                        %  some Newton-Raphson iterations
            pc0 = p(k-1) - e(k-1)*(p(k) - p(k-1))/(e(k) - e(k-1));
            iter = 0;
            success = 0;
            while success == 0
                iter = iter + 1;
                [SAc0,CTc0] = stp_interp([SA(k-1),SA(k)],[CT(k-1),CT(k)],[p(k-1),p(k)],pc0);
                [sigl,sigu] = sig_vals(SA0,CT0,p0,SAc0,CTc0,pc0);
%                  pmid = 0.5*(p0 + pc0);
%                  sig1 = gsw_rho_CT(SA0,CT0,pmid) - 1000;
%                  sig2 = gsw_rho_CT(SAc0,CTc0,pmid) - 1000;
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
                    SA = SA(:);
                    CT = CT(:);
                    p = p(:);
                    data = SA.*CT.*p;
                    indsp = find(isfinite(data));
                    [SAns,CTns,pns,niter] = e_solve(SA(indsp),CT(indsp),p(indsp),e(indsp),k,SA0,CT0,p0);
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
                        SAns = SAc0; 
                        CTns = CTc0; 
                        pns = pc0;
                        success = 1; 
                        niter = iter;
                    elseif iter > 10
                        [SAns,CTns,pns,niter] = e_solve(SA,CT,p,e,k,SA0,CT0,p0);
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
        SAns = SA(k); 
        CTns = CT(k); 
        pns = p(k);
    end
end

%  multiple and no crossings

if ncr == 0
    if e(1) > 0                                 %  outcropping
        SAns = -99.1;
        CTns = -99.1;
        pns = -99.1;
    else                                      %  undercropping
        SAns = -99.2;
        CTns = -99.2;
        pns = -99.2;
    end
elseif ncr >= 2                               %  multiple crossings
    SAns = -99.3;
    CTns = -99.3;
    pns = -99.3;
end


return