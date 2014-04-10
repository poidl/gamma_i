function [k_east,r_east] = gamma_intersections_east(SA,CT,p,ocean,n)

%==========================================================================
% intersections_east  -  easterly isopycnal intersection statistics
%
% Usage:
%  [k_east, r_east] = gamma_intersections_east(SA,CT,p,g,ocean,n,longs,lats)
%
% Description:
%  Sweep through data making easterly neutral tangent plane calculations
%
% Input:
%  SA      - Absolute Salinity	         [ g/kg ]            {nz,ny,nx}
%  CT      - Conservative Temperature    [ deg C ]           {nz,ny,nx}
%  p       - pressure                    [ dbar ]            {nz,ny,nx}
%  ocean   -             {ny,nx}
%  n       -             {ny,nx}
%
%  handles.zbc
%  handles.quad
%
% Output:
%  k_east  - index of easterly intersection                  {nz,ny,nx}
%  r_east  - interpolation distance from index               {nz,ny,nx}
%
% Author:
%  David Jackett   06/04/06
%
%==========================================================================

% Define handles
handles.zbc = 0;
handles.quad = 0;

% initialize
[nz,ny,nx] = size(SA);
k_east = nan(size(SA));
r_east = k_east;

% Calculate intersections in zonal bands
for j0 = 1:ny
    Idata = find(isfinite(SA(1,j0,:)));
    nn = length(Idata);
    
    if nn > 0
        %  meridional sweeps ...
        for i = 1:nn
            i0 = Idata(i);
            
            if i0 < nx
                i0_east = i0 + 1;
            else
                if handles.zbc == 1
                    i0_east = 1;
                else
                    i0_east = nan;
                end
            end
            
            if isfinite(i0_east) & isfinite(SA(1,j0,i0_east))
                % Set the central cast
                n0 = n(j0,i0);
                SA0 = SA(1:n0,j0,i0);
                CT0 = CT(1:n0,j0,i0);
                p0 = p(1:n0,j0,i0);
                
                % Set the eastern cast
                n0_east = n(j0,i0_east);
                SA_east = SA(1:n0_east,j0,i0_east);
                CT_east = CT(1:n0_east,j0,i0_east);
                p_east = p(1:n0_east,j0,i0_east);
                
                % Do the eastward depth_ns calculation
                SAdns_east = nan(n0,1);
                CTdns_east = SAdns_east;
                pdns_east = SAdns_east;
                
                for k0 = 1:n0
                    if ~isempty(SA_east)
%                          [SAdns_east(k0),CTdns_east(k0),pdns_east(k0)] = ...
%                              depth_ntp_iter(SA0(k0),CT0(k0),p0(k0),SA_east,CT_east,p_east);
                         [SAdns_east(k0),CTdns_east(k0),pdns_east(k0)] = ...
                             depth_ntp_orig(SA0(k0),CT0(k0),p0(k0),SA_east,CT_east,p_east);
                        
                    end
                end
                % Clean up the eastward depth_ns pressures, ensuring NaN's
                % were not over written with the depth_ns calculations.
                tmp = pdns_east;
                pdns_east = nan(nz,1);
                pdns_east(1:n0) = tmp;
                
                % intersection vitals
                for k0 = 1:n0
                    
                    k_east(k0,j0,i0) = indx_orig(p(1:n0_east,j0,i0_east),pdns_east(k0));
                    if pdns_east(k0) == p(n0_east,j0,i0_east)
                        k_east(k0,j0,i0) = n0_east - 1;
                    end
                    kk = k_east(k0,j0,i0);
                    r_east(k0,j0,i0) = (pdns_east(k0) - p(kk,j0,i0_east))/ ...
                                         (p(kk+1,j0,i0_east) - p(kk,j0,i0_east));
                    
%                     if handles.quad == 1
%                         r_ntp = r_east(k0,j0,i0);
%                         pmid = 0.5*(p(kk+1,j0,i0_east) + p(kk,j0,i0_east));
%                         [rho_u,alpha_u,beta_u] = gsw_rho_alpha_beta(SA(kk,j0,i0_east),CT(kk,j0,i0_east),pmid);
%                         [rho_l,alpha_l,beta_l] = gsw_rho_alpha_beta(SA(kk+1,j0,i0_east),CT(kk+1,j0,i0_east),pmid);
%                         part1 = (alpha_u + 0.5*r_ntp*(alpha_l - alpha_u))*(CT(kk+1,j0,i0_east) - CT(kk,j0,i0_east));
%                         part2 = (beta_u + 0.5*r_ntp*(beta_l - beta_u))*(SA(kk+1,j0,i0_east) - SA(kk,j0,i0_east));
%                         r_east(k0,j0,i0) = 0.5*r_ntp*(part2 - part1)*(rho_l + rho_u)/(rho_l - rho_u);
%                     end   
                    
                    
                    
                    oce_test = ocean(j0,i0)*ocean(j0,i0_east);                  
                    if oce_test == 3 | oce_test==5
                        k_east(k0,j0,i0) = -1;
                        r_east(k0,j0,i0) = nan;
                    end
                    
%                     % Stop the talking between the Pacific Ocean and the
%                     % Gulf of Mexico (across Panama), and the Pacific Ocean
%                     % and the North Indian Ocean (across Asia).
%                     oce_test = ocean(j0,i0)*ocean(j0,i0_east);
%                     if oce_test == 5
%                         k_east(k0,j0,i0) = -1;
%                         r_east(k0,j0,i0) = nan;
%                     end
%                     % Restrict the talking though the Indonesian Throughflow
%                     % Allow it to occur only when the depth is less than 1200 dbar.
%                     if oce_test == 32
%                         [I] = find(pdns_east < 1200);
%                         if ~isempty(I)
%                             k_east(I,j0,i0) = -1;
%                             r_east(I,j0,i0) = nan;
%                         end
%                     end
                end
            end
        end
    end
end

end
