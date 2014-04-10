function [k_east,r_east] = intersections_east_orig(s,t,p,ocean,n)

%%%   intersections_east  -  easterly isopycnal intersection statistics
%%%
%%%   Usage:         [k_east,r_east] = intersections_east(s,t,p,g,ocean,n,longs,lats)
%%%
%%%   Input:         s       - salinity	                                {nz,ny,nx}
%%%                  t       - temperature (ITS90)                      {nz,ny,nx}
%%%                  p       - pressure (db)                            {nz,ny,nx}
%%%                  longs   - longitude (degE)                         {ny,nx}
%%%                  lats    - latitude (degN)                          {ny,nx}
%%%
%%%   Output:        k_east  - index of easterly intersection           {nz,ny,nx}
%%%                  r_east  - interpolation distance from index        {nz,ny,nx}
%%%
%%%   Method:        Sweep through data making easterly neutral tangent plane calculations
%%%
%%%   Author:        David Jackett
%%%
%%%   Date:       06/04/06
%%%

%%			initialize
[nz,ny,nx] = size(s);
%nz2 = round(nz/2); %nz2 = 40;

k_east = nan*ones(size(s));
r_east = k_east;

%% zonal bands 
for j0 = 1:ny
    
    inds = find(isfinite(s(1,j0,:)));
    nn = length(inds);
    if nn>0
        
        %% and meridional sweeps 
        for i = 1:nn
            
            i0 = inds(i);
            
            if i0<nx
                i0_east = i0 + 1;
            else
                i0_east = nan;
            end
            
            if isfinite(i0_east) & isfinite(s(1,j0,i0_east))
                
                %%			the central cast
                n0 = n(j0,i0);
                s0 = s(1:n0,j0,i0);
                t0 = t(1:n0,j0,i0);
                p0 = p(1:n0,j0,i0);
                
                %%          the eastern cast
                n0_east = n(j0,i0_east);
                s_east = s(1:n0_east,j0,i0_east);
                t_east = t(1:n0_east,j0,i0_east);
                p_east = p(1:n0_east,j0,i0_east);
                
                %%          depth_ns calculation
                sdns_east = nan*ones(n0,1);
                tdns_east = sdns_east;
                pdns_east = sdns_east;
                
                for k0 = 1:n0
%                    [sdns_east(k0),tdns_east(k0),pdns_east(k0)] = ...
%                        depth_ntp_orig(s0(k0),t0(k0),p0(k0),s_east,t_east,p_east); 
                    [sdns_east(k0),tdns_east(k0),pdns_east(k0)] = ...
                             depth_ntp_iter(s0(k0),t0(k0),p0(k0),s_east,t_east,p_east);
                end
                % Clean up the eastward depth_ns pressures, ensuring NaN's
                % were not over written with the depth_ns calculations.
                tmp = pdns_east;
                pdns_east = nan(nz,1);
                pdns_east(1:n0) = tmp;
                
                %%          intersection vitals
                for k0 = 1:n0
                    if ~isnan(pdns_east(k0))
                        
                    if abs(pdns_east(k0) + 99) <= 1e-2
                        k_east(k0,j0,i0) = nan;
                        r_east(k0,j0,i0) = nan;
                    elseif abs(pdns_east(k0) + 99.1) <= 1e-2
                        k_east(k0,j0,i0) = 1;
                        r_east(k0,j0,i0) = nan;
                    elseif abs(pdns_east(k0) + 99.2) <= 1e-2
                        k_east(k0,j0,i0) = n0;
                        r_east(k0,j0,i0) = nan;
                    elseif abs(pdns_east(k0) + 99.3) <= 1e-2
                        % not_ok = 'a triple'
                        k_east(k0,j0,i0) = 0;
                        r_east(k0,j0,i0) = nan;
                    else
                        
                    %    k_east(k0,j0,i0) = indx_orig(p(1:n0_east,j0,i0_east),pdns_east(k0));
                        k_east(k0,j0,i0) = gamma_indx(p(1:n0_east,j0,i0_east),pdns_east(k0));
                        
                        if pdns_east(k0) == p(n0_east,j0,i0_east)
                            k_east(k0,j0,i0) = n0_east - 1;
                        end
                        kk = k_east(k0,j0,i0);
                        
                        r_east(k0,j0,i0) = (pdns_east(k0) - p(kk,j0,i0_east))/  ...
                                            (p(kk+1,j0,i0_east) - p(kk,j0,i0_east));
                        
                        %                         if handles.quad==1
                        %
                        %                             r_ntp = r_east(k0,j0,i0); %[k0,j0,i0,r_ntp]
                        %
                        %                             pmid = (p(kk+1,j0,i0_east)+p(kk,j0,i0_east))/2;
                        %
                        %                             ct_u = gsw_CT_from_t(s(kk,j0,i0_east),t(kk,j0,i0_east),p(kk,j0,i0_east));
                        %                             [rho_u,alfa_u,beta_u] = gsw_rho_alpha_beta(s(kk,j0,i0_east),ct_u,p(kk,j0,i0_east));
                        %  %                           rho_local_u = gsw_rho(s(kk,j0,i0_east),ct_u,pmid);
                        %
                        %                             ct_l = gsw_CT_from_t(s(kk+1,j0,i0_east),t(kk+1,j0,i0_east),p(kk+1,j0,i0_east));
                        %                             [rho_l,alfa_l,beta_l] = gsw_rho_alpha_beta(s(kk+1,j0,i0_east),ct_l,p(kk+1,j0,i0_east));
                        % %                            rho_local_l = gsw_rho(s(kk+1,j0,i0_east),ct_l,pmid);
                        %
                        %                             term1 = (alfa_u + r_ntp*(alfa_l - alfa_u)/2)*(ct_l - ct_u)/(rho_l - rho_u);
                        %                             term2 = (beta_u + r_ntp*(beta_l - beta_u)/2)*(s(kk+1,j0,i0_east) - s(kk,j0,i0_east))/(rho_l - rho_u);
                        %                             r_east(k0,j0,i0) = r_ntp*(rho_l + rho_u)*(-term1 + term2)/2; %r_quad = r_east(k0,j0,i0), dj_pause(0)
                        %
                        %                         end
                        
                    end
                    
                    %%       switch off between oceans
                    oce_test = ocean(j0,i0)*ocean(j0,i0_east);                  
                    if oce_test == 3 | oce_test==5
                        k_east(k0,j0,i0) = -1;
                        r_east(k0,j0,i0) = nan;
                    end
                    end
                end
            end
        end
        
    end
    
end


return


