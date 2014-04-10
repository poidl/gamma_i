function [k_north, r_north] = gamma_intersections_north(SA,CT,p,ocean,n,long,lat,handles)

% intersections_north  -  northerly isopycnal intersection statistics
%
% Usage:
%  [k_north,r_north] = gamma_intersections_north(SA,CT,p,g,long,lat,handles)
%
% Method:
%  Sweep through data making northerly isopycnal surface calculations
%
% Input:
%  SA     - Absolute Salinity	[ g/kg ]			        {nz,ny,nx}
%  CT     - Conservative Temperature (ITS90)                {nz,ny,nx}
%  p      - pressure (dbar)                                 {nz,ny,nx}
%  long  - longitude (degE)                                 {ny,nx}
%  lat   - latitude (degN)                                  {ny,nx}
%
% Output:
%  k_north - index of northerly intersection                {nz,ny,nx}
%  r_north - interpolation distance from index              {nz,ny,nx}
%
% Author:
%  David Jackett  6/4/06
%


counter = 0;
nwrite = 10;
check_calcs = 0;

longss = nanmean(long(:,:));
latss = nanmean(lat(:,:)');

% initialize
[nz,ny,nx] = size(SA);
nz2 = round(nz/2); %nz2 = 20;
k_north = nans(size(SA));
r_north = k_north;

% zonal bands
for j0 = 1:ny-1
    j0_north = j0 + 1;
    
    Idata = find(isfinite(SA(1,j0,:)));
    if ~isempty(Idata)
        nn = length(Idata);
        % meridional sweeps
        for i = 1:nn
            i0 = inds(i);
            long0 = long(j0,i0);
            lat0 = lat(j0,i0);
            
            if isfinite(SA(1,j0_north,i0))
                % the central cast
                n0 = n(j0,i0);
                SA0 = SA(1:n0,j0,i0);
                CT0 = CT(1:n0,j0,i0);
                p0 = p(1:n0,j0,i0);
                
%                 if check_calcs == 1
%                     location = [longss(i0),latss(j0)]
%                     CENTRAL_CAST = [s0,t0,p0]
%                 end
                
                % the northern cast
                n0_north = n(j0_north,i0);
                SA_north = SA(1:n0_north,j0_north,i0);
                CT_north = CT(1:n0_north,j0_north,i0);
                p_north = p(1:n0_north,j0_north,i0);
                
%                 if check_calcs == 1
%                     NORTHERN_CAST = [SA_north,CT_north,p_north]
%                 end
                
                %  depth_ntp calculation
                SAdns_north = nans(n0,1);
                CTdns_north = SAdns_north;
                pdns_north = SAdns_north;
                                
                for k0 = 1:n0
                    [SAdns_north(k0),CTdns_north(k0),pdns_north(k0)] = ...
                        gsw_depth_ntp(SA0(k0),CT0(k0),p0(k0),SA_north,CT_north,p_north);
                end
                
                tmp = pdns_north;
                pdns_north = nans(nz,1);
                pdns_north(1:n0) = tmp;
                
%                 if check_calcs == 1
%                     disp('NORTHERN CAST DEPTH_NS COMPUTATIONS')
%                     cast_vals = [pdns_north(:)]
%                 end
                
                % intersection vitals
                for k0 = 1:n0
%                     if abs(pdns_north(k0)+99) <= 1e-2
%                         k_north(k0,j0,i0) = nan;
%                         r_north(k0,j0,i0) = nan;
%                     elseif abs(pdns_north(k0)+99.1) <= 1e-2
%                         k_north(k0,j0,i0) = 1;
%                         r_north(k0,j0,i0) = nan;
%                     elseif abs(pdns_north(k0)+99.2) <= 1e-2
%                         k_north(k0,j0,i0) = n0;
%                         r_north(k0,j0,i0) = nan;
%                     elseif abs(pdns_north(k0)+99.3) <= 1e-2
%                         % not_ok = 'a triple'
%                         k_north(k0,j0,i0) = 0;
%                         r_north(k0,j0,i0) = nan;
%                     else
                        try
                            k_north(k0,j0,i0) = indx(p(1:n0_north,j0_north,i0),pdns_north(k0));
                            if pdns_north(k0) == p(n0_north,j0_north,i0)
                                k_north(k0,j0,i0) = n0_north - 1;
                            end
                            kk = k_north(k0,j0,i0);
                            r_north(k0,j0,i0) = (pdns_north(k0) - p(kk,j0_north,i0))/...
                                (p(kk+1,j0_north,i0) - p(kk,j0_north,i0));
                            
                            if handles.quad == 1
                                r_ntp = r_north(k0,j0,i0); %[k0,j0,i0,r_ntp]
                                pmid = 0.5*(p(kk+1,j0_north,i0) + p(kk,j0_north,i0));
                                [rho_u,alpha_u,beta_u] = gsw_rho_alpha_beta(SA(kk,j0_north,i0),CT(kk,j0_north,i0),pmid);                               
                                [rho_l,alpha_l,beta_l] = gsw_rho_alpha_beta(SA(kk+1,j0_north,i0),CT(kk+1,j0_north,i0),pmid);
                                term1 = (alpha_u + 0.5*r_ntp*(alpha_l - alpha_u))*(CT(kk+1,j0_north,i0) - CT(kk,j0_north,i0));
                                term2 = (beta_u + 0.5*r_ntp*(beta_l - beta_u))*(SA(kk+1,j0_north,i0) - SA(kk,j0_north,i0));
                                r_north(k0,j0,i0) = 0.5*r_ntp*(term2 - term1)*(rho_l + rho_u)/(rho_l - rho_u);
                                %r_quad = r_north(k0,j0,i0);
                            end
                        end
                  %  end
                    % switch off between oceans
                    oce_test = ocean(j0,i0)*ocean(j0_north,i0);
                    if oce_test == 3 | oce_test == 5
                        k_north(k0,j0,i0) = -1;
                        r_north(k0,j0,i0) = nan;
                    end
                end
            end
        end
    end
    
    % ints_north plot
%     Iplot = find(isfinite(k_north));
%     if ~isempty(Iplot)
%         counter = counter + 1;
%         if mod(counter,nwrite)==0
%             
%             figure(2)
%             subplot(2,2,1)
%             zz = squeeze(k_north(nz2,:,:));
%             dj_pltmp(long(1,:),lat(:,1),zz)
%             
%             I_3 = find(zz == 0);
%             if ~isempty(I_3)
%                 [j_triples,i_triples]= ind2sub([ny,nx],I_3);
%                 hold on
%                 plot(longss(i_triples),latss(j_triples),'m*')
%                 hold off
%             end
%             title('k\_north')
%             set(gcf,'Name','north/south intersections','NumberTitle','off','Color',[0.961 0.988 0.965])
%             
%             % figure(2)
%             subplot(2,2,2)
%             dj_pltmp(long(1,:),lat(:,1),squeeze(r_north(nz2,:,:)))
%             title('r\_north')
%             %             figure(gcf);
%             %             dj_pause(1);
%         end
%     end
    
end

% and save final intersections_north file

figure(2)
subplot(2,2,1)
dj_pltmp(long(1,:),lat(:,1),squeeze(k_north(nz2,:,:)))
title('k\_north')

subplot(2,2,2)
dj_pltmp(long(1,:),lat(:,1),squeeze(r_north(nz2,:,:)))
title('r\_north')

% triples
I_3 = find(k_north==0);
if ~isempty(I_3)
    [k_triples,j_triples,i_triples]= ind2sub([nz,ny,nx],I_3);
    hold on
    plot(longss(i_triples),latss(j_triples),'m*')
    hold off
end

% figure(gcf);
% dj_toc;
% dj_pause(1)

end
