function [k_south,r_south] = gamma_intersections_south(SA,CT,p,ocean,n)

% intersections_south  -  southerly isopycnal intersection statistics
%==========================================================================
%
% Usage:      
%  [k_south, r_south] = intersections_south(SA,CT,p,ocean,n)
%
% Method:
%  Sweep through data making southerly isopycnal surface calculations
%
% Input:
%  SA      - Absolute Salinity	         [ g/kg ]            {nz,ny,nx}
%  CT      - Conservative Temperature    [ deg C ]           {nz,ny,nx}
%  p       - pressure                    [ dbar ]            {nz,ny,nx}
%  ocean   -             {ny,nx}
%  n       -             {ny,nx}
%
% Output:
%  k_south - index of southerly intersection                 {nz,ny,nx}
%  r_south - interpolation distance from index               {nz,ny,nx}
%
% Author:     David Jackett
%
% Date:       06/04/06
%
%==========================================================================

% Define handles
handles.zbc = 0;
handles.quad = 0;

% initialize
[nz,ny,nx] = size(SA);
% nz2 = round(0.5*nz); %nz2 = 20;
k_south = nan(size(SA));
r_south = k_south;

% zonal bands ...
for j0 = 2:ny
    j0_south = j0 - 1;
    
    Idata = find(isfinite(SA(1,j0,:)));
    if ~isempty(Idata)
        nn = length(Idata);
        % meridional sweeps
        for i = 1:nn
            i0 = Idata(i);
            
            if isfinite(SA(1,j0_south,i0))
                % the central cast
                n0 = n(j0,i0);
                SA0 = SA(1:n0,j0,i0);
                CT0 = CT(1:n0,j0,i0);
                p0 = p(1:n0,j0,i0);
                
                % the southern cast
                n0_south = n(j0_south,i0);
                SA_south = SA(1:n0_south,j0_south,i0);
                CT_south = CT(1:n0_south,j0_south,i0);
                p_south = p(1:n0_south,j0_south,i0);
                
                % depth_ntp calculation
                SAdns_south = nan(n0,1);
                CTdns_south = SAdns_south;
                pdns_south = SAdns_south;
                
                for k0 = 1:n0
                    [SAdns_south(k0),CTdns_south(k0),pdns_south(k0)] = ...
                        gsw_depth_ntp(SA0(k0),CT0(k0),p0(k0),SA_south,CT_south,p_south);
                end
                
                tmp = pdns_south;
                pdns_south = nan(size(p));
                pdns_south(1:n0) = tmp;
                
                % intersection vitals
                for k0 = 1:n0
                    try
                        k_south(k0,j0,i0) = gamma_indx(p(1:n0_south,j0_south,i0),pdns_south(k0));
                        if pdns_south(k0) == p(n0_south,j0_south,i0)
                            k_south(k0,j0,i0) = n0_south - 1;
                        end
                        kk = k_south(k0,j0,i0);
                        r_south(k0,j0,i0) = (pdns_south(k0) - p(kk,j0_south,i0))/ ...
                            (p(kk+1,j0_south,i0) - p(kk,j0_south,i0));
                        
                        if handles.quad == 1
                            r_ntp = r_south(k0,j0,i0); %[k0,j0,i0,r_ntp]
                            pmid = 0.5*(p(kk+1,j0_south,i0) + p(kk,j0_south,i0));
                            [rho_u,alpha_u,beta_u] = gsw_rho_alpha_beta(SA(kk,j0_south,i0),CT(kk,j0_south,i0),pmid);
                            [rho_l,alpha_l,beta_l] = gsw_rho_alpha_beta(SA(kk+1,j0_south,i0),CT(kk+1,j0_south,i0),pmid);
                            term1 = (alpha_u + 0.5*r_ntp*(alpha_l - alpha_u))*(CT(kk+1,j0_south,i0) - CT(kk,j0_south,i0));
                            term2 = (beta_u + 0.5*r_ntp*(beta_l - beta_u))*(SA(kk+1,j0_south,i0) - SA(kk,j0_south,i0));
                            r_south(k0,j0,i0) = 0.5*r_ntp*(term2 - term1)*(rho_l + rho_u)/(rho_l - rho_u);
                            %r_quad = r_south(k0,j0,i0), dj_pause(0)
                        end
                    end
                    
                    % switch off between oceans
                    oce_test = ocean(j0,i0)*ocean(j0_south,i0);
                    if oce_test == 3 | oce_test == 5
                        k_south(k0,j0,i0) = -1;
                        r_south(k0,j0,i0) = nan;
                    end
                end
            end
        end
    end
    
    %     %			ints_south plot
    %     Iplot = find(isfinite(k_south));
    %     if ~isempty(Iplot)
    %         counter = counter + 1;
    %         if mod(counter,nwrite) == 0
    %             %             na = length(find(isfinite(r_south(:,1:j0,:))));
    %             %             ns = length(find(isfinite(s(:,1:j0,:))));
    %             %             ok = [j0,100*j0/ny,100*na/ns]
    %
    %             figure(3)
    %             subplot(2,2,3)
    %             zz = squeeze(k_south(nz2,:,:));
    %             dj_pltmp(long(1,:),lat(:,1),zz)
    %
    %             I_3 = find(zz==0);
    %             if ~isempty(I_3)
    %                 [j_triples,i_triples]= ind2sub([ny,nx],I_3);
    %                 hold on
    %                 plot(longss(i_triples),latss(j_triples),'m*')
    %                 hold off
    %             end
    %             title('k\_south')
    %             set(gcf,'Name','north/south intersections','NumberTitle','off','Color',[0.961 0.988 0.965])
    %
    %             subplot(2,2,4)
    %             dj_pltmp(long(1,:),lat(:,1),squeeze(r_south(nz2,:,:)))
    %             title('r\_south')
    %         end
    %     end
    
end


%			and save final intersections_south file

% figure(3)
% subplot(2,2,3)
% dj_pltmp(long(1,:),lat(:,1),squeeze(k_south(nz2,:,:)))
% title('k\_south')
%
% subplot(2,2,4)
% dj_pltmp(long(1,:),lat(:,1),squeeze(r_south(nz2,:,:)))
% title('r\_south')
%
% % triples
% inds_3 = find(k_south==0);
% if length(inds_3)>0
%     [k_triples,j_triples,i_triples]= ind2sub([nz,ny,nx],inds_3);
%     hold on
%     plot(longss(i_triples),latss(j_triples),'m*')
%     hold off
% end

end
