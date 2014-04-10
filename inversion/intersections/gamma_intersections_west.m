function [k_west,r_west] = gamma_intersections_west(SA,CT,p,ocean,n)

%   intersections_west  -  westerly isopycnal intersection statistics
%==========================================================================
%
% Usage:      
%  [k_west,r_west] = gamma_intersections_west(SA,CT,p,ocean,n)
%
% Method: 
%  westerly isopycnal surface calculations
%
% Input:
%  SA      - Absolute Salinity	         [ g/kg ]            {nz,ny,nx}
%  CT      - Conservative Temperature    [ deg C ]           {nz,ny,nx}
%  p       - pressure                    [ dbar ]            {nz,ny,nx}
%  ocean   -             {ny,nx}
%  n       -             {ny,nx}
%
% Output:
%  k_west - index of westerly intersection      {nz,ny,nx}
%  r_west - interpolation distance from index   {nz,ny,nx}
%
%   Author:     David Jackett
%
%   Date:       06/04/06
%
%==========================================================================

% Define handles
handles.zbc = 0;
handles.quad = 0;

% initialize
[nz,ny,nx] = size(SA);
%z2 = round(0.5*nz); %z2 = 20;
k_west = nan(size(SA));
r_west = k_west;

%			zonal bands ...
for j0 = 1:ny
    
    Idata = find(isfinite(SA(1,j0,:)));
    nn = length(Idata);
    if ~isempty(Idata)
 % meridional sweeps
        for i = 1:nn
            
            i0 = Idata(i);
            
            if i0 == 1
                if handles.zbc==1
                    i0_west = nx;
                else
                    i0_west = nan;
                end
            else
                i0_west = i0-1;
            end
            
            if isfinite(i0_west) & isfinite(SA(1,j0,i0_west))
                
% the central cast               
                n0 = n(j0,i0);
                SA0 = SA(1:n0,j0,i0);
                CT0 = CT(1:n0,j0,i0);
                p0 = p(1:n0,j0,i0);
                                
% the western cast                
                n0_west = n(j0,i0_west);
                SA_west = SA(1:n0_west,j0,i0_west);
                CT_west = CT(1:n0_west,j0,i0_west);
                p_west = p(1:n0_west,j0,i0_west);
                
% depth_ntp calculation                
                SAdns_west = nan(n0,1);
                CTdns_west = SAdns_west;
                pdns_west = SAdns_west;
                                
                for k0 = 1:n0
                    [SAdns_west(k0),CTdns_west(k0),pdns_west(k0)] = ...
                        gsw_depth_ntp(SA0(k0),CT0(k0),p0(k0),SA_west,CT_west,p_west);
                end

                tmp = pdns_west;
                pdns_west = nan(nz,1);
                pdns_west(1:n0) = tmp;
                
                %          intersection vitals
                for k0 = 1:n0
                        try
                            k_west(k0,j0,i0) = gamma_indx(p(1:n0_west,j0,i0_west),pdns_west(k0));
                            if pdns_west(k0) == p(n0_west,j0,i0_west)
                                k_west(k0,j0,i0) = n0-1;
                            end
                            kk = k_west(k0,j0,i0);
                            r_west(k0,j0,i0) = (pdns_west(k0)-p(kk,j0,i0_west))/ ...
                                (p(kk+1,j0,i0_west)-p(kk,j0,i0_west));
                            
                            if handles.quad == 1
                                r_ntp = r_west(k0,j0,i0); %[k0,j0,i0,r_ntp]
                                pmid = 0.5*(p(kk+1,j0,i0_west) + p(kk,j0,i0_west));                                
                                [rho_u,alpha_u,beta_u] = gsw_rho_alpha_beta(SA(kk,j0,i0_west),CT(kk,j0,i0_west),pmid);
                                [rho_l,alpha_l,beta_l] = gsw_rho_alpha_beta(SA(kk+1,j0,i0_west),CT(kk+1,j0,i0_west),pmid);                               
                                term1 = (alpha_u + 0.5*r_ntp*(alpha_l - alpha_u))*(CT(kk+1,j0,i0_west) - CT(kk,j0,i0_west));
                                term2 = (beta_u + 0.5*r_ntp*(beta_l - beta_u))*(SA(kk+1,j0,i0_west) - SA(kk,j0,i0_west));
                                r_west(k0,j0,i0) = 0.5*r_ntp*(term2 - term1)*(rho_l + rho_u)/(rho_l - rho_u);
                                %r_quad = r_west(k0,j0,i0)
                            end
                        end
                    %       switch off between oceans
                    oce_test = ocean(j0,i0)*ocean(j0,i0_west);
                    if oce_test == 3 | oce_test == 5
                        k_west(k0,j0,i0) = -1;
                        r_west(k0,j0,i0) = nan;
                    end
                end
            end
        end
    end
    
%     %ints_west plot
%     inds_write = find(isfinite(k_west));
%     if ~isempty(inds_write)
%         counter = counter + 1;
%         if mod(counter,nwrite)==0
%             %         na = length(find(isfinite(r_west(:,1:j0,:))));
%             %         ns = length(find(isfinite(s(:,1:j0,:))));
%             %         ok = [j0,100*j0/ny,100*na/ns]
%             
%             figure(4)
%             subplot(2,2,3)
%             zz = squeeze(k_west(nz2,:,:));
%             dj_pltmp(long(1,:),lat(:,1),zz)
%             
%             I_3 = find(zz==0);
%             if ~isempty(I_3)
%                 [j_triples,i_triples]= ind2sub([ny,nx],I_3);
%                 hold on
%                 plot(longss(i_triples),latss(j_triples),'m*')
%                 hold off
%             end
%             title('k\_west')
%             set(gcf,'Name','east/west intersections','NumberTitle','off','Color',[0.961 0.988 0.965])
%             
%             subplot(2,2,4)
%             dj_pltmp(long(1,:),lat(:,1),squeeze(r_west(nz2,:,:)))
%             title('r\_west')
% %             figure(gcf);
% %             dj_pause(1); %dj_toc
%         end
%     end
    
end

% save final intersections_west file

% figure(4)
% subplot(2,2,3)
% dj_pltmp(long(1,:),lat(:,1),squeeze(k_west(nz2,:,:)))
% title('k\_west')
% 
% subplot(2,2,4)
% dj_pltmp(long(1,:),lat(:,1),squeeze(r_west(nz2,:,:)))
% title('r\_west')
% 
% % triples
% I_3 = find(k_west==0);
% if ~isempty(I_3)
%     [k_triples,j_triples,i_triples]= ind2sub([nz,ny,nx],I_3);
%     hold on
%     plot(longss(i_triples),latss(j_triples),'m*')
%     hold off
%     figure(gcf)
% end

end
