function [k_east,r_east] = intersections_east(SA,CT,p,ocean,n,long,lat,handles)

% intersections_east  -  easterly isopycnal intersection statistics
%
% Usage:
%  [k_east,r_east] = intersections_east(SA,CT,p,g,ocean,n,longs,lats)
%
% Description:
%  Sweep through data making easterly neutral tangent plane calculations
%
% Input:
%  SA      - Absolute Salinity	         [ g/kg ]            {nz,ny,nx}
%  CT      - Conservative Temperature    [ deg C ]           {nz,ny,nx}
%  p       - pressure                    [ dbar ]            {nz,ny,nx}
%  longs   - longitude                   [ deg E ]           {ny,nx}
%  lats    - latitude                    [ deg N ]           {ny,nx}
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
handles.zbc = 0;
handles.quad = 0;

counter = 0;
nwrite = 1;

longss = nanmean(long(:,:));
latss = nanmean(lat(:,:)');

check_calcs = 0;
% whos

% initialize
[nz,ny,nx] = size(SA);
nz2 = round(0.5*nz); %nz2 = 40;
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
            long0 = long(j0,i0);
            lat0 = lat(j0,i0);
            
            if i0 < nx
                i0_east = i0+1;
            else
                if handles.zbc == 1
                    i0_east = 1;
                else
                    i0_east = nan;
                end
            end
            
            if isfinite(i0_east) & isfinite(SA(1,j0,i0_east))
                % the central cast
                n0 = n(j0,i0);
                SA0 = SA(1:n0,j0,i0);
                CT0 = CT(1:n0,j0,i0);
                p0 = p(1:n0,j0,i0);
                
                if check_calcs == 1
                    location_central = [long(j0,i0),lat(j0,i0)]
                    CENTRAL_CAST = [SA0,CT0,p0]
                end
                
                % the eastern cast
                n0_east = n(j0,i0_east);
                SA_east = SA(1:n0_east,j0,i0_east);
                CT_east = CT(1:n0_east,j0,i0_east);
                p_east = p(1:n0_east,j0,i0_east);
                
                if check_calcs == 1
                    EASTERN_CAST = [SA_east,CT_east,p_east]
                end
                
                % depth_ns calculation
                SAdns_east = nan(n0,1);
                CTdns_east = SAdns_east;
                pdns_east = SAdns_east;
                
                %                 if strcmp(handles.eos,'eos80_t')
                %                     eosno = 1;
                %                     for k0 = 1:n0
                %                         [sdns_east(k0),tdns_east(k0),pdns_east(k0)] = ...
                %                             depth_ntp(s0(k0),t0(k0),p0(k0),s_east,t_east,p_east);
                %                         %                          	       quick_depth_ntp(s0(k0),t0(k0),p0(k0),s_east,t_east,p_east,eosno);
                %                     end
                %                 elseif strcmp(handles.eos,'eos05_th')
                %                     error('not yet implemented in intersections_east')
                %                 elseif strcmp(handles.eos,'eos05_ct')
                %                     for k0 = 1:n0
                %                         [sdns_east(k0),tdns_east(k0),pdns_east(k0)] = ...
                %                             depth_ntp(s0(k0),t0(k0),p0(k0),s_east,t_east,p_east);
                %                     end
                %                 end
                
                for k0 = 1:n0
                    [SAdns_east(k0),CTdns_east(k0),pdns_east(k0)] = ...
                        gsw_depth_ntp(SA0(k0),CT0(k0),p0(k0),SA_east,CT_east,p_east);
                end
                
                tmp = pdns_east;
                pdns_east = nan(nz,1);
                pdns_east(1:n0) = tmp;
                
                if check_calcs == 1
                    disp('EASTERN CAST DEPTH_NS COMPUTATIONS')
                    cast_vals = [pdns_east(:)]
                end
                
                % intersection vitals
                for k0 = 1:n0
                    if abs(pdns_east(k0)+99) <= 1e-2
                        k_east(k0,j0,i0) = nan;
                        r_east(k0,j0,i0) = nan;
                    elseif abs(pdns_east(k0)+99.1) <= 1e-2
                        k_east(k0,j0,i0) = 1;
                        r_east(k0,j0,i0) = nan;
                    elseif abs(pdns_east(k0)+99.2) <= 1e-2
                        k_east(k0,j0,i0) = n0;
                        r_east(k0,j0,i0) = nan;
                    elseif abs(pdns_east(k0)+99.3) <= 1e-2
                        % not_ok = 'a triple'
                        k_east(k0,j0,i0) = 0;
                        r_east(k0,j0,i0) = nan;
                    else
                        try
                            k_east(k0,j0,i0) = indx(p(1:n0_east,j0,i0_east),pdns_east(k0));
                            if pdns_east(k0) == p(n0_east,j0,i0_east)
                                k_east(k0,j0,i0) = n0_east - 1;
                            end
                            kk = k_east(k0,j0,i0);
                            r_east(k0,j0,i0) = (pdns_east(k0) - p(kk,j0,i0_east))/...
                                (p(kk+1,j0,i0_east) - p(kk,j0,i0_east));
                            
                            if handles.quad == 1
                                r_ntp = r_east(k0,j0,i0); %[k0,j0,i0,r_ntp]
                                pmid = 0.5*(p(kk+1,j0,i0_east) + p(kk,j0,i0_east));
                                
                                %                             if strcmp(handles.eos,'eos80_t')
                                %                                 ct_u = ct_from_t(s(kk,j0,i0_east),t(kk,j0,i0_east),p(kk,j0,i0_east));
                                %                                 [rho_u,rho_u_s,rho_u_t,rho_u_p] = eosall_from_ct(s(kk,j0,i0_east),ct_u,p(kk,j0,i0_east));
                                %                                 rho_local_u = rho_from_ct(s(kk,j0,i0_east),ct_u,pmid);
                                %                             else
                                %                                 'not yet implemented in intersections_east'
                                %                             end
                                %                             alfa_u = -rho_u_t/rho_u;
                                %                             beta_u = rho_u_s/rho_u;
                                [rho_u,alpha_u,beta_u] = gsw_rho_alpha_beta(SA(kk,j0,i0_east),CT(kk,j0,i0_east),pmid);
                                %rho_local_u = gsw_rho(SA(kk,j0,i0_east),CT(kk,j0,i0_east),pmid);
                                
                                %                             if strcmp(handles.eos,'eos80_t')
                                %                                 ct_l = ct_from_t(s(kk+1,j0,i0_east),t(kk+1,j0,i0_east),p(kk+1,j0,i0_east));
                                %                                 [rho_l,rho_l_s,rho_l_t,rho_l_p] = eosall_from_ct(s(kk+1,j0,i0_east),ct_l,p(kk+1,j0,i0_east));
                                %                                 rho_local_l = rho_from_ct(s(kk+1,j0,i0_east),ct_l,pmid);
                                %                             else
                                %                                 'not yet implemented in intersections_east'
                                %                             end
                                %                             alfa_l = -rho_l_t/rho_l;
                                %                             beta_l = rho_l_s/rho_l;
                                [rho_l,alpha_l,beta_l] = gsw_rho_alpha_beta(SA(kk+1,j0,i0_east),CT(kk+1,j0,i0_east),pmid);
                                %rho_local_l = gsw_rho(SA(kk+1,j0,i0_east),CT(kk+1,j0,i0_east),pmid);
                                
                                term1 = (alpha_u + 0.5*r_ntp*(alpha_l - alpha_u))*(CT(kk+1,j0,i0_east) - CT(kk,j0,i0_east));
                                term2 = (beta_u + 0.5*r_ntp*(beta_l - beta_u))*(SA(kk+1,j0,i0_east) - SA(kk,j0,i0_east));
                                r_east(k0,j0,i0) = 0.5*r_ntp*(term2 - term1)*(rho_l + rho_u)/(rho_l - rho_u);
                                %r_quad = r_east(k0,j0,i0)
                                %dj_pause(0)
                                
                            end
                        end
                    end
                    
                    % switch off between oceans
                    oce_test = ocean(j0,i0)*ocean(j0,i0_east);
                    
                    if (oce_test == 3) | (oce_test == 5)
                        k_east(k0,j0,i0) = -1;
                        r_east(k0,j0,i0) = nan;
                    end
                end
            end
        end
    end
    
    % intersections_east plot
    Iplot = find(isfinite(k_east));
    if ~isempty(Iplot)
        counter = counter + 1;
        if mod(counter,nwrite) == 0
            % disp(' ');
            % done = [j0,i0,i]
            % na = length(find(isfinite(r_east(:,1:j0,:))));
            % ns = length(find(isfinite(SA(:,1:j0,:))));
            % ok = [j0,lats(j0),100*j0/ny,100*na/ns]
            
            figure(1)
            subplot(2,2,1)
            zz = squeeze(k_east(nz2,:,:));
            dj_pltmp(long(1,:),lat(:,1),zz)
            I_3 = find(zz == 0);
            %length(I_3)
            if ~isempty(I_3)
                [j_triples,i_triples] = ind2sub([ny,nx],I_3);
                hold on
                plot(longss(i_triples),latss(j_triples),'m*')
                hold off
            end
            title('k\_east')
            set(gcf,'Name','east/west intersections','NumberTitle','off','Color',[0.961 0.988 0.965])
            
            % figure(1)
            subplot(2,2,2)
            dj_pltmp(longss,latss,squeeze(r_east(nz2,:,:)))
            title('r\_east')
            figure(gcf);
            %dj_pause(0);
            %dj_toc
        end
    end
end

% save final intersections_east file

figure(1)
subplot(2,2,1)
dj_pltmp(long(1,:),lat(:,1),squeeze(k_east(nz2,:,:)))
title('k\_east')

subplot(2,2,2)
dj_pltmp(long(1,:),lat(:,1),squeeze(r_east(nz2,:,:)))
title('r\_east')

I_3 = find(k_east==0);
if ~isempty(I_3)
    [k_triples,j_triples,i_triples] = ind2sub([nz,ny,nx],I_3);
    hold on
    plot(longss(i_triples),latss(j_triples),'m*')
    hold off
end

% figure(gcf);
% dj_toc;
% dj_pause(0)

end
%return
