function [k_west,r_west] = intersections_west(s,t,p,ocean,n,longs,lats,handles)

%%%   intersections_west  -  westerly isopycnal intersection statistics
%%%
%%%   Usage:      [k_west,r_west] = intersections_west(s,t,p,g,longs,lats)
%%%
%%%   Input:      s      - salinity	(PSU)			        {nz,ny,nx}
%%%   		      t      - conservative temperature (ITS90) {nz,ny,nx}
%%%               p      - pressure (db)                    {nz,ny,nx}
%%%               longs  - longitude (degE)                 {ny,nx}
%%%               lats   - latitude (degN)                  {ny,nx}
%%%
%%%   Output:     k_west - index of westerly intersection      {nz,ny,nx}
%%%               r_west - interpolation distance from index   {nz,ny,nx}
%%%
%%%   Method:     Sweep through data making westerly isopycnal surface calculations
%%%
%%%   Author:     David Jackett
%%%
%%%   Date:       06/04/06
%%%


counter = 0; nwrite = 10;

longss = nanmean(longs(:,:)); latss = nanmean(lats(:,:)');

check_cals = 0;

%%			initialize

[nz,ny,nx] = size(s); nz2 = round(nz/2); %z2 = 20;

k_west = nan*ones(size(s)); r_west = k_west;


%%			zonal bands ...

for j0 = 1:ny
    
    inds = find(isfinite(s(1,j0,:))); nn = length(inds);
    
    if nn>0
        
        %%			and meridional sweeps
        
        for i = 1:nn
            
            i0 = inds(i); long0 = longs(j0,i0);  lat0 = lats(j0,i0);
            
            if i0==1
                if handles.zbc==1
                    i0_west = nx;
                else
                    i0_west = nan;
                end
            else
                i0_west = i0-1;
            end
            
            if isfinite(i0_west) & isfinite(s(1,j0,i0_west))
                
                %%			the central cast
                
                n0 = n(j0,i0);
                
                s0 = s(1:n0,j0,i0);
                t0 = t(1:n0,j0,i0);
                p0 = p(1:n0,j0,i0);
                
                if check_cals == 1
                    location = [longss(i0),latss(j0)]
                    CENTRAL_CAST = [s0,t0,p0];
                end
                
                %%          the western cast
                
                n0_west = n(j0,i0_west);
                
                s_west = s(1:n0_west,j0,i0_west);
                t_west = t(1:n0_west,j0,i0_west);
                p_west = p(1:n0_west,j0,i0_west);
                
                if check_cals == 1
                    WESTERN_CAST = [s_west,t_west,p_west];
                end
                
                %%          depth_ntp calculation
                
                sdns_west = nan*ones(n0,1);
                tdns_west = sdns_west;
                pdns_west = sdns_west;
                
                if strcmp(handles.eos,'eos80_t')
                    eosno = 1;
                    for k0 = 1:n0
                        [sdns_west(k0),tdns_west(k0),pdns_west(k0)] = ...
                            depth_ntp(s0(k0),t0(k0),p0(k0),s_west,t_west,p_west);
                        %                                quick_depth_ntp(s0(k0),t0(k0),p0(k0),s_west,t_west,p_west,eosno);
                    end
                elseif strcmp(handles.eos,'eos05_th')
                    error('not yet implemented in intersections_west')
                elseif strcmp(handles.eos,'eos05_ct')
                    for k0 = 1:n0
                        [sdns_west(k0),tdns_west(k0),pdns_west(k0)] = ...
                            depth_ntp(s0(k0),t0(k0),p0(k0),s_west,t_west,p_west);
                    end
                end
                
                tmp = pdns_west; pdns_west = nan*ones(nz,1);
                pdns_west(1:n0) = tmp;
                
                if check_cals == 1
                    disp('WESTERN CAST DEPTH_NS COMPUTATIONS')
                    cast_vals = [pdns_west(:)]
                end
                
                
                %%          intersection vitals
                
                for k0 = 1:n0
                    if abs(pdns_west(k0)+99)<=1e-2
                        k_west(k0,j0,i0) = nan;
                        r_west(k0,j0,i0) = nan;
                    elseif abs(pdns_west(k0)+99.1)<=1e-2
                        k_west(k0,j0,i0) = 1;
                        r_west(k0,j0,i0) = nan;
                    elseif abs(pdns_west(k0)+99.2)<=1e-2
                        k_west(k0,j0,i0) = n0;
                        r_west(k0,j0,i0) = nan;
                    elseif abs(pdns_west(k0)+99.3)<=1e-2
                        %                    not_ok = 'a triple'
                        k_west(k0,j0,i0) = 0;
                        r_west(k0,j0,i0) = nan;
                    else
                        k_west(k0,j0,i0) = indx(p(1:n0_west,j0,i0_west),pdns_west(k0));
                        if pdns_west(k0)==p(n0_west,j0,i0_west)
                            k_west(k0,j0,i0) = n0-1;
                        end
                        kk = k_west(k0,j0,i0);
                        r_west(k0,j0,i0) = (pdns_west(k0)-p(kk,j0,i0_west))/...
                            (p(kk+1,j0,i0_west)-p(kk,j0,i0_west));
                        
                        if handles.quad==1
                            
                            r_ntp = r_west(k0,j0,i0); %[k0,j0,i0,r_ntp]
                            
                            pmid = (p(kk+1,j0,i0_west)+p(kk,j0,i0_west))/2;
                            
                            if strcmp(handles.eos,'eos80_t')
                                ct_u = ct_from_t(s(kk,j0,i0_west),t(kk,j0,i0_west),p(kk,j0,i0_west));
                                %                        [rho_u,rho_u_s,rho_u_t,rho_u_p] = eosall_from_ct(s(kk,j0,i0_west),ct,p(kk,j0,i0_west));
                                [rho_u,rho_u_s,rho_u_t,rho_u_p] = eosall_from_ct(s(kk,j0,i0_west),ct_u,pmid);
                                %                        rho_u_local = rho_from_ct(s(kk,j0,i0_west),ct,pmid);
                            else
                                'not yet implemented in intersections_west'
                            end
                            alfa_u = -rho_u_t/rho_u; beta_u = rho_u_s/rho_u;
                            if strcmp(handles.eos,'eos80_t')
                                ct_l = ct_from_t(s(kk+1,j0,i0_west),t(kk+1,j0,i0_west),p(kk+1,j0,i0_west));
                                %                        [rho_l,rho_l_s,rho_l_t,rho_l_p] = eosall_from_ct(s(kk+1,j0,i0_west),ct,p(kk+1,j0,i0_west));
                                [rho_l,rho_l_s,rho_l_t,rho_l_p] = eosall_from_ct(s(kk+1,j0,i0_west),ct_l,pmid);
                                %                        rho_l_local = rho_from_ct(s(kk+1,j0,i0_west),ct,pmid);
                            else
                                'not yet implemented in intersections_west'
                            end
                            
                            alfa_l = -rho_l_t/rho_l; beta_l = rho_l_s/rho_l;
                            term1 = (alfa_u+r_ntp*(alfa_l-alfa_u)/2)*(ct_l-ct_u)/(rho_l-rho_u);
                            term2 = (beta_u+r_ntp*(beta_l-beta_u)/2)*(s(kk+1,j0,i0_west)-s(kk,j0,i0_west))/(rho_l-rho_u);
                            r_west(k0,j0,i0) = r_ntp*(rho_l+rho_u)*(-term1+term2)/2; %r_quad = r_west(k0,j0,i0), dj_pause(0)
                            
                        end
                        
                    end
                    
                    %%       switch off between oceans
                    
                    oce_test = ocean(j0,i0)*ocean(j0,i0_west);
                    
                    if oce_test==3 | oce_test==5
                        k_west(k0,j0,i0) = -1;
                        r_west(k0,j0,i0) = nan;
                    end
                end
            end
        end
    end
    
    
    %%			ints_west plot
    
    inds_write = find(isfinite(k_west));
    if length(inds_write)>0
        counter = counter+1;
        if mod(counter,nwrite)==0
            na = length(find(isfinite(r_west(:,1:j0,:))));
            ns = length(find(isfinite(s(:,1:j0,:))));
            ok = [j0,100*j0/ny,100*na/ns]
            zz = squeeze(k_west(nz2,:,:));
            figure(2),clf, [hax,hax_noxlabel,hax_noylabel,hax_left,hax_bottom] = nfigaxes([1 2],[0.08 0.08],[0.08 0.92],[0.08 0.92]);
            axes(hax(1)) %subplot(2,2,3)
            dj_pltmp(longs(1,:),lats(:,1),zz,0)
            inds_3 = find(zz==0);
            if length(inds_3)>0
                [j_triples,i_triples]= ind2sub([ny,nx],inds_3);
                hold on
                plot(longss(i_triples),latss(j_triples),'m*')
                hold off
            end
            title('k\_west')
            set(gcf,'Name','West intersections','NumberTitle','off','Color',[0.961 0.988 0.965])
            
            zz = squeeze(r_west(nz2,:,:));
            figure(2), axes(hax(2))%subplot(2,2,4)
            dj_pltmp(longs(1,:),lats(:,1),zz,0)
            title('r\_west')
            figure(2); dj_pause(1); %dj_toc
        end
    end
    
end


%%			and save final intersections_west file

zz = squeeze(k_west(nz2,:,:));
figure(2)
%clf
% [hax,hax_noxlabel,hax_noylabel,hax_left,hax_bottom] = nfigaxes([2 2],[0.08 0.08],[0.08 0.92],[0.08 0.92]);
axes(hax(1))
hold on
%figure(1), subplot(2,2,3)
dj_pltmp(longs(1,:),lats(:,1),zz,0)
title('k\_west')

zz = squeeze(r_west(nz2,:,:));
figure(2), axes(hax(2))%subplot(2,2,4)
hold on
dj_pltmp(longs(1,:),lats(:,1),zz,0)
title('r\_west')


%%          triples

inds_3 = find(k_west==0);
if length(inds_3)>0
    [k_triples,j_triples,i_triples]= ind2sub([nz,ny,nx],inds_3);
    hold on
    plot(longss(i_triples),latss(j_triples),'m*')
    hold off
    figure(2)
end

figure(2); dj_toc; dj_pause(1)

return
