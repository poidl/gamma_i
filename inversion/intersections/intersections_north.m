function [k_north,r_north] = intersections_north(s,t,p,ocean,n,longs,lats,handles)

%%%   intersections_north  -  northerly isopycnal intersection statistics
%%%
%%%   Usage:      [k_north,r_north] = intersections_north(s,t,p,g,longs,lats)
%%%
%%%   Input:      s      - salinity	(PSU)			        {nz,ny,nx}
%%%   		      t      - conservative temperature (ITS90) {nz,ny,nx}
%%%               p      - pressure (db)                    {nz,ny,nx}
%%%               longs  - longitude (degE)                 {ny,nx}
%%%               lats   - latitude (degN)                  {ny,nx}
%%%
%%%   Output:     k_north - index of northerly intersection     {nz,ny,nx}
%%%               r_north - interpolation distance from index   {nz,ny,nx}
%%%
%%%   Method:     Sweep through data making northerly isopycnal surface calculations
%%%
%%%   Author:     David Jackett
%%%
%%%   Date:       6/4/06
%%%


counter = 0; nwrite = 10;

longss = nanmean(longs(:,:)); latss = nanmean(lats(:,:)');

check_cals = 0;

%%			initialize

[nz,ny,nx] = size(s); nz2 = round(nz/2); %nz2 = 20;

k_north = nan*ones(size(s)); r_north = k_north;


%%			zonal bands ...

for j0 = 1:ny-1
    
    j0_north = j0+1;
    
    inds = find(isfinite(s(1,j0,:))); nn = length(inds);
    
    if nn>0
        
        %%			and meridional sweeps
        
        for i = 1:nn
            
            i0 = inds(i); long0 = longs(j0,i0);  lat0 = lats(j0,i0);
            
            if isfinite(s(1,j0_north,i0))
                
                %%			the central cast
                
                n0 = n(j0,i0);
                
                s0 = s(1:n0,j0,i0);
                t0 = t(1:n0,j0,i0);
                p0 = p(1:n0,j0,i0);
                
                if check_cals == 1
                    location = [longss(i0),latss(j0)]
                    CENTRAL_CAST = [s0,t0,p0]
                end
                
                %%          the northern cast
                
                n0_north = n(j0_north,i0);
                
                s_north = s(1:n0_north,j0_north,i0);
                t_north = t(1:n0_north,j0_north,i0);
                p_north = p(1:n0_north,j0_north,i0);
                
                if check_cals == 1
                    NORTHERN_CAST = [s_north,t_north,p_north]
                end
                
                %%          depth_ntp calculation
                
                sdns_north = nan*ones(n0,1);
                tdns_north = sdns_north;
                pdns_north = sdns_north;
                
                if strcmp(handles.eos,'eos80_t')
                    eosno = 1;
                    for k0 = 1:n0
                        [sdns_north(k0),tdns_north(k0),pdns_north(k0)] = ...
                            depth_ntp(s0(k0),t0(k0),p0(k0),s_north,t_north,p_north);
                        %                                quick_depth_ntp(s0(k0),t0(k0),p0(k0),s_north,t_north,p_north,eosno);
                    end
                elseif strcmp(handles.eos,'eos05_th')
                    error('not yet implemented in intersections_north')
                elseif strcmp(handles.eos,'eos05_ct')
                    for k0 = 1:n0
                        [sdns_north(k0),tdns_north(k0),pdns_north(k0)] = ...
                            depth_ntp(s0(k0),t0(k0),p0(k0),s_north,t_north,p_north);
                    end
                end
                
                tmp = pdns_north; pdns_north = nan*ones(nz,1);
                pdns_north(1:n0) = tmp;
                
                if check_cals == 1
                    disp('NORTHERN CAST DEPTH_NS COMPUTATIONS')
                    cast_vals = [pdns_north(:)]
                end
                
                
                %%          intersection vitals
                
                for k0 = 1:n0
                    if abs(pdns_north(k0)+99)<=1e-2
                        k_north(k0,j0,i0) = nan;
                        r_north(k0,j0,i0) = nan;
                    elseif abs(pdns_north(k0)+99.1)<=1e-2
                        k_north(k0,j0,i0) = 1;
                        r_north(k0,j0,i0) = nan;
                    elseif abs(pdns_north(k0)+99.2)<=1e-2
                        k_north(k0,j0,i0) = n0;
                        r_north(k0,j0,i0) = nan;
                    elseif abs(pdns_north(k0)+99.3)<=1e-2
                        %                    not_ok = 'a triple'
                        k_north(k0,j0,i0) = 0;
                        r_north(k0,j0,i0) = nan;
                    else
                        k_north(k0,j0,i0) = indx(p(1:n0_north,j0_north,i0),pdns_north(k0));
                        if pdns_north(k0)==p(n0_north,j0_north,i0)
                            k_north(k0,j0,i0) = n0_north-1;
                        end
                        kk = k_north(k0,j0,i0);
                        r_north(k0,j0,i0) = (pdns_north(k0)-p(kk,j0_north,i0))/...
                            (p(kk+1,j0_north,i0)-p(kk,j0_north,i0));
                        
                        if handles.quad==1
                            
                            r_ntp = r_north(k0,j0,i0); %[k0,j0,i0,r_ntp]
                            
                            pmid = (p(kk+1,j0_north,i0)+p(kk,j0_north,i0))/2;
                            
                            if strcmp(handles.eos,'eos80_t')
                                ct_u = ct_from_t(s(kk,j0_north,i0),t(kk,j0_north,i0),p(kk,j0_north,i0));
                                %                        [rho_u,rho_u_s,rho_u_t,rho_u_p] = eosall_from_ct(s(kk,j0_north,i0),ct,p(kk,j0_north,i0));
                                [rho_u,rho_u_s,rho_u_t,rho_u_p] = eosall_from_ct(s(kk,j0_north,i0),ct_u,pmid);
                                %                        rho_u_local = rho_from_ct(s(kk,j0_north,i0),ct,pmid);
                            else
                                'not yet implemented in intersections_north'
                            end
                            alfa_u = -rho_u_t/rho_u; beta_u = rho_u_s/rho_u;
                            if strcmp(handles.eos,'eos80_t')
                                ct_l = ct_from_t(s(kk+1,j0_north,i0),t(kk+1,j0_north,i0),p(kk+1,j0_north,i0));
                                %                        [rho_l,rho_l_s,rho_l_t,rho_l_p] = eosall_from_ct(s(kk+1,j0_north,i0),ct,p(kk+1,j0_north,i0));
                                [rho_l,rho_l_s,rho_l_t,rho_l_p] = eosall_from_ct(s(kk+1,j0_north,i0),ct_l,pmid);
                                %                        rho_l_local = rho_from_ct(s(kk+1,j0_north,i0),ct,pmid);
                            else
                                'not yet implemented in intersections_north'
                            end
                            
                            alfa_l = -rho_l_t/rho_l; beta_l = rho_l_s/rho_l;
                            term1 = (alfa_u+r_ntp*(alfa_l-alfa_u)/2)*(ct_l-ct_u)/(rho_l-rho_u);
                            term2 = (beta_u+r_ntp*(beta_l-beta_u)/2)*(s(kk+1,j0_north,i0)-s(kk,j0_north,i0))/(rho_l-rho_u);
                            r_north(k0,j0,i0) = r_ntp*(rho_l+rho_u)*(-term1+term2)/2; %r_quad = r_north(k0,j0,i0), dj_pause(0)
                            
                        end
                        
                    end
                    
                    %%       switch off between oceans
                    
                    oce_test = ocean(j0,i0)*ocean(j0_north,i0);
                    
                    if oce_test==3 | oce_test==5
                        k_north(k0,j0,i0) = -1;
                        r_north(k0,j0,i0) = nan;
                    end
                    
                end
            end
        end
    end
    
    
    %%			ints_north plot
    
    inds_write = find(isfinite(k_north));
    if length(inds_write)>0
        counter = counter+1;
        if mod(counter,nwrite)==0
            %	    disp(' '); done = [j0,i0,i]
            na = length(find(isfinite(r_north(:,1:j0,:))));  ns = length(find(isfinite(s(:,1:j0,:))));
            ok = [j0,lats(j0),100*j0/ny,100*na/ns]
            zz = squeeze(k_north(nz2,:,:));
            figure(3),clf, [hax,hax_noxlabel,hax_noylabel,hax_left,hax_bottom] = nfigaxes([1 2],[0.08 0.08],[0.08 0.92],[0.08 0.92]);
            axes(hax(1)) %subplot(2,2,1)
            dj_pltmp(longs(1,:),lats(:,1),zz,0)
            inds_3 = find(zz==0);
            if length(inds_3)>0
                [j_triples,i_triples]= ind2sub([ny,nx],inds_3);
                hold on
                plot(longss(i_triples),latss(j_triples),'m*')
                hold off
            end
            title('k\_north')
            set(gcf,'Name','North intersections','NumberTitle','off','Color',[0.961 0.988 0.965])
            
            zz = squeeze(r_north(nz2,:,:));
            figure(3), axes(hax(2))%subplot(2,2,2)
            dj_pltmp(longs(1,:),lats(:,1),zz,0)
            title('r\_north')
            figure(3); dj_pause(1);
        end
    end
    
end


%%			and save final intersections_north file

zz = squeeze(k_north(nz2,:,:));
figure(3), axes(hax(1))%subplot(2,2,1)
dj_pltmp(longs(1,:),lats(:,1),zz,0)
title('k\_north')

zz = squeeze(r_north(nz2,:,:));
figure(3), axes(hax(2))%subplot(2,2,2)
dj_pltmp(longs(1,:),lats(:,1),zz,0)
title('r\_north')

%%          triples

inds_3 = find(k_north==0);
if length(inds_3)>0
    [k_triples,j_triples,i_triples]= ind2sub([nz,ny,nx],inds_3);
    hold on
    plot(longss(i_triples),latss(j_triples),'m*')
    hold off
end

figure(3); dj_toc; dj_pause(1)

return
