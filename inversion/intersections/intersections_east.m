function [k_east,r_east] = intersections_east(s,t,p,ocean,n,longs,lats,handles)

%%%   intersections_east  -  easterly isopycnal intersection statistics
%%%
%%%   Usage:         [k_east,r_east] = intersections_east(s,t,p,g,ocean,n,longs,lats)
%%%
%%%   Input:         s       - salinity	(PSU)                           {nz,ny,nx}
%%%                  t       - in situ temperature (ITS90)              {nz,ny,nx}
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


counter = 0; nwrite = 1;

longss = nanmean(longs(:,:)); latss = nanmean(lats(:,:)');

check_cals = 0; whos


%%			initialize

[nz,ny,nx] = size(s); nz2 = round(nz/2); %nz2 = 40;

k_east = nan*ones(size(s)); r_east = k_east;


%%			zonal bands ...

for j0 = 1:ny
    
    inds = find(isfinite(s(1,j0,:))); nn = length(inds);
    if nn>0
        
        %%			and meridional sweeps ...
        
        for i = 1:nn
            
            i0 = inds(i); long0 = longs(j0,i0);  lat0 = lats(j0,i0);
            
            if i0<nx
                i0_east = i0+1;
            else
                if handles.zbc==1
                    i0_east = 1;
                else
                    i0_east = nan;
                end
            end
            
            if isfinite(i0_east) & isfinite(s(1,j0,i0_east))
                
                %%			the central cast
                
                n0 = n(j0,i0);
                
                s0 = s(1:n0,j0,i0);
                t0 = t(1:n0,j0,i0);
                p0 = p(1:n0,j0,i0);
                
                if check_cals == 1
                    location_central = [longs(j0,i0),lats(j0,i0)]
                    CENTRAL_CAST = [s0,t0,p0]
                end
                
                %%          the eastern cast
                
                n0_east = n(j0,i0_east);
                
                s_east = s(1:n0_east,j0,i0_east);
                t_east = t(1:n0_east,j0,i0_east);
                p_east = p(1:n0_east,j0,i0_east);
                
                if check_cals == 1
                    EASTERN_CAST = [s_east,t_east,p_east]
                end
                
                %%          depth_ns calculation
                
                sdns_east = nan*ones(n0,1);
                tdns_east = sdns_east;
                pdns_east = sdns_east;
                
                if strcmp(handles.eos,'eos80_t')
                    eosno = 1;
                    for k0 = 1:n0
                        [sdns_east(k0),tdns_east(k0),pdns_east(k0)] = ...
                            depth_ntp(s0(k0),t0(k0),p0(k0),s_east,t_east,p_east);
                        %                          	       quick_depth_ntp(s0(k0),t0(k0),p0(k0),s_east,t_east,p_east,eosno);
                        
                    end
                elseif strcmp(handles.eos,'eos05_th')
                    error('not yet implemented in intersections_east')
                elseif strcmp(handles.eos,'eos05_ct')
                    for k0 = 1:n0
                        [sdns_east(k0),tdns_east(k0),pdns_east(k0)] = ...
                            depth_ntp(s0(k0),t0(k0),p0(k0),s_east,t_east,p_east);
                    end
                end
                
                tmp = pdns_east; pdns_east = nan*ones(nz,1);
                pdns_east(1:n0) = tmp;
                
                if check_cals == 1
                    disp('EASTERN CAST DEPTH_NS COMPUTATIONS')
                    cast_vals = [pdns_east(:)]
                end
                
                
                %%          intersection vitals
                for k0 = 1:n0
                    if abs(pdns_east(k0)+99)<=1e-2
                        k_east(k0,j0,i0) = nan;
                        r_east(k0,j0,i0) = nan;
                    elseif abs(pdns_east(k0)+99.1)<=1e-2
                        k_east(k0,j0,i0) = 1;
                        r_east(k0,j0,i0) = nan;
                    elseif abs(pdns_east(k0)+99.2)<=1e-2
                        k_east(k0,j0,i0) = n0;
                        r_east(k0,j0,i0) = nan;
                    elseif abs(pdns_east(k0)+99.3)<=1e-2
                        %                    not_ok = 'a triple'
                        k_east(k0,j0,i0) = 0;
                        r_east(k0,j0,i0) = nan;
                    else
                        k_east(k0,j0,i0) = indx(p(1:n0_east,j0,i0_east),pdns_east(k0));
                        if pdns_east(k0)==p(n0_east,j0,i0_east)
                            k_east(k0,j0,i0) = n0_east-1;
                        end
                        kk = k_east(k0,j0,i0);
                        r_east(k0,j0,i0) = (pdns_east(k0)-p(kk,j0,i0_east))/...
                            (p(kk+1,j0,i0_east)-p(kk,j0,i0_east));
                        
                        if handles.quad==1
                            
                            r_ntp = r_east(k0,j0,i0); %[k0,j0,i0,r_ntp]
                            
                            pmid = (p(kk+1,j0,i0_east)+p(kk,j0,i0_east))/2;
                            
                            if strcmp(handles.eos,'eos80_t')
                                ct_u = ct_from_t(s(kk,j0,i0_east),t(kk,j0,i0_east),p(kk,j0,i0_east));
                                [rho_u,rho_u_s,rho_u_t,rho_u_p] = eosall_from_ct(s(kk,j0,i0_east),ct_u,p(kk,j0,i0_east));
                                rho_local_u = rho_from_ct(s(kk,j0,i0_east),ct_u,pmid);
                            else
                                'not yet implemented in intersections_east'
                            end
                            alfa_u = -rho_u_t/rho_u; beta_u = rho_u_s/rho_u;
                            
                            if strcmp(handles.eos,'eos80_t')
                                ct_l = ct_from_t(s(kk+1,j0,i0_east),t(kk+1,j0,i0_east),p(kk+1,j0,i0_east));
                                [rho_l,rho_l_s,rho_l_t,rho_l_p] = eosall_from_ct(s(kk+1,j0,i0_east),ct_l,p(kk+1,j0,i0_east));
                                rho_local_l = rho_from_ct(s(kk+1,j0,i0_east),ct_l,pmid);
                            else
                                'not yet implemented in intersections_east'
                            end
                            
                            alfa_l = -rho_l_t/rho_l; beta_l = rho_l_s/rho_l;
                            term1 = (alfa_u+r_ntp*(alfa_l-alfa_u)/2)*(ct_l-ct_u)/(rho_l-rho_u);
                            term2 = (beta_u+r_ntp*(beta_l-beta_u)/2)*(s(kk+1,j0,i0_east)-s(kk,j0,i0_east))/(rho_l-rho_u);
                            r_east(k0,j0,i0) = r_ntp*(rho_l+rho_u)*(-term1+term2)/2; %r_quad = r_east(k0,j0,i0), dj_pause(0)
                            
                        end
                        
                    end
                    
                   %%       switch off between oceans
                    
                    oce_test = ocean(j0,i0)*ocean(j0,i0_east);
                    
                    if oce_test==3 | oce_test==5
                        k_east(k0,j0,i0) = -1;
                        r_east(k0,j0,i0) = nan;
                    end
                    
                end
            end
        end
    end
    
    
    %%			ints_east plot
    
    inds_write = find(isfinite(k_east));
    if length(inds_write)>0
        counter = counter+1;
        if mod(counter,nwrite)==0
            %	    disp(' '); done = [j0,i0,i]
            na = length(find(isfinite(r_east(:,1:j0,:))));
            ns = length(find(isfinite(s(:,1:j0,:))));
            %        ok = [j0,lats(j0),100*j0/ny,100*na/ns]
            zz = squeeze(k_east(nz2,:,:));
            figure(1)
            clf
            [hax,hax_noxlabel,hax_noylabel,hax_left,hax_bottom] = nfigaxes([1 2],[0.08 0.08],[0.08 0.92],[0.08 0.92]);
            axes(hax(1))
            
            %figure(1), subplot(2,2,1)
            dj_pltmp(longs(1,:),lats(:,1),zz,0)
            inds_3 = find(zz==0); %length(inds_3)
            if length(inds_3)>0
                [j_triples,i_triples]= ind2sub([ny,nx],inds_3);
                hold on
                plot(longss(i_triples),latss(j_triples),'m*')
                hold off
            end
            title('k\_east')
            set(gcf,'Name','east intersections','NumberTitle','off','Color',[0.961 0.988 0.965])
            
            zz = squeeze(r_east(nz2,:,:));
            figure(1), axes(hax(2))%subplot(2,2,2)
            dj_pltmp(longss,latss,zz,0)
            title('r\_east')
            figure(1); %dj_pause(0); %dj_toc
        end
    end
    
end


%%			and save final intersections_east file

zz = squeeze(k_east(nz2,:,:));
figure(1)
axes(hax(1))
hold on
%figure(1), subplot(2,2,1)
dj_pltmp(longs(1,:),lats(:,1),zz,0)
title('k\_east')

zz = squeeze(r_east(nz2,:,:));
figure(1), axes(hax(2)) %subplot(2,2,2)
hold on
dj_pltmp(longs(1,:),lats(:,1),zz,0)
title('r\_east')

inds_3 = find(k_east==0);
if length(inds_3)>0
    [k_triples,j_triples,i_triples]= ind2sub([nz,ny,nx],inds_3);
    hold on
    plot(longss(i_triples),latss(j_triples),'m*')
    hold off
end

figure(1); dj_toc; dj_pause(1)

return


