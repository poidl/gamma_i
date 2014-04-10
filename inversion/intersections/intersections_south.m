function [k_south,r_south] = intersections_south(s,t,p,ocean,n,longs,lats,handles)

%%%   intersections_south  -  southerly isopycnal intersection statistics
%%%
%%%   Usage:      [k_south,r_south] = intersections_south(s,t,p,g,longs,lats)
%%%
%%%   Input:      s       - salinity	(PSU)			                    {nz,ny,nx}
%%%               t       - conservative temperature (ITS90)                {nz,ny,nx}
%%%               p       - pressure (db)                                   {nz,ny,nx}
%%%               longs   - longitude (degE)                                {nx}
%%%               lats    - latitude (degN)                                 {ny}
%%%
%%%   Output:     k_south - index of southerly intersection                 {nz,ny,nx}
%%%               r_south - interpolation distance from index               {nz,ny,nx}
%%%
%%%   Method:     Sweep through data making southerly isopycnal surface calculations
%%%
%%%   Author:     David Jackett
%%%
%%%   Date:       06/04/06
%%%


counter = 0; nwrite = 10;

longss = nanmean(longs(:,:)); latss = nanmean(lats(:,:)');

check_cals = 0;

%%			initialize

[nz,ny,nx] = size(s); nz2 = round(nz/2); %nz2 = 20;

k_south = nan*ones(size(s)); r_south = k_south;


%%			zonal bands ...

for j0 = 2:ny
    
    j0_south = j0-1;
    
    inds = find(isfinite(s(1,j0,:))); nn = length(inds);
    
    if nn>0
        
        %%			and meridional sweeps
        
        for i = 1:nn
            
            i0 = inds(i); long0 = longs(j0,i0);  lat0 = lats(j0,i0);
            
            if isfinite(s(1,j0_south,i0))
                
                %%			the central cast
                
                n0 = n(j0,i0);
                
                s0 = s(1:n0,j0,i0);
                t0 = t(1:n0,j0,i0);
                p0 = p(1:n0,j0,i0);
                
                if check_cals == 1
                    location = [longss(i0),latss(j0)]
                    CENTRAL_CAST = [s0,t0,p0]
                end
                
                %%          the southern cast
                
                n0_south = n(j0_south,i0);
                
                s_south = s(1:n0_south,j0_south,i0);
                t_south = t(1:n0_south,j0_south,i0);
                p_south = p(1:n0_south,j0_south,i0);
                
                if check_cals == 1
                    SOUTHERN_CAST = [s_south,t_south,p_south]
                end
                
                %%          depth_ntp calculation
                
                sdns_south = nan*ones(n0,1);
                tdns_south = sdns_south;
                pdns_south = sdns_south;
                
                if strcmp(handles.eos,'eos80_t')
                    eosno = 1;
                    for k0 = 1:n0
                        [sdns_south(k0),tdns_south(k0),pdns_south(k0)] = ...
                            depth_ntp(s0(k0),t0(k0),p0(k0),s_south,t_south,p_south);
                        %                          	       quick_depth_ntp(s0(k0),t0(k0),p0(k0),s_south,t_south,p_south,eosno);
                    end
                elseif strcmp(handles.eos,'eos05_th')
                    error('not yet implemented in intersections_south')
                elseif strcmp(handles.eos,'eos05_ct')
                    for k0 = 1:n0
                        [sdns_south(k0),tdns_south(k0),pdns_south(k0)] = ...
                            depth_ntp(s0(k0),t0(k0),p0(k0),s_south,t_south,p_south);
                    end
                end
                
                tmp = pdns_south; pdns_south = nan*ones(size(p));
                pdns_south(1:n0) = tmp;
                
                %                        pdns_south = change(pdns_south,'<=',-99,nan);
                
                if check_cals == 1
                    disp('SOUTHERN CAST DEPTH_NTP COMPUTATIONTP')
                    cast_vals = [pdns_south(:)]
                end
                
                
                %%          intersection vitals
                
                for k0 = 1:n0
                    if abs(pdns_south(k0)+99)<=1e-2
                        k_south(k0,j0,i0) = nan;
                        r_south(k0,j0,i0) = nan;
                    elseif abs(pdns_south(k0)+99.1)<=1e-2
                        k_south(k0,j0,i0) = 1;
                        r_south(k0,j0,i0) = nan;
                    elseif abs(pdns_south(k0)+99.2)<=1e-2
                        k_south(k0,j0,i0) = n0;
                        r_south(k0,j0,i0) = nan;
                    elseif abs(pdns_south(k0)+99.3)<=1e-2
                        k_south(k0,j0,i0) = 0;
                        r_south(k0,j0,i0) = nan;
                    else
                        k_south(k0,j0,i0) = indx(p(1:n0_south,j0_south,i0),pdns_south(k0));
                        if pdns_south(k0)==p(n0_south,j0_south,i0)
                            k_south(k0,j0,i0) = n0_south-1;
                        end
                        kk = k_south(k0,j0,i0);
                        r_south(k0,j0,i0) = (pdns_south(k0)-p(kk,j0_south,i0))/...
                            (p(kk+1,j0_south,i0)-p(kk,j0_south,i0));
                        
                        if handles.quad==1
                            
                            r_ntp = r_south(k0,j0,i0); %[k0,j0,i0,r_ntp]
                            
                            pmid = (p(kk+1,j0_south,i0)+p(kk,j0_south,i0))/2;
                            
                            if strcmp(handles.eos,'eos80_t')
                                ct_u = ct_from_t(s(kk,j0_south,i0),t(kk,j0_south,i0),p(kk,j0_south,i0));
                                %                        [rho_u,rho_u_s,rho_u_t,rho_u_p] = eosall_from_ct(s(kk,j0_south,i0),ct,p(kk,j0_south,i0));
                                [rho_u,rho_u_s,rho_u_t,rho_u_p] = eosall_from_ct(s(kk,j0_south,i0),ct_u,pmid);
                                %                        rho_u_local = rho_from_ct(s(kk,j0_south,i0),ct,pmid);
                            else
                                'not yet implemented in intersections_south'
                            end
                            alfa_u = -rho_u_t/rho_u; beta_u = rho_u_s/rho_u;
                            if strcmp(handles.eos,'eos80_t')
                                ct_l = ct_from_t(s(kk+1,j0_south,i0),t(kk+1,j0_south,i0),p(kk+1,j0_south,i0));
                                %                        [rho_l,rho_l_s,rho_l_t,rho_l_p] = eosall_from_ct(s(kk+1,j0_south,i0),ct,p(kk+1,j0_south,i0));
                                [rho_l,rho_l_s,rho_l_t,rho_l_p] = eosall_from_ct(s(kk+1,j0_south,i0),ct_l,pmid);
                                %                        rho_l_local = rho_from_ct(s(kk+1,j0_south,i0),ct,pmid);
                            else
                                'not yet implemented in intersections_south'
                            end
                            
                            alfa_l = -rho_l_t/rho_l; beta_l = rho_l_s/rho_l;
                            term1 = (alfa_u+r_ntp*(alfa_l-alfa_u)/2)*(ct_l-ct_u)/(rho_l-rho_u);
                            term2 = (beta_u+r_ntp*(beta_l-beta_u)/2)*(s(kk+1,j0_south,i0)-s(kk,j0_south,i0))/(rho_l-rho_u);
                            r_south(k0,j0,i0) = r_ntp*(rho_l+rho_u)*(-term1+term2)/2; %r_quad = r_south(k0,j0,i0), dj_pause(0)
                            
                        end
                        
                    end
                    
                    %%       switch off between oceans
                    
                    oce_test = ocean(j0,i0)*ocean(j0_south,i0);
                    
                    if oce_test==3 | oce_test==5
                        k_south(k0,j0,i0) = -1;
                        r_south(k0,j0,i0) = nan;
                    end
                end
            end
        end
    end
    
    
    %%			ints_south plot
    
    inds_write = find(isfinite(k_south));
    if length(inds_write)>0
        counter = counter+1;
        if mod(counter,nwrite)==0
            na = length(find(isfinite(r_south(:,1:j0,:))));
            ns = length(find(isfinite(s(:,1:j0,:))));
            ok = [j0,100*j0/ny,100*na/ns]
            zz = squeeze(k_south(nz2,:,:));
            figure(4),clf, [hax,hax_noxlabel,hax_noylabel,hax_left,hax_bottom] = nfigaxes([1 2],[0.08 0.08],[0.08 0.92],[0.08 0.92]);
            axes(hax(1)) %subplot(2,2,3)
            dj_pltmp(longs(1,:),lats(:,1),zz,0)
            inds_3 = find(zz==0);
            if length(inds_3)>0
                [j_triples,i_triples]= ind2sub([ny,nx],inds_3);
                hold on
                plot(longss(i_triples),latss(j_triples),'m*')
                hold off
            end
            title('k\_south')
            set(gcf,'Name','South intersections','NumberTitle','off','Color',[0.961 0.988 0.965])
            
            zz = squeeze(r_south(nz2,:,:));
            figure(4), axes(hax(2))%subplot(2,2,4)
            dj_pltmp(longs(1,:),lats(:,1),zz,0)
            title('r\_south')
            figure(4); dj_pause(1);
        end
    end
    
end


%%			and save final intersections_south file

zz = squeeze(k_south(nz2,:,:));
figure(4), axes(hax(1))%subplot(2,2,3)
hold on
dj_pltmp(longs(1,:),lats(:,1),zz,0)
title('k\_south')

zz = squeeze(r_south(nz2,:,:));
figure(4), axes(hax(2))%subplot(2,2,4)
hold on
dj_pltmp(longs(1,:),lats(:,1),zz,0)
title('r\_south')


%%          triples

inds_3 = find(k_south==0);
if length(inds_3)>0
    [k_triples,j_triples,i_triples]= ind2sub([nz,ny,nx],inds_3);
    hold on
    plot(longss(i_triples),latss(j_triples),'m*')
    hold off
    figure(4)
end

figure(4); dj_pause(1);


return
