function [b_n,v_n,b_s,v_s] = dj_burden_2h(longs,lats,p,c,surface)

%%%    dj_burden		burden (and volume) of a tracer above a surface
%%%
%%%    Usage:			   [burden,volume] = dj_burden(ongs,lats,p,c,surface)
%%%
%%%    Input:        longs   	- vector of longitudes
%%%                  lats    	- vector of latitudes
%%%                  p			  - vector of pressures
%%%                  c       	- 3D (nz,ny,nx) tracer array
%%%                  surface	- 2d (ny,nx) surface array
%%%
%%%    Output:       burden   - area integral of c above surface
%%%                  volume 	- volume above surface
%%%
%%%    Author:       David Jackett
%%%
%%%    Date:         8/7/98
%%%

global brdn0 vlm0 brdn vlm



%%									identical code to dj_burden.m


ifig = 1; jmod = 112;

nx = length(longs); ny = length(lats); nz = length(p);

[nz_c,ny_c,nx_c] = size(c);

if nx~=nx_c | ny~=ny_c | nz~=nz_c,
				 error('ERROR in function dj_burden'); end


dx = 111.2e3*mean(diff(longs)); dy = 111.2e3*mean(diff(lats));


burden = 0; volume = 0;

brdn = NaN*ones(ny,nx); vlm = NaN*ones(ny,nx); 

for j = 1:ny

	 dA = dx*cos(pi*lats(j)/180)*dy;   

   inds = find(~isnan(c(1,j,:))); n = length(inds);
   
   if n > 0,     
      for i = 1:n
      	cc = c(:,j,inds(i)); phi = surface(j,inds(i));
			  if finite(phi),
			     brdn(j,inds(i)) = dj_quadi(p(:),cc(:),p(1),phi)*dA;
           burden = burden+brdn(j,inds(i));
		       vlm(j,inds(i)) = (phi-p(1))*dA;
		       volume = volume+vlm(j,inds(i));					 
					 end      
      end
   end
               
	if ifig==1 & (mod(j,jmod) == 0 | j == ny)
		h = gcf; zz = (brdn-brdn0)./change((vlm-vlm0),'==',0,eps);
		if volume~=0,
			 cave = nansum(nansum(brdn-brdn0))/nansum(nansum(vlm-vlm0));
	  else,
			 cave = 0;
	  end     
    dj_pltmp(longs,lats,zz);
		legend = ['volume averaged c: ', num2str(cave,3)];
		title(legend)       
    dj_pause(1); dj_toc; close(h);        
  end

end


%% 								now split into hemispheres


inds = find(lats>=0); ind0 = inds(1);

inds_n = ind0:ny; inds_s = 1:ind0-1;

brdn_n = brdn; brdn_n(inds_s,:) = nan; b_n = nansum(nansum(brdn_n));
vlm_n = vlm; vlm_n(inds_s,:) = nan; v_n = nansum(nansum(vlm_n));

brdn_s = brdn; brdn_s(inds_n,:) = nan; b_s = nansum(nansum(brdn_s));
vlm_s = vlm; vlm_s(inds_n,:) = nan; v_s = nansum(nansum(vlm_s));


%zz = brdn_n./change(vlm_n,'==',0,eps);
%dj_pltmp(longs,lats,zz);

%zz = brdn_s./change(vlm_s,'==',0,eps);
%dj_pltmp(longs,lats,zz);

%dj_toc; dj_pause(5)


return

