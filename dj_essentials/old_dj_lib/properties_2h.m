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

global s1 th1 p1


%%									identical code to dj_burden.m


jmod = 12; method = 1; pr0 = 0; g = 9.8;

nx = length(longs); ny = length(lats); nz = length(p);

[nz_c,ny_c,nx_c] = size(c);

if nx~=nx_c | ny~=ny_c | nz~=nz_c,
				 error('ERROR in function dj_burden'); end


dx = 111.2e3*mean(diff(longs)); dy = 111.2e3*mean(diff(lats));


burden = 0; volume = 0; oo = ones(nz,1);

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
					 if method == 1,
		         vlm(j,inds(i)) = dj_quadi(p(:),oo(:),p(1),phi)*dA;
		         volume = volume+vlm(j,inds(i));					 
					 elseif method == 2,
						 ss = s1(:,j,inds(i)); pr0 = 0;
						 tt = sw_ptmp(ss,th1(:,j,inds(i)),pr0,p); 
						 oo = ones(size(tt))./sw_dens(ss,tt,p);
		         vlm(j,inds(i)) = 1.e4*dj_quadi(p(:),oo(:),p(1),phi)*dA/g;
		         volume = volume+vlm(j,inds(i));
					 end
			  end      
      end
   end
               
	if mod(j,jmod) == 0 | j == ny
		h = gcf; zz = brdn./change(vlm,'==',0,eps);
		if volume~=0, cave = burden/volume;  else, cave = 0; end     
    dj_pltmp(longs,lats,zz);
		legend = ['volume averaged c: ', num2str(cave,3)];
		title(legend)       
    dj_pause(1); dj_toc; close(h);        
  end

end


%% 								now split into hemispheres


inds = find(lats>=0); ind0 = inds(1);

inds_n = ind0:ny; inds_s = 1:ny-1;

brdn_n = brdn; brdn_n(inds_s) = nan;
vlm_n = vlm; vlm_n(inds_s) = nan;

brdn_s = brdn_n;
vlm_s = vlm_n;

zz = brdn_n./vlm_n;

dj_pltmp(longs,lats,zz);
dj_pause(1); dj_toc;


return

