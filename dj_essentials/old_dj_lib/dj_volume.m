function volume = dj_volume(longs,lats,p1,surface)

%%%    dj_volume		volume of a tracer above a surface
%%%
%%%    Usage:			   volume = dj_colume(longs,lats,p,c,surface)
%%%
%%%    Input:        longs   	- vector of longitudes
%%%                  lats    	- vector of latitudes
%%%                  p			  - vector of pressures
%%%                  surface	- 2d (ny,nx) surface array
%%%
%%%    Output:       volume 	- volume above surface
%%%
%%%    Author:       David Jackett
%%%
%%%    Date:         17/12/98
%%%

global s th p

jmod = 8; method = 1; pr0 = 0; g = 9.8; ifig = 0;

nx = length(longs); ny = length(lats); nz = length(p1);

dx = 111.2e3*mean(diff(longs)); dy = 111.2e3*mean(diff(lats));


volume = 0; oo = ones(nz,1);

vlm = NaN*ones(ny,nx); 

for j = 1:ny

	 dA = dx*cos(pi*lats(j)/180)*dy;   

   inds = find(~isnan(th(1,j,:))); n = length(inds);
   
   if n > 0,     
      for i = 1:n
      	phi = surface(j,inds(i));
			  if finite(phi),
					 if method == 1,
		         vlm(j,inds(i)) = dj_quadi(p1(:),oo(:),p1(1),phi)*dA;
		         volume = volume+vlm(j,inds(i));					 
					 elseif method == 2,
						 ss = s(:,j,inds(i)); pr0 = 0;
						 tt = sw_ptmp(ss,th(:,j,inds(i)),pr0,p1); 
						 oo = ones(size(tt))./sw_dens(ss,tt,p1);
		         vlm(j,inds(i)) = 1.e4*dj_quadi(p1(:),oo(:),p1(1),phi)*dA/g;
		         volume = volume+vlm(j,inds(i));
					 end
			  end      
      end
   end
               
	if ifig == 1 & (mod(j,jmod) == 0 | j == ny)
		h = gcf; dj_pltmp(longs,lats,vlm);
    dj_pause(1); dj_toc; close(h);        
  end

end

return

