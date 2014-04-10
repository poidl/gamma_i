function [b_n,a_n,b_s,a_s] = dj_burden_2hA(longs,lats,p,c,s,slim)

%%%    dj_burden2        2D burden (and area) of a tracer above an s limit
%%%
%%%    Usage:			[burden,volume] = dj_burden2(longs,lats,p,c,s,slim)
%%%
%%%    Input:        longs   	- vector of longitudes
%%%                  lats    	- vector of latitudes
%%%                  p        - vector of pressures
%%%                  c       	- 2D (ny,nx) tracer array
%%%                  s        - surface s distribution
%%%                  slim   	- surface s limit
%%%
%%%    Output:       burden   - area integral of c above limit
%%%                  area	    - surface area above s limit
%%%
%%%    Author:       David Jackett
%%%
%%%    Date:         25/11/98
%%%

global brdn

jmod = 124; ifig = 0;

nx = length(longs); ny = length(lats);

[ny_c,nx_c] = size(c);

if nx~=nx_c | ny~=ny_c,
				 error('ERROR in function dj_burden'); end


dx = 111.2e3*mean(diff(longs)); dy = 111.2e3*mean(diff(lats));


burden = 0; area = 0; 

brdn = NaN*ones(ny,nx); ar = NaN*ones(ny,nx); 

for j = 1:ny
   
	dA = dx*cos(pi*lats(j)/180)*dy;

   inds = find(s(j,:)<=slim); n = length(inds);
   
   if n > 0,     
      for i = 1:n
      	cc = c(j,inds(i));
			brdn(j,inds(i)) = cc*dA;
         burden = burden+brdn(j,inds(i));
		   ar(j,inds(i)) = dA;
		   area = area+dA;
		end      
   end
               
	if ifig == 1 &  j == ny  | mod(j,jmod) == 0
		h = gcf; zz = NaN*ones(ny,nx);
		if area~=0
			zz = brdn./ar; cave = burden/area;
		else
  			cave = 0;
	 	end     
    	dj_pltmp(longs,lats,zz);
		legend = ['area averaged c: ', num2str(cave,3)];
		title(legend)       
    	dj_toc; dj_pause(1); close(h);        
  	end

end



%% 								now split into hemispheres


inds = find(lats>=0); ind0 = inds(1);

inds_n = ind0:ny; inds_s = 1:ind0-1;

brdn_n = brdn; brdn_n(inds_s,:) = nan; b_n = nansum(nansum(brdn_n));
ar_n = ar; ar_n(inds_s,:) = nan; a_n = nansum(nansum(ar_n));

brdn_s = brdn; brdn_s(inds_n,:) = nan; b_s = nansum(nansum(brdn_s));
ar_s = ar; ar_s(inds_n,:) = nan; a_s = nansum(nansum(ar_s));

if ifig == 2
   
	zz = brdn_n./change(ar_n,'==',0,eps);
	dj_pltmp(longs,lats,zz);

	zz = brdn_s./change(ar_s,'==',0,eps);
	dj_pltmp(longs,lats,zz);

	dj_pause(1)

end


return

