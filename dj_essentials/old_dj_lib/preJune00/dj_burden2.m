function [burden,area] = dj_burden2(longs,lats,p,c,s,slim)

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

jmod = 4; ifig = 0;

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
               
	if ifig == 1 &  j == ny  %| mod(j,jmod) == 0
 		h = gcf; zz = NaN*ones(ny,nx);
		if area~=0, zz = brdn./ar; cave = burden/area;
	  else  			cave = 0;
	  end     
    dj_pltmp(longs,lats,zz);
		legend = ['area averaged c: ', num2str(cave,3)]; title(legend)       
    dj_toc; dj_pause(1); close(h);        
  end

end


return

