function ocean = dj_ocean(longs,lats,z)

%%%    dj_ocean      ocean of longs/lats observations
%%%
%%%    Usage:     	ocean = dj_ocean(longs,lats,z)
%%%
%%%    Input:     	longs - vector (nx) of longitudes [0,360] 
%%%               	lats  - vector (ny) of latitudes [-90,90] 
%%%               	z     - array (ny,nx) describing ocean extent
%%%
%%%    Output:    	ocean - array (ny,nx) of ocean values
%%%                          0    land,
%%%                          1-6  main oceans,
%%%                          >6   Arctic & marginal seas
%%%
%%%    Author:    	David Jackett
%%%
%%%    Date:      	7/7/98
%%%


nx = length(longs); ny = length(lats); ocean = NaN*ones(ny,nx);

for j = 1:ny
   
   for i = 1:nx
      if ~isnan(z(j,i))
			ocean(j,i) = dj_ocean0(longs(i),lats(j));
		end   
   end
   
   if mod(j,1) == 0 | j == ny  
      dj_pltmp(longs,lats,ocean,0);      
   	figure(gcf), dj_pause(1)       
   end

end

dj_toc


return
