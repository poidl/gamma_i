function ocean = gsw_ocean_3d(long,lat,z)

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
%%%                          1    main oceans,
%%%                          7    Arctic
%%%                          8    Med
%%%
%%%    Author:    	David Jackett
%%%
%%%    Date:      	7/7/98
%%%

nx = length(long); 
ny = length(lat); 
ocean = nan(ny,nx);

for j = 1:ny
   for i = 1:nx
      if ~isnan(z(j,i))
        ocean(j,i) = dj_ocean_2d(long(i),lat(j));
      end   
   end
end

end