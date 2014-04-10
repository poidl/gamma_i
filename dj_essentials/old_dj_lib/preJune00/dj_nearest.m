function  inds = dj_nearest(longs,lats,long0,lat0,ocean,n)

%%%    dj_nearest       nearest neighbour in longs/lats space
%%%
%%%    Usage:        inds = dj_arave(longs,lats,long0,lat0)
%%%
%%%    Input:        longs   - longitudes [0,360]
%%%                  lats    - latitudes of observations [-90,90]
%%%                  long0   - longitude
%%%                  lat0    - latitude
%%%                  c       - 2D tracer array
%%%                  n       - # of closest points required
%%%
%%%    Output:       inds    - indices, [1,nx*ny]
%%%
%%%    Author:       David Jackett
%%%
%%%    Date:         10/12/98
%%%


nx = length(longs); ny = length(lats);

[LONGS,LATS] = meshgrid(longs,lats); LONGS0 = LONGS;

iinds = [1:nx]; jinds = [1:ny};

[IINDS,JINDS] = meshgrid(iinds,jinds);


inds = find(~isnan(c));

lo_alive = LONGS(inds); la_alive = LATS(inds);


dist = []; nn = length(lo_alive);

for i = 1:nn
	dx = abs(long0-lo_alive(i)); dy = abs(lat0-la_alive(i));
	dist0 = dx+dy
	dist = [dist,dist0];
end

[dist_srtd,i_srtd] = sort(dist);

ii = inds(i_srtd(1:10));

lo_alive(ii)

la_alive(ii)


return

