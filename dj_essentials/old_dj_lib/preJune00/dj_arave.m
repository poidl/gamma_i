function ave_value = dj_arave(longs,lats,z)

%%%    dj_arave		area average value
%%%
%%%    Usage:			ave = dj_arave(longs,lats,z)
%%%
%%%    Input:        longs - longitudes of observations (nx)
%%%                  lats  - latitudes of observations (ny)
%%%                  z     - matrix (ny-by-nx)
%%%
%%%    Output:       ave_value   - area average value of array z
%%%
%%%    Author:       David Jackett
%%%
%%%    Date:         17/10/96
%%%


nx = length(longs); ny = length(lats);

[ny_z,nx_z] = size(z);

if nx~=nx_z|ny~=ny_z, error('ERROR in function dj_arave'); end
	
dx = mean(diff(longs)); dy = mean(diff(lats));

areas = dx*cos(pi*lats/180)*dy;

areas = areas(:)*ones(1,nx);

inds = find(isnan(z)); areas(inds) = nan*ones(size(inds));

ave_value = nansum(nansum(z'.*areas'))/nansum(nansum(areas'));



%fpcolor(longs,lats,areas); colormap(jet); colorbar('horizontal')

%disp(['mean surface temperature: ',num2str(ave_value),' degrees C'])

%figure(gcf); dj_pause(0)


return

