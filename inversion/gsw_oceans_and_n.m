function [ocean, n]=gsw_oceans_and_n(SA,CT,p,long,lat)



% the ocean array
dim_SA = size(SA);
dim_lat = size(lat);
dim_long = size(long);

if length(dim_SA) == 3
    if length(dim_lat) == 1 & length(dim_long) == 1
        ocean = gsw_ocean_3d(long,lat,squeeze(SA(1,:,:)));
    elseif length(dim_lat) == 1 & length(dim_long) == 2
        ocean = gsw_ocean_3d(squeeze(long(1,:)),squeeze(lat(:,1)),squeeze(SA(1,:,:)));
    elseif length(dim_lat) == 1 & length(dim_long) == 1
        ocean = gsw_ocean_3d(squeeze(long(1,1,:)),squeeze(lat(1,1,:)),squeeze(SA(1,:,:)));
    else
        ocean = gsw_ocean_3d(unique(sort(long(:))),unique(sort(lat(:))),squeeze(SA(1,:,:)));
    end 
elseif length(dim_SA)==2 && dim_SA(2)>1
    nx = dim_SA(2);
    ocean = nan(1,nx);
    for i = 1:nx
        ocean(i) = gsw_ocean_2d(long(i),lat(i));
    end
else
    error('****    not a section or ocean of data    ****')
end

% the n array
n = ones(size(SA)); 
data = SA.*CT.*p;
Inan = find(isnan(data));
if ~isempty(Inan)
    n(Inan) = nan;
end
n = reshape(nansum(n),size(ocean));

end
