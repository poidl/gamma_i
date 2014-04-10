

%      the ocean array

dims = size(s);
dim_lat = size(lats);
dim_long = size(longs);

if length(dims)==3
    if length(dim_lat) == 1 & length(dim_long) == 1
        ocean = dj_ocean(longs,lats,squeeze(s(1,:,:)));
    elseif length(dim_lat) == 1 & length(dim_long) == 2
        ocean = dj_ocean(squeeze(longs(1,:)),squeeze(lats(:,1)),squeeze(s(1,:,:)));
    elseif length(dim_lat) == 1 & length(dim_long) == 1
        ocean = dj_ocean(squeeze(longs(1,1,:)),squeeze(lats(1,1,:)),squeeze(s(1,:,:)));
    else
        lat = unique(sort(lats(:)));
        long = unique(sort(longs(:)));
        ocean = dj_ocean(long,lat,squeeze(s(1,:,:)));
    end
    
    figure
    set(gcf,'Name','ocean','NumberTitle','off')
	dj_pltmp(longs,lats,ocean)
    title('ocean')
%     fpcolor(ocean)
elseif length(dims)==2 && dims(2)>1
    nx = dims(2);
    ocean = nan*ones(1,nx);
    for i = 1:nx
        ocean(i) = dj_ocean0(longs(i),lats(i));
    end
else
    error('****    not a section or ocean of data    ****')
end


%      the n array

n = ones(size(s)); 
stp = s.*t.*p;
Inan = find(isnan(stp));
if ~isempty(Inan)
    n(Inan) = nan;
end
n = reshape(nansum(n),size(ocean));
if length(dims) == 3
    figure(2)
    set(gcf,'Name','ocean depth index','NumberTitle','off','Color',[0.961 0.988 0.965])
	dj_pltmp(longs,lats,n)
    title('ocean depth index')
%     fpcolor(n)
    colormap(flipud(jet(64)))
end
