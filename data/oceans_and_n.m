

%%      the ocean array

dims = size(s);

if length(dims)==3
    ocean = dj_ocean(longs(1,:),lats(:,1),squeeze(s(1,:,:)));
    figure, set(gcf,'Name','ocean','NumberTitle','off')
	dj_pltmp(longs(1,:),lats(:,1),ocean), title('ocean')
%     fpcolor(ocean)
elseif length(dims)==2 && dims(2)>1
    nx = dims(2); ocean = nan*ones(1,nx);
    for i = 1:nx, ocean(i) = dj_ocean0(longs(i),lats(i)); end
else
    error('****    not a section or ocean of data    ****')
end



%%      the n array

n = ones(size(s)); 
inds = find(isnan(s)); n(inds) = nan;
n = reshape(nansum(n),size(ocean));
if length(dims)==3
    figure(2)
    set(gcf,'Name','ocean depth index','NumberTitle','off','Color',[0.961 0.988 0.965])
	dj_pltmp(longs(1,:),lats(:,1),n), title('ocean depth index')
%     fpcolor(n)
    colormap(flipud(jet(64)))
end
