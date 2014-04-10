figure

%%      the ocean array

ocean = dj_ocean(longs,lats,squeeze(s(1,:,:)));

dj_pltmp(longs,lats,ocean)

%%      the n array

n = ones(size(s));
[nz,ny,nx] = size(s);
inds = find(isnan(s)); n(inds) = nan;
n = reshape(nansum(n),ny,nx);
figure
dj_pltmp(longs,lats,n)

