figure

%%      the ocean array
tic
ocean = dj_ocean(longs,lats,squeeze(s(1,:,:)));
toc
dj_pltmp(longs,lats,ocean)

%%      the n array

dj_tic
    n = ones(size(s));
    [nz,ny,nx] = size(s);
    inds = find(isnan(s)); n(inds) = nan;
    n = reshape(nansum(n),ny,nx);
    figure
    dj_pltmp(longs,lats,n)
    dj_toc

