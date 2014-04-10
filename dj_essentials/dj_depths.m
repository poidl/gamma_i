function depths = dj_depths(s,t,p,longs,lats)

[nz,ny,nx] = size(s);
depths = reshape(p*ones(1,nx*ny),nz,ny,nx);
inds = find(isnan(s)); depths(inds) = nan;
depths = reshape(nanmax(depths),ny,nx);
dj_pltmp(longs,lats,depths)

zoom on

return