outcrop = reshape(burden(1,:),ny,nx);
dj_pltmp(longs,lats,outcrop);
dj_pause(0)

for k = 1:nlevels

	inds = find(~isnan(burden(k,:)));
	outcrop = reshape(burden(k+1,:),ny,nx);
	outcrop(inds) = NaN;
	dj_pltmp(longs,lats,outcrop);

	dj_pause(0)

end