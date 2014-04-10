factor = 0; figure(gcf);

z = colormap; r = z(:,1); green = z(:,2); b = z(:,3); nr = length(r);

smid = 0; soln = w0;

s0 = nanmin(nanmin(soln)); s1 = nanmax(nanmax(soln));

ds1 = (smid-s0)/31;

salt = s0:ds1:smid;

v1 = (salt-s0)/(s1-s0);

ds1 = (s1-smid)/32;

salt = smid+ds1:ds1:s1;

v2 = (salt-s0)/(s1-s0);

v = [v1,v2];

[cm] = interpcolormap(r, green, b, nr, v);

colormap(cm)

