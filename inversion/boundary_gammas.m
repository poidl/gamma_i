
[nz,ny,nx] = size(g); 
    
longs3 = reshape(ones(nz,1)*longs(:)',nz,ny,nx);

lats3 = reshape(ones(nz,1)*lats(:)',nz,ny,nx);

inds_g = find(isfinite(g)); ng = length(inds_g);

inds_bg = find(isfinite(g)&170<=longs3&longs3<=240&-10<=lats3&lats3<=10); n_bg = length(inds_bg);

g_bdry = g(inds_bg);

if n_bg<=20
    inds_bg_w = 1:n_bg;
else
    inds_bg_w = 1:floor(n_bg/10):n_bg;
end

some_boundary_gammas = [inds_bg(inds_bg_w),g_bdry(inds_bg_w)]

boundary_percent_of_data = 100*n_bg/ng

figure

longss = nanmean(longs); latss = nanmean(lats');

dj_pltmp(longss,latss,n), hold on

plot(longs3(inds_bg(:)),lats3(inds_bg(:)),'m*'), hold off

clear longs3 lats3

dj_pause(1)