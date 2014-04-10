function [Ibg,gamma_bdry] = gamma_boundary_gammas(gamma,long,lat)

[nz,ny,nx] = size(gamma); 
    
long3 = reshape(ones(nz,1)*long(:)',nz,ny,nx);
lat3 = reshape(ones(nz,1)*lat(:)',nz,ny,nx);

% Ig = find(isfinite(g));
% ng = length(Ig);
Ibg = find(isfinite(gamma) & 170<=long3 & long3<=270 & -10<=lat3 & lat3<=10);
% n_bg = length(Ibg);
gamma_bdry = gamma(Ibg);

% if n_bg<=20
%     Ibg_w = [1:n_bg];
% else
%     Ibg_w = [1:floor(n_bg/10):n_bg];
% end

%some_boundary_gammas = [Ibg(Ibg_w),g_bdry(Ibg_w)];

%boundary_percent_of_data = 100*n_bg/ng;

% figure
% longss = nanmean(long);
% latss = nanmean(lat');
% dj_pltmp(longss,latss,n)
% hold on
% plot(longs3(Ibg(:)),lats3(Ibg(:)),'m*')
% hold off
%clear longs3 lats3

%dj_pause(1)

end