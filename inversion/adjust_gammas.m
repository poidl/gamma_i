%  get lower and upper medians
load intersections/gamma_medians

median_gammas_required = gi_adj(1:2);
gmu = gi_adj(3);

Idata = find(isfinite(g));
gg = g(Idata); 
nn = length(Idata);

[gg_srtd, kk_srtd] = sort(gg); 

nn05 = round(nn*0.05); 
nn95 = round(nn*0.95); 

median_gammas = gg_srtd([nn05; nn95])

sums_of_squares = sum((median_gammas_required-median_gammas).^2)

%  the linear transformation
g0bar = median_gammas_required(1);
g1bar = median_gammas_required(2);

g0 = median_gammas(1);
g1 = median_gammas(2);

ggbar = nanmean(g(:)); 

alpha_num = g0*g0bar + g1*g1bar - g0*ggbar - g1*ggbar - ...
                            g0bar*gmu - g1bar*gmu + 2*ggbar*gmu;
                        
alpha_den = g0^2 + g1^2 - 2*g0*gmu - 2*g1*gmu + 2*gmu^2;

alpha = alpha_num/alpha_den

beta = ggbar - alfa*gmu


%  make the adjustment and check
g0 = g; 
g = alfa*g + beta;

gmeans = [nanmean(g0(:)), nanmean(g(:))]

Idata = find(isfinite(g));
gg = g(Idata);
nn = length(Idata);

[gg_srtd, kk_srtd] = sort(gg); 

nn05 = round(nn*0.05);
nn95 = round(nn*0.95); 

median_gammas_achieved = gg_srtd([nn05; nn95])

sums_of_squares = sum((median_gammas_required - median_gammas_achieved).^2)
