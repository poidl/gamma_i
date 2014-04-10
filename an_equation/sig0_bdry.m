global s ct t p g ocean n longs lats o3 lats3

global h_numerator h_denominator h_normalise

cmd = ['copyfile(''rfunc_', h_numerator, '_', h_denominator, '.dat'', ''rfunc.dat'')']

eval(cmd)

load rfunc.dat

s0 = 0:0.1:42; ct0 = -5:0.2:40; ns = length(s0); nct = length(ct0);

[ss,ctt] = meshgrid(s0,ct0);  ss = ss(:); ctt = ctt(:);

pp = zeros(size(ss)); sig00 = rho_from_ct(ss,ctt,pp)-1000;

if h_normalise==1
    ss = ss/40; ctt = ctt/30;
end

cmd = ['g0 = rfunc_', h_numerator, '_', h_denominator, '(ss,ctt,rfunc);'];, eval(cmd)

g0 = reshape(g0,nct,ns); sig0 = reshape(sig00,nct,ns);


if h_normalise==1
    g0 = 30*g0;
end


figure

    fpcolor(s0,ct0,g0), colorbar
    
    inds = find(finite(s(:))); ss = s(inds); ctt = ct(inds); pp = p(inds); gg = g(inds);
    
    ctt_max = max(ctt);
    
    hold on, plot(ss,ctt,'w.')
    
    contour(s0,ct0,g0,20,'k'), %contour(s0,ct0,sig0,20,'w'),hold off
    
    [gg_srtd, isrtd] = sort(gg); ng = length(gg);
    
    n100 = floor(ng/100); 
    
    sss = ss(isrtd(1:n100)); cttt = ctt(isrtd(1:n100)); ggg = gg(isrtd(1:n100));
    
    inds1 = find(pp(isrtd(1:n100))==0);
    
    hold on, plot(sss,cttt,'k.'), plot(sss(inds1),cttt(inds1),'b.'), hold off
    
    ss(isrtd(1)),ctt(isrtd(1)), pp(isrtd(1)), max(pp(isrtd(1:n100)))
    
    
    
%%          add leftmost and upper gammas

s0 = 0:1:42; ct0 = -5:1:40; ns = length(s0); nct = length(ct0);

[ss,ctt] = meshgrid(s0,ct0);  ss = ss(:); ctt = ctt(:);

pp = zeros(size(ss)); sig00 = rho_from_ct(ss,ctt,pp)-1000;

sig0 = rho_from_ct(sss(1)-1,cttt(1),0)-1000



inds = find(sig00>=sig0&ctt<ctt_max+3); ss(inds) = nan; ctt(inds) = nan; pp(inds) = nan;

hold on, plot(ss,ctt,'m.'), hold off

inds = find(finite(ss)); ss_sigbdry = ss(inds); ctt_sigbdry = ctt(inds); sig00 = sig00(inds);

gg_sigbdry = ggg(1) - (sig0-sig00);

no_sig0_bdry = length(ss_sigbdry)

save sig0_bdry ss_sigbdry ctt_sigbdry gg_sigbdry

    
    
