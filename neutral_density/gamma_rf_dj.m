
%%      label with gamma_rf

dj_tic

inds = find(s<0); s(inds) = nan;

ss = s(inds); ctt = ct(inds); pp = p(inds); 

%cd gfunction
    gg = gfunc(ss,ctt);
%cd ..

g = nan*ones(size(s)); g(inds) = gg;

if is_a_section(s,ct) == 1
    g(:,:,1) = g(:,:,2); g(:,:,3) = g(:,:,2);
    ok = 'section fixed'
end



%%      and plot

smin = nanmin(s(:)); smax = nanmax(s(:)); sby = (smax-smin)/100;

ctmin = nanmin(ct(:)); ctmax = nanmax(ct(:)); ctby = (ctmax-ctmin)/100;

s0 = smin:sby:smax; ct0 = ctmin:ctby:ctmax; ns = length(s0); nct = length(ct0);

[ss,ctt] = meshgrid(s0,ct0); ss = ss(:); ctt = ctt(:);

cd gfunction
    g0 = gfunc(ss,ctt); g0 = reshape(g0,size(ss));
cd ..

g0 = reshape(g0,nct,ns); s0, ct0, g0

figure(5)
    fpcolor(s0,ct0,g0), colorbar
    inds = find(isfinite(s(:))); ss = s(inds); ctt = ct(inds);
    hold on, plot(ss,ctt,'w.')
    contour(s0,ct0,g0,20,'k'), hold off
    


