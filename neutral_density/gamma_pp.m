
%%      label with gamma_pp

indss = find(finite(s));

if h_normalise==1
    ss = s(indss)/40; ctt = ct(indss)/30; pp = p(indss); 
else
    ss = s(indss); ctt = ct(indss); pp = p(indss); 
end
    
h_normalise

cd d:/neutrals/ness8/an_equation
    cmd = ['copyfile(''rfunc_', num2str(h_numerator), '_', num2str(h_denominator), '.dat'', ''rfunc.dat'')']
    eval(cmd)
    load rfunc.dat
    rfunc
    cmd = ['g_pr = rfunc_', num2str(h_numerator), '_', num2str(h_denominator), '(ss,ctt,rfunc);']
    eval(cmd)
cd ../neutral_density
        
g(indss) = g_pr;

if h_normalise==1, g = 30*g; end

if mean(g_pr)>1000, g = g-1000; end


%%      and plot

smin = nanmin(s(:)); smax = nanmax(s(:)); sby = (smax-smin)/100;

ctmin = nanmin(ct(:)); ctmax = nanmax(ct(:)); ctby = (ctmax-ctmin)/100;

s0 = smin:sby:smax; ct0 = ctmin:ctby:ctmax; ns = length(s0); nct = length(ct0);

[ss,ctt] = meshgrid(s0,ct0); ss = ss(:); ctt = ctt(:);

if h_normalise==1
    ss = ss/40; ctt = ctt/30;
end

cd ../an_equation
cmd = ['g0 = rfunc_', h_numerator, '_', h_denominator, '(ss,ctt,rfunc);'];, eval(cmd)
cd ../neutral_density

g0 = reshape(g0,nct,ns);

if h_normalise==1
    g0 = 30*g0;
end

figure(5)
    fpcolor(s0,ct0,g0), colorbar
    inds = find(finite(s(:))); ss = s(inds); ctt = ct(inds);
    hold on, plot(ss,ctt,'w.')
    contour(s0,ct0,g0,20,'k'), hold off

dj_pause(1)



