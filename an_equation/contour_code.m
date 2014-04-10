global s ct p g ocean n longs lats 

global h_numerator h_denominator 

if eval(h_denominator)~=0
    cmd = ['copyfile(''rfunc_', handles.numerator, '_', handles.denominator, '.dat'', ''gamma.dat'')']
else
    cmd = ['copyfile(''gamma_p', handles.numerator, '.dat'', ''gamma.dat'')']
end
eval(cmd)

load gamma.dat

smin = nanmin(s(:)); smax = nanmax(s(:)); sby = (smax-smin)/100;

ctmin = nanmin(ct(:)); ctmax = nanmax(ct(:)); ctby = (ctmax-ctmin)/100;


s0 = smin:sby:smax; ct0 = ctmin:ctby:ctmax; ns = length(s0); nct = length(ct0);

[ss,ctt] = meshgrid(s0,ct0); ss = ss(:); ctt = ctt(:);

if handles.normalise==1
    ss = ss/40; ctt = ctt/30;
end
   
if eval(h_denominator)~=0
    cmd = ['g0 = rfunc_', handles.numerator, '_', handles.denominator, '(ss,ctt,gamma);'];
else
    cmd = ['g0 = gamma_p', handles.numerator, '(ss,ctt,gamma);'];, 
end
eval(cmd)

g0 = reshape(g0,nct,ns);

if handles.normalise==1
    g0 = 30*g0;
end

if nanmean(g0)>1000
    g0 = g0-1000;
end

figure(4), subplot(2,2,1)
    fpcolor(s0,ct0,g0), colorbar
    inds = find(finite(s(:))); ss = s(inds); ctt = ct(inds);
    hold on, plot(ss,ctt,'w.')
    contour(s0,ct0,g0,20,'k'), hold off

dj_pause(1)



s0 = 0:0.1:42; ct0 = -5:0.2:40; ns = length(s0); nct = length(ct0);

[ss,ctt] = meshgrid(s0,ct0);  ss = ss(:); ctt = ctt(:);

pp = zeros(size(ss)); sig00 = rho_from_ct(ss,ctt,pp)-1000;

if handles.normalise==1
    ss = ss/40; ctt = ctt/30;
end

if eval(h_denominator)~=0
    cmd = ['g0 = rfunc_', handles.numerator, '_', handles.denominator, '(ss,ctt,gamma);'];
else
    cmd = ['g0 = gamma_p', handles.numerator, '(ss,ctt,gamma);'];, 
end
eval(cmd)

g0 = reshape(g0,nct,ns); sig0 = reshape(sig00,nct,ns);

if handles.normalise==1
    g0 = 30*g0;
end

if nanmean(g0)>1000
    g0 = g0-1000;
end

figure(4), subplot(2,2,4)
    fpcolor(s0,ct0,g0), colorbar
    inds = find(finite(s(:))); ss = s(inds); ctt = ct(inds);
    hold on, plot(ss,ctt,'w.')
    contour(s0,ct0,g0,20,'k'), contour(s0,ct0,sig0,20,'m'),hold off
    
    
    
