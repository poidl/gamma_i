global s ct p g ocean n longs lats 

global h_numerator 

smin = 20; smax = 42; sby = (smax-smin)/100;

ctmin = -5; ctmax = 40; ctby = (ctmax-ctmin)/100;

s0 = smin:sby:smax; ct0 = ctmin:ctby:ctmax; ns = length(s0); nct = length(ct0)

[ss,ctt] = meshgrid(s0,ct0); ss = ss(:); ctt = ctt(:);

cmd = ['copyfile(''gamma_pp', h_numerator, '.dat'', ''gamma_pp.dat'')']
eval(cmd)

load gamma_pp.dat

n_numerator = eval(h_numerator);

nregions = length(gamma_pp)/n_numerator; figure

for k = 1:nregions
    
    coeffs = gamma_pp((k-1)*n_numerator+1: k*n_numerator);
    
    cmd = ['g0 = gamma_p', h_numerator, '(ss,ctt,coeffs);'];
    eval(cmd)
    
    g0 = reshape(g0,nct,ns);
    
    fpcolor(s0,ct0,g0), colorbar
    inds = find(finite(s(:))); sw = s(inds); ctw = ct(inds);
    hold on, plot(sw,ctw,'w.')
    contour(s0,ct0,g0,20,'k')
    
    pp0 = zeros(size(ss)); sig00 = rho_from_ct(ss,ctt,pp0)-1000;

    g0 = reshape(g0,nct,ns); sig0 = reshape(sig00,nct,ns);
    
    contour(s0,ct0,sig0,20,'m'),hold off
    
    if k==1
        region = 'NA north';
    elseif k==2
        region = 'NA south';
    elseif k==3
        region = 'SA west';
    elseif k==4
        region = 'SAI';
    end
    
    title(['region ', num2str(k), '   ', region])
    
    figure(gcf), dj_pause(0)

end


