function gamma_i = gamma_3d(SA,CT,p,long,lat)

% gamma_3d                             gamma 3d "Neutral Density" (gamma_i)
%==========================================================================
%
% USAGE:
%  gamma_i = gamma_3d(SA,CT,p,lat,long)
% 
% DESCRIPTION:
%  Calculates the Neutral Density of 3D gridded data, all points are 
%  labeled
%
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in IOC et al. (2013).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA   =  Absolute Salinity.                                      [ g/kg ]
%  CT   =  Conservative Temperature (ITS-90)                      [ deg C ]
%  p    =  sea pressure                                            [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%  long =  longitude
%  lat  =  latitude
%
% OUTPUT:
%  gamma_i  =  Neutral Density based on Jackett and McDougall     [ g/kg ]
%
% REFERENCES:
%
%  McDougall T.J., P.M. Barker, R. Feistel and D.R. Jackett, 2011:  A 
%   computationally efficient 48-term expression for the density of 
%   seawater in terms of Conservative Temperature, and related properties
%   of seawater.  To be submitted to Ocean Science Discussions. 
%
%==========================================================================



[ocean, n] = gamma_ocean_and_n(SA,CT,p,long,lat);

%Initial estimate of the neutral surface - it would be better to use a
%locally referenced density surface.
gamma_initial = gamma_rf(SA,CT); % (SA,CT,p,long,lat);

[nz,ny,nx] = size(SA);

% called_ggrads = 0;

handles.inv_meth = 2;
% if handles.inv_meth == 1
%     method = 'direct'
% else
    method = 'iterative';
% end

vertical_eqns = 0;
bxz_eqns = 0;
gn_eqns = 1;

handles.zbc = 0;
handles.quad = 0;
handles.linb = 0;
handles.h2vwt = 1;
handles.helwt = 1e-6;
handles.bdrywt = 1e-5;
handles.maxits = 100;
handles.ntp = 1;
handles.ns = 1;
handles.gfunc = 0;
handles.sigp = 0;

% load equation set-up information

%load gk_interp_int_east
[k_east,r_east] = gamma_intersections_east(SA,CT,p,ocean,n);

%load gk_interp_int_west
[k_west,r_west] = gamma_intersections_west(SA,CT,p,ocean,n);

%load gk_interp_int_north
[k_north,r_north] = gamma_intersections_north(SA,CT,p,ocean,n);

%load gk_interp_int_south
[k_south,r_south] = gamma_intersections_south(SA,CT,p,ocean,n);

%load gk_interp_gamma_boundary
[I_bg, gamma_bdry] = gamma_boundary_gammas(gamma_initial,long,lat);

% load helicities

%load gk_interp_vert_vitals
[k_vert,r_vert,gprod_vert,no_eqs] = gamma_vertical_vitals(SA,CT,p,gamma_initial,n,r_north,r_east,r_south,r_west);

% load bxz_equations

%==========================================================================
% change this to interpolate SA and CT, not N2.
% determine N2 for weighting
Idata = find(isfinite(SA(1,:)));
nn = length(Idata);
ss = SA(:,Idata);
tt = CT(:,Idata);
pp = p(:,Idata);
N2 = nan(size(SA));
for kk = 1:nn
    Iz = find(isfinite(ss(:,kk)));
    [N2_dummy, p_mid_dummy] = gsw_Nsquared(ss(Iz,kk),tt(Iz,kk),pp(Iz,kk));
    N2(Iz,Idata(kk)) = gsw_interp_Nsquared(N2_dummy,p_mid_dummy,pp(Iz,kk));
end
%==========================================================================

% determine # equations and # unknowns
I = find(isfinite(r_east));
no_h_eqns = length(I);

I = find(isfinite(r_west));
no_h_eqns = no_h_eqns + length(I);

I = find(isfinite(r_north));
no_h_eqns = no_h_eqns + length(I);

I = find(isfinite(r_south));
no_h_eqns = no_h_eqns + length(I);

no_v_eqns = 0;
if vertical_eqns == 1
    no_v_eqns = no_eqs;
end

no_b_eqns = length(I_bg);

no_eqns = no_h_eqns + no_v_eqns + no_b_eqns;

I_gamma = find(isfinite(gamma_initial));
ng = length(I_gamma);

no_equations = no_eqns;
no_unknowns = ng;

ginds = nan(size(SA));
ginds(I_gamma) = [1:ng];

gamma_old = gamma_initial(I_gamma);

g_old_min = min(gamma_old);
g_old_max = max(gamma_old);

kvert = 0;
kxzbeq = 0;

ns = 3*(no_h_eqns + no_v_eqns) + no_b_eqns;

neq = 0;
ieq = 0;
s1 = nan(ns,1);
s2 = s1;
s3 = s1;
b = zeros(no_eqns,1);

ks = nan(no_eqns,1);
js = ks;
is = ks;
wts = ks;

% set up lateral equations.
for kg = 1:ng
    [k0,j0,i0] = ind2sub([nz,ny,nx],I_gamma(kg));
    helwt = handles.helwt;
    kwt = 1;
    % wt = (1+p(k0,j0,i0)/kwt)^helwt;
    wt = 1/(helwt + N2(k0,j0,i0));
    wts(kg) = wt;
    ks(kg) = k0;
    js(kg) = j0;
    is(kg) = i0;
      
    % lateral equations
    if isfinite(r_east(I_gamma(kg)))        
        if i0 < nx
            i0_east = i0 + 1;
        else
            i0_east = 1;
        end
        
        % eastern equations        
        if isfinite(gamma_initial(k_east(k0,j0,i0),j0,i0_east)) & ...
                isfinite(gamma_initial(k_east(k0,j0,i0)+1,j0,i0_east))
            neq = neq + 1;
            ieq = ieq + 1;
            ks(neq) = k0;
            js(neq) = j0;
            is(neq) = i0;
            s1(ieq) = neq;
            s2(ieq) = kg;
            s3(ieq) = 1*wt;
            
            ind_gu = ginds(k_east(k0,j0,i0),j0,i0_east);
            c_gu = -(1-r_east(k0,j0,i0));
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ind_gu;
            s3(ieq) = c_gu*wt;
            
            ind_gl = ginds(k_east(k0,j0,i0)+1,j0,i0_east);
            c_gl = -r_east(k0,j0,i0);
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ind_gl;
            s3(ieq) = c_gl*wt;
        end
    end
    
    %  northern equations
    if isfinite(r_north(I_gamma(kg)))  
        if j0 < ny
            j0_north = j0 + 1;
        else
            error('this north not possible!')
        end
        
        if j0 < ny & ...
            isfinite(gamma_initial(k_north(k0,j0,i0),j0_north,i0)) & ...
              isfinite(gamma_initial(k_north(k0,j0,i0)+1,j0_north,i0))
            neq = neq + 1;
            ieq = ieq + 1;
            ks(neq) = k0;
            js(neq) = j0;
            is(neq) = i0;
            s1(ieq) = neq;
            s2(ieq) = kg;
            s3(ieq) = 1*wt;
            
            ind_gu = ginds(k_north(k0,j0,i0),j0_north,i0);
            c_gu = -(1 - r_north(k0,j0,i0));
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ind_gu;
            s3(ieq) = c_gu*wt;
            
            ind_gl = ginds(k_north(k0,j0,i0)+1,j0_north,i0);
            c_gl = -r_north(k0,j0,i0);
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ind_gl;
            s3(ieq) = c_gl*wt;
        end
    end
    
    %  western equations
    if isfinite(r_west(I_gamma(kg)))
        if i0 > 1
            i0_west = i0-1;
        else
            i0_west = nx;
        end
        
        if isfinite(gamma_initial(k_west(k0,j0,i0),j0,i0_west)) & ...
                isfinite(gamma_initial(k_west(k0,j0,i0)+1,j0,i0_west))
            neq = neq+1;
            ieq = ieq+1;
            ks(neq) = k0;
            js(neq) = j0;
            is(neq) = i0;
            s1(ieq) = neq;
            s2(ieq) = kg;
            s3(ieq) = 1*wt;
            
            ind_gu = ginds(k_west(k0,j0,i0),j0,i0_west);
            c_gu = -(1 - r_west(k0,j0,i0));
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ind_gu;
            s3(ieq) = c_gu*wt;
            
            ind_gl = ginds(k_west(k0,j0,i0)+1,j0,i0_west);
            c_gl = -r_west(k0,j0,i0);
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ind_gl;
            s3(ieq) = c_gl*wt;
            
        end
    end
    
    % southern equations
    if isfinite(r_south(I_gamma(kg)))
        if j0 > 1
            j0_south = j0 - 1;
        else
            error('this south not possible!')
        end
        
        if j0 > 1 & ...
           isfinite(gamma_initial(k_south(k0,j0,i0),j0_south,i0)) & ...
           isfinite(gamma_initial(k_south(k0,j0,i0)+1,j0_south,i0))
       
            neq = neq + 1;
            ieq = ieq + 1;
            ks(neq) = k0;
            js(neq) = j0;
            is(neq) = i0;
            s1(ieq) = neq;
            s2(ieq) = kg;
            s3(ieq) = 1*wt;
            
            ind_gu = ginds(k_south(k0,j0,i0),j0_south,i0);
            c_gu = -(1 - r_south(k0,j0,i0));
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ind_gu;
            s3(ieq) = c_gu*wt;
            
            ind_gl = ginds(k_south(k0,j0,i0)+1,j0_south,i0);
            c_gl = -r_south(k0,j0,i0);
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ind_gl;
            s3(ieq) = c_gl*wt;
        end
    end
end

ieq_h = ieq;
neq_h = neq;

%dj_toc

%  vertical b equations
if vertical_eqns == 1
    %dj_disp('setting up vertical equations ... ')
    for kg = 1:ng
        [k0,j0,i0] = ind2sub([nz,ny,nx],I_gamma(kg));
        wt = wts(kg);
        kvert = kvert + 1;
        
        if kvert > no_v_eqns | k0 ~= k_vert(kvert)
            error('*** vertical equations & data mis-match ***')
        end
        rr = r_vert(kvert);
        gprod = gprod_vert(kvert);
        if k0 == 1
            % dgamma = g(k0+1,j0,i0)-g(k0,j0,i0); 
            % dp = p(k0+1,j0,i0)-p(k0,j0,i0);
            ind_gu = ginds(k0+1,j0,i0);
            ind_gl = ginds(k0+2,j0,i0);
        elseif k0 == n(j0,i0)
            % dgamma = g(k0,j0,i0)-g(k0-1,j0,i0); 
            % dp = p(k0,j0,i0)-p(k0-1,j0,i0);
            ind_gu = ginds(k0-2,j0,i0);
            ind_gl = ginds(k0-1,j0,i0);
        else
            % dgamma = (g(k0+1,j0,i0)-g(k0-1,j0,i0))/2; 
            % dp = (p(k0+1,j0,i0)-p(k0-1,j0,i0))/2;
            ind_gu = ginds(k0-1,j0,i0);
            ind_gl = ginds(k0+1,j0,i0);
        end
        
        neq = neq + 1;
        ieq = ieq + 1;
        if neq > no_eqns
            disp('neq > no_eqns')
        end
        if ieq > ns,
            disp('ieq > ns')
        end
        ks(neq) = k0; 
        js(neq) = j0; 
        is(neq) = i0;
        s1(ieq) = neq; 
        s2(ieq) = kg; 
        s3(ieq) = wt*gprod;
        
        c_gu = -rr;
        ieq = ieq + 1;
        s1(ieq) = neq; 
        s2(ieq) = ind_gu; 
        s3(ieq) = wt*c_gu;
        
        c_gl = -(gprod - rr);
        ieq = ieq + 1;
        s1(ieq) = neq;
        s2(ieq) = ind_gl;
        s3(ieq) = wt*c_gl;
        
        b(neq) = 0; 
    end
end

ieq_v = ieq; 
neq_v = neq;

%dj_toc
% the z-x b equations
if bxz_eqns == 1
    for kg = 1:ng
        [k0,j0,i0] = ind2sub([nz,ny,nx],I_gamma(kg));
        wt = p(k0,j0,i0)^kwt;
        if isnan(wt)
            wt = 1; 
        end
        
        if isfinite(gc1(kg))
            kxzbeq = kxzbeq + 1;
            if kxzbeq > nxz_beqs
                error('*** x-z b equations data mis-match ***')
            end
            if k0 == 1
                dj_disp('error dj1')
            elseif k0 == n(j0,i0)
                error_dj2 = [kg,I_gamma(kg),k0,j0,i0,n(j0,i0),gc1(kg)]
                error('   *****')
            else
                ind_gu = ginds(k0-1,j0,i0);
                ind_gl = ginds(k0+1,j0,i0);
            end
            
            if i0 == 1
                dj_disp('error dj3')
            elseif k0 == n(j0,i0)
                dj_disp('error dj4')
            else
                ind_ge = ginds(k0,j0,i0+1);
                ind_gw = ginds(k0,j0,i0-1);
            end
            
            neq = neq + 1;
            ieq = ieq + 1;
            
            ks(neq) = k0;
            js(neq) = j0;
            is(neq) = i0;
            
            s1(ieq) = neq;
            s2(ieq) = kg;
            s3(ieq) = wt*gc1(kg);
            
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ind_gu;
            s3(ieq) = wt*gc2(kg);
            
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ind_gl;
            s3(ieq) = wt*gc3(kg);
            
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ind_ge;
            s3(ieq) = wt*gc4(kg);
            
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ind_gw;
            s3(ieq) = wt*gc5(kg);
            
            b(neq) = wt*0;
            
        end
        
        if mod(kg,floor(0.2*ng)) == 0
            thus_far_inversion = [kg, round(100*kg/ng)] 
            %dj_toc
        end
    end
end

ieq_vxz = ieq; 
neq_vxz = neq;

% gamma_n equations
wt = handles.bdrywt;

if gn_eqns == 1
    if no_b_eqns<20
        inds_bg_w = [1:no_b_eqns];
    else
        n10 = floor(no_b_eqns/10); 
        inds_bg_w = [1:n10:no_b_eqns];
    end
    
    %[I_bg(inds_bg_w(:)),gamma_initial(I_bg(inds_bg_w(:))),gamma_bdry(inds_bg_w(:))]
    
    for kk = 1:no_b_eqns
        if I_bg(kk) > 0
            [k0,j0,i0] = ind2sub([nz,ny,nx],I_bg(kk));
            if gamma_bdry(kk) > 0
                neq = neq + 1; 
                ieq = ieq + 1;
                ks(neq) = k0; 
                js(neq) = j0; 
                is(neq) = i0;
                kg = ginds(k0,j0,i0);
                s1(ieq) = neq; 
                s2(ieq) = kg; 
                s3(ieq) = wt;
                b(neq) = wt*gamma_bdry(kk);
                %[neq, s2(ieq), b(neq)]
                if ~isfinite(neq) | ~isfinite(kg)
                    an_error = [kk,neq,kg,gamma_bdry(kk)]
                end
            else
                error('boundary gamma value not possible')
            end
        end
    end  
end

ieq_g = ieq; 
neq_g = neq;

% h_eq_stats = [ieq_h,neq_h]
% v_eq_stats = [ieq_v,neq_v]
% vxz_eq_stats = [ieq_vxz,neq_vxz]
% g_eq_stats = [ieq_g,neq_g]

s1 = s1(1:ieq); 
s2 = s2(1:ieq); 
s3 = s3(1:ieq); 
b = b(1:neq);

ks = ks(1:neq); 
js = js(1:neq); 
is = is(1:neq);

neqs = neq;

actual_no_equations = neqs

% final weights
horizontal_max = max([s1(1:ieq_h),s2(1:ieq_h),s3(1:ieq_h)])
horizontal_b_max = max(b(1:neq_h))

wt_h = handles.h2vwt/max(s3(1:ieq_h));
s3(1:ieq_h) = s3(1:ieq_h)*wt_h;
b(1:neq_h) = b(1:neq_h)*wt_h;

vertical_z_max = max([s1(ieq_h+1:ieq_v),s2(ieq_h+1:ieq_v),s3(ieq_h+1:ieq_v)])
vertical_z_b_max = max(b(neq_h+1:neq_v))

if max(s3(ieq_h+1:ieq_v)) == 0
    error('dj5_error')
end

wt_v = 1/max(s3(ieq_h+1:ieq_v));

s3(ieq_h+1:ieq_v) = s3(ieq_h+1:ieq_v)*wt_v;
b(neq_h+1:neq_v) = b(neq_h+1:neq_v)*wt_v;

vertical_xz_max = max([s1(ieq_v+1:ieq_vxz),s2(ieq_v+1:ieq_vxz),s3(ieq_v+1:ieq_vxz)])
vertical_xz_b_max = max(b(neq_v+1:neq_vxz))

if max(s3(ieq_v+1:ieq_vxz)) == 0
    error('dj6_error')
end

wt_vxz = 1/max(s3(ieq_v+1:ieq_vxz));

s3(ieq_v+1:ieq_vxz) = s3(ieq_v+1:ieq_vxz)*wt_vxz;
b(neq_v+1:neq_vxz) = b(neq_v+1:neq_vxz)*wt_vxz;

no_isolated_gamma = length(setdiff((1:ng),unique(s2(:))))

%dj_pause(0)

A = sparse(s1,s2,s3,neqs,ng);

clear s1 s2 s3
%dj_toc
ok = ['inverting ', method, 'ly ...']

% if strcmp(method,'direct')    
%     g1 = A\b;
%     
%     ok = 'done',    
%     gamma_i = nans(size(s));
%     gamma_i(inds_g) = g1;
%     
%     adjust_gammas
%     
%     stats = [min(g(inds_g)), mean(g(inds_g)), max(g(inds_g))]
%     
%     save  out s t ct p  gamma_i longs lats ocean n
% 
%     resids = b - A*g1(:);
%     size(resids)
%     size(g1)
    %    figure, plot(g1,resids,'.'), grid on, dj_pause(5)
    
%     figure
%     hist(resids,500)
%     grid on
%     dj_pause(3)
%     
%     figure
%     plot(resids,ks,'.')
%     grid on
%     set(gca,'ydir','reverse')
%     dj_pause(3)
    
%     gg = nan(neq,1);
%     for keq = 1:neq
%         if gamma_i(ks(keq),js(keq),is(keq))>1
%             gg(keq) = gamma_i(ks(keq),js(keq),is(keq));
%         end
%     end
%     
%     save
    
%     figure
%     plot(resids,gg,'.')
%     grid on
%     set(gca,'ydir','reverse')
    
% figure(gcf)
% dj_pause(0) 
% dj_toc
% cleanup
    
    %        Veronis error
%     cd ../neutral_density
%     
%     figure
% % dj_pause(1)
%     [ave0,perc0] = get_gradients(handles);
%     
%     global logD
%     max_logD = nanmax(logD(:))
    
%     cd ../inversion
    %  iterative method
    
% elseif strcmp(method,'iterative')

    tol = 1e-5;
    maxits = handles.maxits;
%     itlo = 0;
%     iters = 0;
    
    for kk = 1:10000    
       % iter = itlo + kk*maxits;
        g1 = lsqr(A,b,tol,maxits,[],[],gamma_old);
        g1_in_concrete = g1;
        gamma_i = nan(size(s));
        gamma_i(I_gamma) = g1;
        
%         %        adjust_gammas
%         %        stats = [iter,min(g(inds_g)), mean(g(inds_g)), max(g(inds_g))]
%         %        Veronis error
% %         cd ../neutral_density
%         [ave0,perc0] = get_gradients(handles);
% %         cd ../inversion
%         %         set(gcf,'Name','current profile','NumberTitle','off')
%         ave = [ave; ave0];
%         perc = [perc; perc0];
%         iters = [iters; iter];
%         
%         gamma_mean = [gamma_mean; gmean];
%         gamma_std = [gamma_std; gstd];
        
%         if mod(kk,5) == 0
%             figure(3)
%             subplot(2,2,1)
%             plot(iters,ave)
%             grid on
%             title(['\mu(D_v)   =   ', num2str(ave0,3)])
%             
%             subplot(2,2,2)
%             plot(iters,perc)
%             grid on, 
%             title(['v_5%   =   ', num2str(perc0,3)])
%             
%             subplot(2,2,3)
%             plot(iters,gamma_mean)
%             grid on
%             title(['\mu(\gamma)   =   ', num2str(gmean,3)])
%             
%             subplot(2,2,4)
%             plot(iters,gamma_std)
%             grid on
%             title(['\sigma(\gamma)   =   ', num2str(gstd,3)])
%             set(gcf,'Position',[708   406   569   400],'Name','iterations','NumberTitle','off')
%             
%             if no_b_eqns<20
%                 inds_bg_w = [1:no_b_eqns];
%             else
%                 n10 = floor(no_b_eqns/10);
%                 inds_bg_w = [1:n10:no_b_eqns];
%             end
%             boundary_gammas = [inds_bg(inds_bg_w(:)),g(inds_bg(inds_bg_w(:))),g_bdry(inds_bg_w(:))]
%             dj_pause(1)
%         end      
        gamma_old = g1_in_concrete;    
%         dj_toc
    end 
%end


end
