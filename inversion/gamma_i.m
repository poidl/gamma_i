function g = gamma_i(SA,CT,p,g,lat,long,ocean,n)

%global s t ct p z g longs lats ocean n
% global called_ggrads gmean gstd

% clc

[nz,ny,nx] = size(s);

% called_ggrads = 0;

handles.inv_meth = 2;
% if handles.inv_meth == 1
%     method = 'direct'
% else
    method = 'iterative'
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

% initial error profile
% cd ../neutral_density
% [ave, perc] = get_gradients(handles);
% gamma_mean = gmean;
% gamma_std = gstd;
% handles.ntp = 0;
% set(gcf,'Position',[4   405   606   401])

% load equation set-up information

%[k_east,r_east] = intersections_east(SA,CT,p,ocean,n,long,lat,)
load gk_interp_int_east

%[k_west,r_west] = intersections_west(SA,CT,p,ocean,n,long,lat)
load gk_interp_int_west

%[k_north,r_north] = intersections_north(SA,CT,p,ocean,n,long,lat)
load gk_interp_int_north

%[k_south,r_south] = intersections_south(SA,CT,p,ocean,n,long,lat)
load gk_interp_int_south

%[inds_bg, g_bdry] = boundary_gammas(g,long,lat,n);
load gk_interp_gamma_boundary    

% load helicities

%[k_vert,r_vert,gprod_vert,no_eqs] = vertical_vitals(SA,CT,p,g,n,r_north,r_east,r_south,r_west)
load gk_interp_vert_vitals

% load bxz_equations

% determine # equations and # unknowns
inds = find(isfinite(r_east));
no_h_eqns = length(inds);

inds = find(isfinite(r_west));
no_h_eqns = no_h_eqns + length(inds);

inds = find(isfinite(r_north));
no_h_eqns = no_h_eqns + length(inds);

inds = find(isfinite(r_south));
no_h_eqns = no_h_eqns + length(inds);

no_v_eqns = 0;
if vertical_eqns == 1
    no_v_eqns = no_eqs;
end

no_b_eqns = length(inds_bg);

no_eqns = no_h_eqns + no_v_eqns + no_b_eqns;

inds_g = find(isfinite(g));
ng = length(inds_g);

no_equations = no_eqns;
no_unknowns = ng;

% determine N2 for weighting
% dj_disp('computing N2 ... ')
% dj_tic
inds = find(isfinite(s(1,:)));
nn = length(inds);
ss = s(:,inds);
tt = t(:,inds);
pp = p(:,inds);
N2 = nan(size(s));
for kk = 1:nn
    Iz = find(isfinite(ss(:,kk)));
    %N2(Iz,inds(kk)) = bfrq1(ss(Iz,kk),tt(Iz,kk),pp(Iz,kk),'ct');
    [dummyN2,dummyp_mid] = gsw_Nsquared(ss(Iz,kk),tt(Iz,kk),pp(Iz,kk));
    N2(Iz,inds(kk)) = gsw_interp_Nsquared(dummyN2,dummyp_mid,pp(Iz,kk));
end
% dj_toc

ginds = nan(size(s));
ginds(inds_g) = [1:ng];
g0 = g(inds_g);

g0_min = min(g0);
g0_max = max(g0);

kvert = 0;
kxzbeq = 0;

ns = 3*(no_h_eqns + no_v_eqns) + no_b_eqns;

neq = 0;
ieq = 0;
s1 = nan(ns,1);
s2 = s1;
s3 = s1;
b = zeros(no_eqns,1);

ks = nan*ones(no_eqns,1);
js = ks;
is = ks;
wts = ks;

% dj_disp('setting up lateral equations ... ')
for kg = 1:ng
    [k0,j0,i0] = ind2sub([nz,ny,nx],inds_g(kg));
    helwt = handles.helwt;
    kwt = 1;
    % wt = (1+p(k0,j0,i0)/kwt)^helwt;
    wt = 1/(helwt+N2(k0,j0,i0));
    wts(kg) = wt;
    ks(kg) = k0;
    js(kg) = j0;
    is(kg) = i0;
      
    % lateral equations
    if isfinite(r_east(inds_g(kg)))        
        if i0 < nx
            i0_east = i0 + 1;
        else
            i0_east = 1;
        end
        
        % eastern equations        
        if isfinite(g(k_east(k0,j0,i0),j0,i0_east)) & isfinite(g(k_east(k0,j0,i0)+1,j0,i0_east))
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
    if isfinite(r_north(inds_g(kg)))  
        if j0 < ny
            j0_north = j0 + 1;
        else
            error('this north not possible!')
        end
        
        if j0 < ny & ...
            isfinite(g(k_north(k0,j0,i0),j0_north,i0)) & ...
              isfinite(g(k_north(k0,j0,i0)+1,j0_north,i0))
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
    if isfinite(r_west(inds_g(kg)))
        if i0 > 1
            i0_west = i0-1;
        else
            i0_west = nx;
        end
        
        if isfinite(g(k_west(k0,j0,i0),j0,i0_west)) & isfinite(g(k_west(k0,j0,i0)+1,j0,i0_west))
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
    if isfinite(r_south(inds_g(kg)))
        if j0 > 1
            j0_south = j0 - 1;
        else
            error('this south not possible!')
        end
        
        if j0 > 1 & ...
           isfinite(g(k_south(k0,j0,i0),j0_south,i0)) & ...
           isfinite(g(k_south(k0,j0,i0)+1,j0_south,i0))
       
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
        [k0,j0,i0] = ind2sub([nz,ny,nx],inds_g(kg));
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
        [k0,j0,i0] = ind2sub([nz,ny,nx],inds_g(kg));
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
                error_dj2 = [kg,inds_g(kg),k0,j0,i0,n(j0,i0),gc1(kg)]
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
    
    [inds_bg(inds_bg_w(:)),g(inds_bg(inds_bg_w(:))),g_bdry(inds_bg_w(:))]
    
    for kk = 1:no_b_eqns
        if inds_bg(kk) > 0
            [k0,j0,i0] = ind2sub([nz,ny,nx],inds_bg(kk));
            if g_bdry(kk) > 0
                neq = neq + 1; 
                ieq = ieq + 1;
                ks(neq) = k0; 
                js(neq) = j0; 
                is(neq) = i0;
                kg = ginds(k0,j0,i0);
                s1(ieq) = neq; 
                s2(ieq) = kg; 
                s3(ieq) = wt;
                b(neq) = wt*g_bdry(kk);
                %[neq, s2(ieq), b(neq)]
                if ~isfinite(neq) | ~isfinite(kg)
                    an_error = [kk,neq,kg,g_bdry(kk)]
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

%         final weights
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
%     g = nans(size(s));
%     g(inds_g) = g1;
%     
%     adjust_gammas
%     
%     stats = [min(g(inds_g)), mean(g(inds_g)), max(g(inds_g))]
%     
%     save  out s t ct p  g longs lats ocean n
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
%         if g(ks(keq),js(keq),is(keq))>1
%             gg(keq) = g(ks(keq),js(keq),is(keq));
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
        g1 = lsqr(A,b,tol,maxits,[],[],g0);
        g1_in_concrete = g1;
        g = nan(size(s));
        g(inds_g) = g1;
        
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
        g0 = g1_in_concrete;    
%         dj_toc
    end 
%end


end
