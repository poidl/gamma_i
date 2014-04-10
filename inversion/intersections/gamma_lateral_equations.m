function [] = gamma_lateral_equations(SA,gamma_initial,N2,r_east,r_west,r_north,r_south,handles)

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


end