function [k_vert,r_vert,gprod_vert,no_eqs] = vertical_vitals(s,t,p,g,ocean,n,longs,lats,handles)

cd intersections

load intersections_east
load intersections_west
load intersections_north
load intersections_south

cd .. 

[nz,ny,nx] = size(s);

inds_g = find(finite(g)); ng = length(inds_g)

%  inds_e = find(finite(r_east)); neast = length(inds_e); perc_east = 100*neast/ng
%  inds_n = find(finite(r_north)); nnth = length(inds_n); perc_nth = 100*nnth/ng
%  inds_w = find(finite(r_west)); nwest = length(inds_w); perc_west = 100*nwest/ng
%  inds_s = find(finite(r_south)); nsth = length(inds_s); perc_sth = 100*nsth/ng

nr = ones(ng,1); k_vert = nan*ones(ng,1); r_vert = k_vert; gprod_vert = k_vert;

for kg = 1:ng
    
    if finite(r_east(inds_g(kg)))
      nr(kg) = nr(kg)+1;
    end
    
    if finite(r_west(inds_g(kg)))
      nr(kg) = nr(kg)+1;
    end
    
    if finite(r_north(inds_g(kg)))
      nr(kg) = nr(kg)+1;
    end
    
    if finite(r_south(inds_g(kg)))
      nr(kg) = nr(kg)+1;
    end
    
    if nr(kg)>=1
        [k,j,i] = ind2sub([nz,ny,nx],inds_g(kg)); indss = find(finite(s(:,j,i)));
        ss = s(indss,j,i); ctt = t(indss,j,i); pp = p(indss,j,i);
        k_vert(kg) = k;
        if k==1
            pmid = (pp(k)+pp(k+1))/2;
            dsig_u = rho_from_ct(ss(k+1),ctt(k+1),pmid)-rho_from_ct(ss(k),ctt(k),pmid); dsig_u2 = dsig_u*dsig_u;
            pmid = (pp(k)+pp(k+2))/2;
            dsig_m = rho_from_ct(ss(k+2),ctt(k+2),pmid)-rho_from_ct(ss(k),ctt(k),pmid); dsig_m2 = dsig_m*dsig_m;
            pmid = (pp(k+1)+pp(k+2))/2;
            dsig_l = rho_from_ct(ss(k+2),ctt(k+2),pmid)-rho_from_ct(ss(k+1),ctt(k+1),pmid); dsig_l2 = dsig_l*dsig_l;
            gprod =  dsig_l*(dsig_m2+dsig_u2); 
            r = (dsig_u+dsig_l)*dsig_m2+dsig_m*dsig_u2;
        elseif k==n(j,i)
            pmid = (pp(k-2)+pp(k-1))/2;
            dsig_u = rho_from_ct(ss(k-1),ctt(k-1),pmid)-rho_from_ct(ss(k-2),ctt(k-2),pmid); dsig_u2 = dsig_u*dsig_u;
            pmid = (pp(k-2)+pp(k))/2;
            dsig_m = rho_from_ct(ss(k),ctt(k),pmid)-rho_from_ct(ss(k-2),ctt(k-2),pmid); dsig_m2 = dsig_m*dsig_m;
            pmid = (pp(k-1)+pp(k))/2;
            dsig_l = rho_from_ct(ss(k),ctt(k),pmid)-rho_from_ct(ss(k-1),ctt(k-1),pmid); dsig_l2 = dsig_l*dsig_l;
            gprod =  dsig_u*(dsig_l2+dsig_m2); 
            r = (dsig_u-dsig_m)*dsig_l2-dsig_l*dsig_m2;
        else
            pmid = (pp(k-1)+pp(k))/2; pu = pmid;
            dsig_u = rho_from_ct(ss(k),ctt(k),pmid)-rho_from_ct(ss(k-1),ctt(k-1),pmid); dsig_u2 = dsig_u*dsig_u;
            pmid = (pp(k-1)+pp(k+1))/2; pm = pmid;
            dsig_m = rho_from_ct(ss(k+1),ctt(k+1),pmid)-rho_from_ct(ss(k-1),ctt(k-1),pmid); dsig_m2 = dsig_m*dsig_m;
            pmid = (pp(k)+pp(k+1))/2; pl = pmid;
            dsig_l = rho_from_ct(ss(k+1),ctt(k+1),pmid)-rho_from_ct(ss(k),ctt(k),pmid); dsig_l2 = dsig_l*dsig_l;
%            zz = handles.linb
            if handles.linb==1
              gprod =  dsig_m*(dsig_l*(pl-pm)+dsig_u*(pu-pm)); 
              r = dsig_l*dsig_m*(pl-pm)+dsig_l*dsig_u*(pu-pl);
            else
              gprod =  dsig_m*(dsig_l2+dsig_u2); 
              r = (dsig_m-dsig_u)*dsig_l2+dsig_l*dsig_u2;
            end
        end
        nr(kg) = nr(kg)+1;
        zz = [gprod, r, gprod-r];
        [cmax,kmax] = max(abs(zz));                   
        cmax = sign(zz(kmax))*cmax;
        if abs(cmax)>0
            gprod = gprod/cmax; r = r/cmax;
        else
            if k==1
                gprod = 1; r = 1;
            elseif k==n(j,i)
                gprod = 1; r = 0;
            else
                gprod = 1; r = 0.5;
            end
        end
        r_vert(kg) = r; gprod_vert(kg) = gprod;
    end
    
    if mod(kg,floor(ng/5))==0
        done = 100*kg/ng
        dj_toc
    end
            
end


perc_r = 100*length(find(nr>1))/ng;

no_v_eqs = [ng-ng*perc_r/100, length(k_vert), length(r_vert), length(gprod_vert)]

z = nr-ones(ng,1); no_eqs = sum(z(:))

figeuure(3), plot(r_vert,'.'), grid on
 
set(gcf,'Name','vertical vitals','NumberTitle','off','Color',[0.961 0.988 0.965])

return
