function [k_vert,r_vert,gprod_vert,no_eqs] = vertical_vitals(SA,CT,p,g,n,r_north,r_east,r_south,r_west,handles)

% cd intersections
% load intersections_east
% load intersections_west
% load intersections_north
% load intersections_south
% cd .. 

[nz,ny,nx] = size(SA);
I_g_data = find(isfinite(g));
ng = length(I_g_data);
%  inds_e = find(finite(r_east)); neast = length(inds_e); perc_east = 100*neast/ng
%  inds_n = find(finite(r_north)); nnth = length(inds_n); perc_nth = 100*nnth/ng
%  inds_w = find(finite(r_west)); nwest = length(inds_w); perc_west = 100*nwest/ng
%  inds_s = find(finite(r_south)); nsth = length(inds_s); perc_sth = 100*nsth/ng
nr = ones(ng,1);
k_vert = nans(ng,1);
r_vert = k_vert; 
gprod_vert = k_vert;

for kg = 1:ng
    if isfinite(r_east(I_g_data(kg)))
      nr(kg) = nr(kg) + 1;
    end
    
    if isfinite(r_west(I_g_data(kg)))
      nr(kg) = nr(kg) + 1;
    end
    
    if isfinite(r_north(I_g_data(kg)))
      nr(kg) = nr(kg) + 1;
    end
    
    if isfinite(r_south(I_g_data(kg)))
      nr(kg) = nr(kg) + 1;
    end
    
    if nr(kg) >= 1
        [k, j, i] = ind2sub([nz,ny,nx], I_g_data(kg));
        I_cast_data = find(isfinite(SA(:,j,i)));
        SA_cast = SA(I_cast_data,j,i);
        CT_cast = CT(I_cast_data,j,i);
        p_cast = p(I_cast_data,j,i);
        k_vert(kg) = I_g_data(k);
        if k == 1
            pmid = 0.5*(p_cast(k) + p_cast(k+1));
            dsig_u = gsw_rho(SA_cast(k+1),CT_cast(k+1),pmid) - gsw_rho(SA_cast(k),CT_cast(k),pmid);
            dsig_u2 = dsig_u*dsig_u;
            
            pmid = 0.5*(p_cast(k) + p_cast(k+2));
            dsig_m = gsw_rho(SA_cast(k+2),CT_cast(k+2),pmid) - gsw_rho(SA_cast(k),CT_cast(k),pmid);
            dsig_m2 = dsig_m*dsig_m;
            
            pmid = 0.5*(p_cast(k+1) + p_cast(k+2));
            dsig_l = gsw_rho(SA_cast(k+2),CT_cast(k+2),pmid) - gsw_rho(SA_cast(k+1),CT_cast(k+1),pmid);
            %dsig_l2 = dsig_l*dsig_l;
            
            gprod =  dsig_l*(dsig_m2 + dsig_u2); 
            r = (dsig_u + dsig_l)*dsig_m2 + dsig_m*dsig_u2;
        elseif k == n(j,i)
            pmid = 0.5*(p_cast(k-2) + p_cast(k-1));
            dsig_u = gsw_rho(SA_cast(k-1),CT_cast(k-1),pmid) - gsw_rho(SA_cast(k-2),CT_cast(k-2),pmid);
            %dsig_u2 = dsig_u*dsig_u;
            
            pmid = 0.5*(p_cast(k-2) + p_cast(k));
            dsig_m = gsw_rho(SA_cast(k),CT_cast(k),pmid) - gsw_rho(SA_cast(k-2),CT_cast(k-2),pmid);
            dsig_m2 = dsig_m*dsig_m;
            
            pmid = 0.5*(p_cast(k-1) + p_cast(k));
            dsig_l = gsw_rho(SA_cast(k),CT_cast(k),pmid) - gsw_rho(SA_cast(k-1),CT_cast(k-1),pmid);
            dsig_l2 = dsig_l*dsig_l;
            
            gprod =  dsig_u*(dsig_l2 + dsig_m2); 
            r = (dsig_u - dsig_m)*dsig_l2 - dsig_l*dsig_m2;
        else
            %[k,n(j,i),longs(j,i),lats(j,i)];
            pu = 0.5*(p_cast(k-1) + p_cast(k));
            dsig_u = gsw_rho(SA_cast(k),CT_cast(k),pu) - gsw_rho(SA_cast(k-1),CT_cast(k-1),pu);
            dsig_u2 = dsig_u*dsig_u;
            
            pmid = 0.5*(p_cast(k-1) + p_cast(k+1));
            dsig_m = gsw_rho(SA_cast(k+1),CT_cast(k+1),pmid) - gsw_rho(SA_cast(k-1),CT_cast(k-1),pmid);
            %dsig_m2 = dsig_m*dsig_m;
            
            pl = 0.5*(p_cast(k) + p_cast(k+1));
            dsig_l = gsw_rho(SA_cast(k+1),CT_cast(k+1),pl) - gsw_rho(SA_cast(k),CT_cast(k),pl);
            dsig_l2 = dsig_l*dsig_l;
%            zz = handles.linb
            if handles.linb == 1
              gprod =  dsig_m*(dsig_l*(pl - pmid) + dsig_u*(pu - pmid)); 
              r = dsig_l*dsig_m*(pl - pmid) + dsig_l*dsig_u*(pu - pl);
            else
              gprod =  dsig_m*(dsig_l2 + dsig_u2); 
              r = (dsig_m - dsig_u)*dsig_l2 + dsig_l*dsig_u2;
            end
        end
        nr(kg) = nr(kg) + 1;
        zz = [gprod, r, gprod-r];
        [cmax,kmax] = max(abs(zz));                   
        cmax = sign(zz(kmax))*cmax;
        if abs(cmax) > 0
            gprod = gprod/cmax;
            r = r/cmax;
        else
            if k == 1
                gprod = 1;
                r = 1;
            elseif k == n(j,i)
                gprod = 1;
                r = 0;
            else
                gprod = 1;
                r = 0.5;
            end
        end
        r_vert(kg) = r;
        gprod_vert(kg) = gprod;
    end
    
%     if mod(kg,floor(ng/5)) == 0
%         done = 100*kg/ng;
%         %dj_toc
%     end
            
end

%perc_r = 100*length(find(nr>1))/ng;

%no_v_eqs = [ng - ng*perc_r/100, length(k_vert), length(r_vert), length(gprod_vert)];

z = nr - ones(ng,1);
no_eqs = sum(z(:));

% figure(5)
% plot(r_vert,'.') 
% grid on
% set(gcf,'Name','vertical vitals','NumberTitle','off','Color',[0.961 0.988 0.965])

end
