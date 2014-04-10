function g_prf = gamma_prf_7_9(s,ct,g,ocean,lats,c)


c1 = [c(1:7); c(15:23)]; c2 = [c(8:14); c(15:23)];

the_cs = [c1, c2]

[nz,ny,nx] = size(g);

wt1 = nan*ones(size(ocean)); wt2 = wt1;

inds_g = find(finite(g(1,:))); wt1(inds_g) = 0; wt2(inds_g) = 0;

grf = nan*ones(size(g));

lat_lo = 32; lat_hi = 40;

for j = 1:ny
    for i = 1:nx
        if ocean(j,i)==5 & finite(g(1,j,i))
            
            if lats(j,i)>=lat_hi
                wt1(j,i) = 1;
            elseif lats(j,i)<=lat_lo
                wt2(j,i) = 1;
            else
                wt1(j,i) = (lats(j,i)-lat_lo)/(lat_hi-lat_lo);  wt2(j,i) = 1-wt1(j,i);
            end  
            
            c = wt1(j,i)*c1 + wt2(j,i)*c2;
            
            indsz = find(finite(g(:,j,i)));
            
            gg_prf = rfunc_7_9(s(indsz,j,i),ct(indsz,j,i),c);
            
            g_prf(indsz,j,i) = gg_prf;
            
        end
    end
end


return