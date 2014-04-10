function g_pp = gamma_pp15(s,ct,g,ocean,lats,c)


c1 = c(1:15); c2 = c(16:30); c3 = c(31:45); c4 = c(46:60);

load wts1234


[nz,ny,nx] = size(g);

g_pp = nan*ones(size(g));


for j = 1:ny
    for i = 1:nx
        
        if finite(g(1,j,i))
            
            c = wt1(j,i)*c1 + wt2(j,i)*c2 + wt3(j,i)*c3 + wt4(j,i)*c4;
            
            indsg = find(finite(g(:,j,i)));
            
            gg_pp = gamma_p15(s(indsg,j,i),ct(indsg,j,i),c);
            
            g_pp(indsg,j,i) = gg_pp;
            
        end
        
    end
end


return