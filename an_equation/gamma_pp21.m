function g_pp = gamma_pp21(s,ct,g,ocean,lats,c)

npoly = 21;

c1 = c(1:npoly); 
c2 = c(npoly+1:2*npoly); 
c3 = c(2*npoly+1:3*npoly);
c4 = c(3*npoly+1:4*npoly);

load wts1234

[nz,ny,nx] = size(g);

g_pp = nan*ones(size(g));

for j = 1:ny
    for i = 1:nx
        
        if finite(g(1,j,i))
            
            c = wt1(j,i)*c1 + wt2(j,i)*c2 + wt3(j,i)*c3 + wt4(j,i)*c4;
            
            indsg = find(finite(g(:,j,i)));
            
            gg_pp = gamma_p21(s(indsg,j,i),ct(indsg,j,i),c);
            
            g_pp(indsg,j,i) = gg_pp;
            
        end
        
    end
end


return