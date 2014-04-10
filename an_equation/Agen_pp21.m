function [A,b] = Agen_pp21(s,ct,p,g,ocean,longs,lats)


indsg = find(finite(g)); ng = length(indsg);

dj_disp(['number of equations = ', int2str(ng)])

A = nan*ones(ng,84); b = nan*ones(ng,1);

[nz,ny,nx] = size(g);

load wts1234

neq = 0; initialise_waitbar

for j = 1:ny
for i = 1:nx
        
	if finite(g(1,j,i))
                     
        indsg = find(finite(g(:,j,i)));
        
        nn = length(indsg); ss = s(indsg,j,i); ctt = ct(indsg,j,i); gg = g(indsg,j,i);
        
        for k = 1:4
            
            cmd = ['wt = wt', int2str(k), ';']; eval(cmd)
            
            kk = (k-1)*21;
        
            A(neq+(1:nn), kk+1)  = wt(j,i)*ones(nn,1);   
            A(neq+(1:nn), kk+2)  = wt(j,i)*ss;
            A(neq+(1:nn), kk+3)  = wt(j,i)*ctt;
    
            A(neq+(1:nn), kk+4)  = wt(j,i)*ss.*ss;                               %   s^2
            A(neq+(1:nn), kk+5)  = wt(j,i)*ss.*ctt;                              %   s   t
            A(neq+(1:nn), kk+6)  = wt(j,i)*ctt.*ctt;                             %       t^2
    
            A(neq+(1:nn), kk+7)  = wt(j,i)*ss.*A(neq+(1:nn), kk+4);              %   s^3
            A(neq+(1:nn), kk+8)  = wt(j,i)*ss.*A(neq+(1:nn), kk+5);              %   s^2 t
            A(neq+(1:nn), kk+9)  = wt(j,i)*ss.*A(neq+(1:nn), kk+6);              %   s   t^2
            A(neq+(1:nn), kk+10) = wt(j,i)*ctt.*A(neq+(1:nn), kk+6);             %       t^3
    
            A(neq+(1:nn), kk+11) = wt(j,i)*ss.*A(neq+(1:nn), kk+7);              %   s^4 
            A(neq+(1:nn), kk+12) = wt(j,i)*ss.*A(neq+(1:nn), kk+8);              %   s^3 t
            A(neq+(1:nn), kk+13) = wt(j,i)*ss.*A(neq+(1:nn), kk+9);              %   s^2 t^2
            A(neq+(1:nn), kk+14) = wt(j,i)*ss.*A(neq+(1:nn), kk+10);             %   s   t^3
            A(neq+(1:nn), kk+15) = wt(j,i)*ctt.*A(neq+(1:nn), kk+10);            %       t^4
    
            A(neq+(1:nn), kk+16) = wt(j,i)*ss.*A(neq+(1:nn), kk+11);             %   s^5 
            A(neq+(1:nn), kk+17) = wt(j,i)*ss.*A(neq+(1:nn), kk+12);             %   s^4 t
            A(neq+(1:nn), kk+18) = wt(j,i)*ss.*A(neq+(1:nn), kk+13);             %   s^3 t^2
            A(neq+(1:nn), kk+19) = wt(j,i)*ss.*A(neq+(1:nn), kk+14);             %   s^2 t^3
            A(neq+(1:nn), kk+20) = wt(j,i)*ss.*A(neq+(1:nn), kk+15);             %   s   t^4
            A(neq+(1:nn), kk+21) = wt(j,i)*ctt.*A(neq+(1:nn), kk+15);            %       t^5
                        
        end

        b(neq+(1:nn)) = gg;                                                        

        neq = neq+nn;
            
    end 
    
end
x = neq/ng; waitbar(x,handle_waitbar,['generated ...   ', int2str(floor(100*x)), ' %'])
end

dj_disp(['confirm number of equations = ', int2str(neq)])

close(handle_waitbar)

return