function [s,ct] = add_sctp_data(s0,ct0,p0,deltas)


n=25; s = nan*ones(n,1); ct = s; g = s;

rho_p0 = rho_from_ct(s0,ct0,p0); 

for k = 1:n
    
    ss = s0+k*deltas; 
    
    if k==1
        ctt = ct0;
    else
        ctt = ct(k-1);
    end

    eps_ctt = 1; iter = 0;  denom = 1;
    
    while eps_ctt>1e-5 & denom~=0
        
        iter = iter+1;
        
        rho_diff = rho_from_ct(ss,ctt,p0) - rho_p0; 
        
        denom = (rho_from_ct(ss,ctt+0.05,p0) - rho_from_ct(ss,ctt-0.05,p0)) ./ 0.1;
        
        if denom~=0
        
            ctt_new = ctt - rho_diff./denom;
        
            eps_ctt = abs(ctt_new-ctt);

%         errors = [k, iter, ctt, ctt_new, eps_ctt]
        
            ctt_last = ctt; ctt = ctt_new;
            
        end 
       
    end
    
    if ctt>fp_ct(ss,500)-5
        s(k) = ss; ct(k) = ctt;
    else
        k = n;
    end
       
end

inds = find(s>0&finite(ct)); s = s(inds); ct = ct(inds);



