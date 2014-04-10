function [gc1,gc2,gc3,gc4,gc5,nxz_beqs] = bxz_equations(s,ct,p,g,ocean,n,longs,lats,handles)

[nz,ny,nx] = size(s);

inds_g = find(isfinite(g)); ng = length(inds_g)

nr = zeros(ng,1); gc1 = nan*ones(ng,1); gc2 = gc1; gc3 = gc1; gc4 = gc1; gc5 = gc1;

nxz_beqs = 0;
    
for kg = 1:ng
    
  [k,j,i] = ind2sub([nz,ny,nx],inds_g(kg)); 
  
  if i>1 && i<nx
    
    ie = i+1; iw = i-1;
        
    if isfinite(s(k,j,ie)) && isfinite(s(k,j,iw))
        indss = find(isfinite(s(:,j,i)));
        ss = s(indss,j,i); ctt = ct(indss,j,i); pp = p(indss,j,i);
        if k==1

        elseif k==n(j,i)

        else
            
            nxz_beqs = nxz_beqs+1;
            
%       vertical differences

            pmid = (pp(k-1)+pp(k))/2;
            dsigu = rho_from_ct(ss(k),ctt(k),pmid)-rho_from_ct(ss(k-1),ctt(k-1),pmid); 
            
            pmid = (pp(k)+pp(k+1))/2;
            dsigl = rho_from_ct(ss(k+1),ctt(k+1),pmid)-rho_from_ct(ss(k),ctt(k),pmid);

%       horizontal differences

            pmid = (pp(k)+p(k,j,ie))/2;            
            dsige = rho_from_ct(s(k,j,ie),ct(k,j,ie),pmid)-rho_from_ct(ss(k),ctt(k),pmid);
            
            pmid = (p(k,j,iw)+pp(k))/2;            
            dsigw = rho_from_ct(ss(k),ctt(k),pmid)-rho_from_ct(s(k,j,iw),ct(k,j,iw),pmid);
            
            gcoeffs_x
                       
        end
 

    end
    
    if mod(kg,floor(ng/10))==0
        dj_disp(['bxz_equations completed ',int2str(round(100*kg/ng)),'%'])
        dj_toc
    end
    
  end

end


nxz_eqs = [nxz_beqs,ng,100*nxz_beqs/ng]
         
ok = 'done the x-z b equation', dj_toc

mean_of_coefficients = nanmean([gc1,gc2,gc3,gc4,gc5])
 
         
return
