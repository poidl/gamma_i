function [ssfc,tsfc,psfc] = dj_sigsfc(s,t,p,longs,lats,pr0,sig_val,tmpq)

%%%    dj_sigsfc         find a sigma_n surface in a 3D ocean
%%%
%%%    Usage:            psfc = dj_sigsfc(s,t,p,pr0,sig_p)
%%%
%%%    Input:            s      - 3D salt array, (nz,ny,nx)
%%%                      t      - 3D temperature array
%%%                      p      - 3D pressure array
%%%                      pr0    - scalar reference pressure
%%%                      sval   - scalar sig_n value
%%%                      tmpq   - 'tmp' or 'ptmp' to indicate
%%%                                nature of temperature vector
%%%
%%%    Output:           psfc   - 2D surface pressure array, (ny,nx)
%%%
%%%    Author:           David Jackett
%%%
%%%    Date:             24/11/98
%%%


jmod = 12; figq = 1;


if length(size(s+t+p)) ~= 3, error('array not 3D in dj_sigsfc'), end 

[nz,ny,nx] = size(s);

psfc = NaN*ones(ny,nx); ssfc = psfc; tsfc = psfc;

for j = 1:ny
   
   inds = find(~isnan(t(1,j,:))); n = length(inds);
   
   if n > 0,     
      for i = 1:n      
         it = find(~isnan(t(:,j,inds(i))));	
         ss = s(it,j,inds(i));
			tt = t(it,j,inds(i));
			pp = p(it,j,inds(i));
         [ssfc(j,inds(i)),tsfc(j,inds(i)),psfc(j,inds(i))] = ...
											dj_sig0(ss,tt,pp,pr0,sig_val,tmpq);
      end
   end
               
   if figq == 1 & (mod(j,jmod) == 0 | j == ny)        
      h = gcf; dj_pltmp(longs,lats,psfc); dj_pause(1)    
      dj_toc; close(h)       
   end

end

