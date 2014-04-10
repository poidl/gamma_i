function [depth,inds_depth] = dj_depth(longs,lats,p,z)

%%%    dj_depth          ocean bathymetry from 3D array
%%%
%%%    Usage:         [depth,inds_depth] = dj_depth(longs,lats,p,z)
%%%
%%%    Input:         longs, lats, p - vectors 
%%%                   z              - 3D array
%%%
%%%    Output:        depth          -  array of ocean depths
%%%                   inds_depth     -  array of ocean depth indexes
%%%
%%%    Author:        David Jackett
%%%
%%%    Date:          9/7/98
%%%


jmod = 5;

nx = length(longs); ny = length(lats);

depth = NaN*ones(ny,nx); inds_depth = NaN*ones(ny,nx);

for j = 1:ny
   
   inds = find(~isnan(z(1,j,:))); n = length(inds);
   
   if n > 0,
      for i = 1:n  
         iz = find(~isnan(z(:,j,inds(i)))); 
         depth(j,inds(i)) = p(iz(length(iz)));
				 inds_depth(j,inds(i)) = iz(length(iz));
      end
   end
               
   if mod(j,jmod) == 0 | j == ny+1    
		if length(get(gcf,'CurrentAxes'))>0, close(gcf); end
      h = gcf; dj_pltmp(longs,lats,depth);      
      dj_toc; figure(gcf); dj_pause(1)
   end

end



return

