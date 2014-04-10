function [sns,tns,pns,dsns,dtns,dpns] = ...
								dj_nsfcs(s,t,p,g,longs,lats,glevels,tmpq)

%%%    dj_nsfcs          find neutral surfaces in a 3D global ocean	
%%%
%%%    Usage:            [sns,tns,pns,dsns,dtns,dpns] = ...
%%%													dj_nsfcs(s,t,p,g,longs,lats,glevels,tmpq)
%%%
%%%    Input:            s      - 3D salt array, (nz,ny,nx)
%%%                      t      - 3D temperature array
%%%                      p      - 3D pressure array
%%%                      g      - 3D gamma array
%%%                      longs  - vector of longitudes
%%%                      lats   - vector of latitudes
%%%                      tmpq   - 'tmp' or 'ptmp' to indicate
%%%                                nature of temperature vector
%%%
%%%    Output:           sns    - 3D salinity values, (nz,ny,nx)
%%%                      tns    - 3D temperature values
%%%                      pns    - 3D pressure values
%%%                      dsns   - 3D salinity errors
%%%                      dtns   - 3D temperature errors
%%%                      dpns   - 3D pressure errors
%%%
%%%    Author:           David Jackett
%%%
%%%    Date:             3/12/98
%%%


jmod = 1;	pr0 = 0;


if length(size(s+t+p)) ~= 3
	error('arrays not 3D in dj_nsfcs')
end 

[nz,ny,nx] = size(s); nlevels = length(glevels);

level = floor(nlevels/2); if level==0, level=1; end


%%
%%			find neutral surface locations
%%

sns = NaN*ones(nlevels,ny,nx); dsns = sns;
tns = sns;							 dtns = sns;
pns = sns;							 dpns = sns;

i64 = find(lats<=64); nlat64 = length(i64);

for j = 1:nlat64
   
   inds = find(~isnan(g(1,j,:))); n = length(inds); j
   
   if n>0,
      
      ss = reshape(s(:,j,inds),nz,n);
		tt = reshape(t(:,j,inds),nz,n);
		pp = reshape(p(:,j,inds),nz,n);
		gg = reshape(g(:,j,inds),nz,n);

		if strcmp(lower(tmpq),'ptmp'),
			tt = sw_ptmp(ss,tt,pr0,pp);
		end 
      
      indsn = find(isnan(gg)); nn = length(indsn);      
      
      ss(indsn) = NaN; tt(indsn) = NaN;
     
		[sns1,tns1,pns1,dsns1,dtns1,dpns1] = nsfcs(ss,tt,pp,gg,glevels); 
 
      sns(:,j,inds) = sns1; dsns(:,j,inds) = dsns1;
      tns(:,j,inds) = tns1; dtns(:,j,inds) = dtns1;
      pns(:,j,inds) = pns1; dpns(:,j,inds) = dpns1;

      if mod(j,jmod) == 0 | j == nlat64
      	h = gcf-1; z = reshape(pns(level,:),ny,nx);
			z = change(z,'==',-99,NaN);
         dj_pltmp(longs,lats,z);
      	dj_toc; figure(gcf); dj_pause(1);         
      end

   end 
end
 

return


