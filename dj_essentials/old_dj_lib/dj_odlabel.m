function nu = dj_odlabel(s,t,p,longs,lats,tmpq)

%%%
%%%    dj_odlabel        label a 3D global ocean with orthobaric density
%%%
%%%    Usage:            nu = dj_odlabel(s,t,p,longs,lats,tmpq)
%%%
%%%    Input:            s       - 3D salt array               {nz,ny,nx}
%%%                      t       - 3D temperature array        {nz,ny,nx}
%%%                      p       - 1D pressure array           {nz}
%%%                      longs   - vector of longitudes        {nx}
%%%                      lats    - vector of latitudes         {ny}
%%%                      tmpq    - 'tmp' or 'ptmp' to indicate
%%%                                nature of temperature array
%%%                                (optional, defaulting to 'tmp')
%%%
%%%    Output:           nu      - 3D nu value                 {nz,ny,nx}
%%%                      
%%%    Author:           David Jackett
%%%
%%%    Date:             8/8/00
%%%


iplot = 1; jmod = 1;	pr0 = 0;

if nargin==5
	tmpq = 'tmp'; disp(' ') 
	disp('*** temperature assumed to be in-situ ***'); disp(' ')
end

if length(size(s+t)) ~= 3, error('arrays not 3D in dj_odlabel'), end 


[nz,ny,nx] = size(s);

nu = NaN*ones(size(s)); 

 
for j = 1:ny
   
   inds = find(~isnan(t(1,j,:))); n = length(inds);
   
   if n>0;
      
		for i = 1:n

			ss = s(:,j,inds(i));
			tt = t(:,j,inds(i));

			if strcmp(lower(tmpq),'ptmp'),
				tt = sw_ptmp(ss,tt,pr0,p);
			end 

      	    [NU, PHI, PSI, NSQR, DELTAC] = orthobaric(ss,tt,p);
            
            nu(:,j,inds(i)) = ones(size(NU))./NU-1000; 

        end
  
        if iplot==1 & (mod(j,jmod) == 0 | j == ny)
		    nu1 = reshape(nu(1,:),ny,nx); nu1 = change(nu1,'<=',20.0,20.0);
		    dj_pltmp(longs,lats,nu1,0);     
      	    figure(gcf); dj_toc; dj_pause(1)
        end

   end
 
end

return


