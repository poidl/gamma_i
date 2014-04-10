function N2 = dj_bfrq(s,t,p)

%%% STABILIZE:      Compute buoyancy frequency profile for 3D ocean
%%%
%%% USAGE:          N2 = dj_bfrq(s,t,p)
%%%
%%% INPUT:          s       matrix of salinity (each column being a cast)
%%%                 t       matrix of in-situ temperatures 
%%%                 p   	 matrix of pressures
%%%                 N2      vector of N2 lower bounds
%%%
%%%                 NOTE:   missing values denoted by NaN's
%%%
%%% OUTPUT:         ss    	 matrix of adjusted salinities
%%%                 tt      matrix of adjusted in-situ temperatures
%%%
%%% UNITS:          salinity    				psu (IPSS-78)
%%%                 temperature degrees 	deg. C (IPTS-68)
%%%                 pressure    				db
%%%                 N2     					sec^-2
%%%
%%%
%%% AUTHOR:         David Jackett
%%%
%%% CREATED:        May, 2000
%%%
%%% REVISION:       1.1     6/3/97
%%%


iplot = 1; jmod = 1;	pr0 = 0;


%%%
%%%     						some checks
%%%

if nargin ~= 4
  error('invalid input arguments in dj_stabilize')
end

if length(size(s+t+p)) ~= 3, error('arrays not 3D in dj_stabilize'), end 


[nz,ny,nx] = size(s);

ss = NaN*ones(size(s)); tt = ss;

  
for j = 1:ny,
   
   inds = find(~isnan(s(1,j,:))); n = length(inds);
   
   if n>0;
      
		ss = reshape(s(:,j,inds),nz,n);
		tt = reshape(t(:,j,inds),nz,n);
		pp = reshape(p(:,j,inds),nz,n);

      [sss,ttt] = stabilize(ss,tt,pp,N2);
          
      if iplot==1 & (mod(j,jmod) == 0 | j == ny)
			z = reshape(t(1,:),ny,nx);
			dj_pltmp(longs,lats,z,0);     
      	figure(gcf); dj_toc; dj_pause(1)
      end

   end
 
end



return


