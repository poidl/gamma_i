function [s_stab,t_stab] = dj_stabilize1(s,t,p,longs,lats,N2)

%%%
%%% STABILIZE:      Stabilize 3D hydrographic data w.r.t. buoyancy frequency N2
%%%
%%% USAGE:          [ss,tt] = dj_stabilize1(s,t,p,N2)
%%%
%%% INPUT:          s       3D array of salinities (each column being a cast)
%%%                 t       3D array of in-situ temperatures 
%%%                 p   	 3D array of pressures
%%%                 N2      vector of N2 lower bounds
%%%
%%%                 NOTE:   missing values denoted by NaN's
%%%
%%% OUTPUT:         s_stab    	3D array of adjusted salinities
%%%                 t_stab      	3D array of adjusted in-situ temperatures
%%%
%%% UNITS:          salinity    				psu (IPSS-78)
%%%                 temperature 			 	deg. C (IPTS-68)
%%%                 pressure    				db
%%%                 N2     					sec^-2
%%%
%%%
%%% AUTHOR:         David Jackett
%%%
%%% CREATED:        May, 2000
%%%
%%% REVISION:       1.1     22/3/97
%%%

global location

iplot = 1; jmod = 1;	pr0 = 0;


%%%
%%%     						some checks
%%%

if nargin ~= 6
  error('invalid input arguments in dj_stabilize')
end

if length(size(s+t+p)) ~= 3, error('arrays not 3D in dj_stabilize'), end 


[nz,ny,nx] = size(s);

s_stab = NaN*ones(size(s)); t_stab = s_stab; trouble2 = nan*ones(ny,nx);

  
for j = 1:ny,
   
   inds0 = find(~isnan(s(1,j,:))); n = length(inds0);
   
   for i = 1:n;  %  if n > 0

		inds = inds0(i); location = [j,i,inds];

%		ss = reshape(s(:,j,inds),nz,n);
%		tt = reshape(t(:,j,inds),nz,n);
%		pp = reshape(p(:,j,inds),nz,n);

		ss = reshape(s(:,j,inds),nz,1);
		tt = reshape(t(:,j,inds),nz,1);
		pp = reshape(p(:,j,inds),nz,1);

%%				N2 computation

      BV_freq = sw_bfrq(ss,tt,pp);


%%				troublesome ones

		in_trouble = BV_freq(:)<N2 & length(find(isfinite(ss)))>3;

		in_trouble = sum(in_trouble);

		if in_trouble>0
	
			if iplot==2, figure(1), plot(ss,tt,'r'); grid on; hold on; end

			out = [ss,tt,pp]; save eg.dat -ascii out

			[ssss,tttt] = stabilize(ss,tt,pp);

			s_stab(:,j,inds) = ssss;
			t_stab(:,j,inds) = tttt;

			trouble2(j,inds) = in_trouble;

			if iplot == 2, figure(1), plot(ssss,tttt), hold off, dj_pause(0), end

%%				check it has stabilized it
	      BV_freq1 = sw_bfrq(ssss,tttt,pp);

			in_trouble1 = BV_freq1(:)<0*N2 & length(find(isfinite(ss)))>3;

			in_trouble1 = sum(in_trouble1);

         if in_trouble1>0
				dj_disp('***	still in trouble	***')
				figure(2), hold on, plot(longs(inds),lats(j),'m.'), dj_pause(1), hold off
%				[BV_freq,BV_freq1]
%				dj_pause(0)
			end

		end

   end
  
   if iplot==1 & (mod(j-1,jmod) == 0 | j == ny)
			figure(2), z = reshape(trouble2,ny,nx);
			dj_pltmp(longs,lats,z,0);
			title(['j = ', int2str(j), '  :  lat = ', num2str(lats(j))])     
      	figure(gcf); dj_toc; dj_pause(1)
			save st_data s_stab t_stab
   end

end



return


