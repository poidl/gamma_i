function [s_stab,t_stab] = dj_N2(s,t,p,longs,lats,N2)

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

global location s_stab t_stab

global s0 t0 p0

iplot = 2; jmod = 1;	pr0 = 0;


%%%
%%%     						some checks
%%%

if nargin ~= 6
  error('invalid input arguments in dj_stabilize')
end

if length(size(s+t+p)) ~= 3, error('arrays not 3D in dj_stabilize'), end 


[nz,ny,nx] = size(s);

s_stab = NaN*ones(size(s)); t_stab = s_stab; trouble2 = nan*ones(ny,nx);

nit = 0;

for j = 5:ny,
   
   inds0 = find(~isnan(s(1,j,:))); n = length(inds0);
   
   for i0 = 12:n;  %  if n > 0

		i = inds0(i0); location = [j,i0,i]

		ss = reshape(s(:,j,i),nz,1);
		tt = reshape(t(:,j,i),nz,1);
		pp = reshape(p(:,j,i),nz,1);

%%				N2 computation

      BV_freq = sw_bfrq(ss,tt,pp);


%%				troublesome ones

%		[size(BV_freq),size(N2)]

		in_trouble = BV_freq(:)<N2(:) & length(find(isfinite(ss)))>3;

		in_trouble = sum(in_trouble);

		if in_trouble>0
	
   if iplot==2
      figure(1)
		subplot(2,2,1)
      	plot(ss,tt); grid on; zoom on
 			title(['j = ', int2str(j), '  :  i = ', num2str(i)])     
      subplot(2,2,3)
         plot(s0(:,j,i),t0(:,j,i),'r'); grid on; hold on
      	plot(ss,tt); zoom on, hold off	
      subplot(2,2,4)
      	pmid = (pp(1:nz-1)+pp(2:nz))/2;
      	indsp = find(BV_freq>0);
         plot(log10(BV_freq(indsp)),pmid(indsp)/1000,'.')
         hold on, grid on
         plot(log10(BV_freq(indsp)),pmid(indsp)/1000)
         inds2 = find(BV_freq<0); 
         plot(-10*ones(size(inds2)),pmid(inds2)/1000,'*')
         set(gca,'ydir','reverse')
         dj_pause(1)
   end

%			out = [ss,tt,pp]; save eg.dat -ascii out

			[ssss,tttt] = stabilize(ss,tt,pp);

			s_stab(:,j,i) = ssss;
			t_stab(:,j,i) = tttt;

			trouble2(j,i) = in_trouble;

	      BV_freq1 = sw_bfrq(ssss,tttt,pp);

			if iplot == 2
   			figure(1), subplot(2,2,1),hold on
            	plot(ssss,tttt,'k'), hold off, zoom on
					title(['j = ', int2str(j), '  :  i = ', num2str(i)])     
				subplot(2,2,4), hold on
      			indsp = find(BV_freq1>0);
            	plot(log10(BV_freq1(indsp)),pmid(indsp)/1000,'k'), 
            	inds2 = find(BV_freq1<0); 
               plot(-9*ones(size(inds2)),pmid(inds2)/1000,'k*')
               hold off, zoom on, dj_pause(1)
         end
            
%%				check it has stabilized it

			in_trouble1 = BV_freq1(:)<0*N2(:) & length(find(isfinite(ss)))>3;

			in_trouble1 = sum(in_trouble1);

			if in_trouble1>0
   			nit = nit+1;
            dj_disp(['***	still in trouble ',int2str(nit)])
            N2_min = [inds2',nanmin(BV_freq),nanmin(BV_freq1)]
				dj_pause(0)
			end

		end

   end
  
   if iplot>=1 & (mod(j-1,jmod) == 0 | j == ny)
      figure(1), subplot(2,2,2)
      	z = reshape(trouble2,ny,nx);
			dj_pltmp(longs,lats,z,0);
			title(['j = ', int2str(j), '  :  lat = ', num2str(lats(j))])     
      	figure(gcf); dj_toc; dj_pause(1)
   end

end



return


