function [gamma,dgl,dgh] = dj_glabel(s,t,p,longs,lats,tmpq)

%%%    dj_glabel         label a 3D global ocean with gamma
%%%
%%%    Usage:            [gamma,dgl,dgh] = dj_glabel(s,t,p,longs,lats,tmpq)
%%%
%%%    Input:            s      - 3D salt array, (nz,ny,nx)
%%%                      t      - 3D temperature array
%%%                      p      - 3D pressure array
%%%                      longs  - vector of longitudes
%%%                      lats   - vector of latitudes
%%%                      tmpq   - 'tmp' or 'ptmp' to indicate
%%%                                nature of temperature vector
%%%										  (defaults to 'tmp')
%%%
%%%    Output:           gamma  - 3D gamma value, (nz,ny,nx)
%%%                      dgl    - 3D array of lower error tolerances
%%%                      dgh    - 3D array of upper error tolerances
%%%
%%%    Author:           David Jackett
%%%
%%%    Date:             3/12/98
%%%


iplot = 1; jmod = 5;	pr0 = 0;

if nargin==5
	tmpq = 'tmp'; disp(' ') 
	disp('*** temperature assumed to be in-situ ***'); disp(' ')
end

if length(size(s+t+p)) ~= 3, error('arrays not 3D in dj_glabel'), end 


[nz,ny,nx] = size(s);

gamma = NaN*ones(size(s)); dgl = gamma; dgh = gamma;

i80 = find(lats<=-80); nlat80 = length(i80)+1;  
i64 = find(lats<=64); nlat64 = length(i64);
  
for j = nlat80:nlat64,
   
   inds = find(~isnan(t(1,j,:))); n = length(inds);
   
   if n>0;
      
      alongs = longs(inds); alats = lats(j)*ones(n,1);


%%						a fix for Tony's data

%      if ismember(longs(20),alongs) & ismember(lats(27),alats),
%         dj_disp(['fixing dud point ',num2str(longs(20)),'  ', ...
%													num2str(lats(27))]);
%         ind = find(alongs==longs(20)); alongs(ind)=longs(ind+1);
%      end;
    
		ss = reshape(s(:,j,inds),nz,n);
		tt = reshape(t(:,j,inds),nz,n);
		pp = reshape(p(:,j,inds),nz,n);

		if strcmp(lower(tmpq),'ptmp'),
			tt = sw_ptmp(ss,tt,pr0,pp);
		end 

      [gg,dg_l,dg_h] = gamma_n(ss,tt,pp,alongs,alats);
      
   	gg = change(gg,'<=',-99.0,NaN); gamma(:,j,inds) = gg;

   	dg_l = change(dg_l,'<=',-99.0,NaN); dgl(:,j,inds) = dg_l;
   	dg_h = change(dg_h,'<=',-99.0,NaN); dgh(:,j,inds) = dg_h;
          
      if iplot==1 & (mod(j,jmod) == 0 | j == nlat64)
			g1 = reshape(gamma(1,:),ny,nx); g1 = change(g1,'<=',20.0,20.0);
			dj_pltmp(longs,lats,g1,0);     
      	figure(gcf); dj_toc; dj_pause(1)
      end

   end
 
end

return


