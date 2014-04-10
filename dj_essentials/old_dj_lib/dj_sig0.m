function [ssfc,tsfc,psfc] = dj_sig0(s,t,p,pr0,sig_p,tmpq)

%%%    dj_sig0           find a sigma_n surface on a cast
%%%
%%%    Usage:         psfc = dj_sig0(s,t,p,pr0,sig_p)
%%%
%%%    Input:         s 		- salt vector
%%%                   t 		- temperature vector
%%%                   p 		- pressure vector
%%%                   pr0 	- reference pressure
%%%                   sig_p - sig_p value
%%%                   tmpq  - 'tmp' or 'ptmp' to indicate
%%%                           nature of temperature vector
%%%
%%%    Output:        ssfc	- salinity on sig_p surface
%%%    		        	 tsfc	- temperature on sig_p surface
%%%    			       psfc	- pressure on sig_p surface
%%%
%%%    Author:        David Jackett
%%%
%%%    Date:          31/7/98
%%%




n = length(find(finite(s)));


if strcmp(lower(tmpq),'tmp')
  sig = sw_pden(s,t,p,pr0)-1000;
elseif strcmp(lower(tmpq),'ptmp')
  sig = sw_dens(s,t,pr0)-1000;
end  


diff = sig-sig_p;

inds = find(diff>=0); ninds = length(inds);  % inds

if ninds ~= 0 & ninds < n
	k = inds(1)-1;
   psfc = interp1(sig(k:k+1),p(k:k+1),sig_p);
elseif ninds == 0										% bottoms out
	inds1 = find(~isnan(sig)); ninds1 = length(inds1);
	psfc = p(inds1(ninds1));
else 														% tops out
   psfc = p(1);
end


%%					interpolate to get temperature and pressure


if n>1
	ssfc = interp1(p,s,psfc);
	tsfc = interp1(p,t,psfc);
else
	ssfc = s;
	tsfc = t;
end
  


return


