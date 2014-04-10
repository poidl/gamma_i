function N2 = bfrq1(s,t,p,tvar)

%%% STABILIZE:      Compute buoyancy frequency squared profile section of data
%%%
%%% USAGE:          N2 = bfrq1(s,t,p,tvar)
%%%
%%% INPUT:          s       vector of salinities 
%%%                 t       vector of in situ temperatures 
%%%                 p       vector of pressures
%%%                 tvar    temperature variable ('ct' or 'theta', default 'ct')
%%%                 
%%%                 NOTE:   missing values denoted by NaN's
%%%
%%% OUTPUT:         N2      vector of N2 values
%%%
%%% UNITS:          salinity                 psu (IPSS-78)
%%%                 temperature degrees      deg. C (ITS-90)
%%%                 pressure                 db
%%%                 N2                       sec^-2
%%%
%%%
%%% AUTHOR:         David Jackett
%%%
%%% CREATED:        Sept, 2002
%%%
%%% REVISION:       June, 2005       
%%%

g = 9.8;

if nargin == 3
  disp('conservative temperature is default'), tvar = 'ct';
elseif (nargin<=2|nargin>4)
  stop('**** invalid input arguments in bfrq ****')
end

if strcmp(tvar,'ct')
    t = ct_from_t(s,t,p);
else
    pr0 = zeros(size(p));
    t = theta_from_t(s,t,p,pr0);
    tvar = 'theta';
end

inds = find(isfinite(s));
n = length(inds);

s_mid = 0.5*(s(1:n-1,:)+s(2:n,:)); %s_mid = s(2:n,:);
t_mid = 0.5*(t(1:n-1,:)+t(2:n,:)); %t_mid = t(2:n,:);
p_mid = 0.5*(p(1:n-1,:)+p(2:n,:)); %p_mid = p(2:n,:);


cmd = ['[rho,rho_s,rho_t,rho_p] = eosall_from_', tvar, '(s_mid,t_mid,p_mid);'];
eval(cmd)

alfa_mid = -rho_t./rho;
beta_mid = rho_s./rho;

alfa_t_z = -alfa_mid.*diff(t)./diff(p);
beta_s_z = -beta_mid.*diff(s)./diff(p);

N2 = g*(alfa_t_z - beta_s_z);

inds1 = find(abs(N2) <= eps);
N2(inds1) = 0;

%   interpolate back onto original grid
N2_h = N2(1);
N2_l = N2(n-1);

N2(2:n-1) = interp1(p_mid,N2(1:n-1),p(2:n-1));

if n > 3
  N2(1) = N2_h-(N2(3)-N2(2))/2; 
  N2(n) = N2_l+(N2(n-1)-N2(n-2))/2;   %  fix for top and bottom
else
  N2(1) = N2(2); 
  N2(3) = N2(2);
end



return