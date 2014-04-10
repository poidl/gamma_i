function N2 = dj_bfrq(s,t,p)

%%% STABILIZE:      Compute buoyancy frequency profile for cast of data
%%%
%%% USAGE:          N2 = dj_bfrq(s,t,p)
%%%
%%% INPUT:          s       vector of salinities 
%%%                 t       vector of in-situ temperatures 
%%%                 p       vector of pressures
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
%%% REVISION:       
%%%

g = 9.8;

if nargin ~= 3
  error('invalid  # input arguments in dj_bfrq')
end

[n,ncasts] = size(s);

ct = ct_from_t(s,t,p);

s_mid = 0.5*(s(1:n-1,:)+s(2:n,:));
ct_mid = 0.5*(ct(1:n-1,:)+ct(2:n,:));
p_mid = 0.5*(p(1:n-1,:)+p(2:n,:));

[rho,rho_s,rho_ct,rho_p] = eosall_from_ct(s_mid,ct_mid,p_mid);

alfa_mid = -rho_ct./rho;
beta_mid = rho_s./rho;

alfa_ct_z = -alfa_mid.*diff(ct)./diff(p);
beta_s_z = -beta_mid.*diff(s)./diff(p);

N2 = g*(alfa_ct_z - beta_s_z);


return