function [N2,p_mid] = N2_from_ct(s,ct,p)

%%% STABILIZE:      Compute buoyancy frequency profile for cast of data
%%%
%%% USAGE:          N2 = N2_from_ct(s,ct,p)
%%%
%%% INPUT:          s       section of salinities 
%%%                 ct      section of conservative temperatures 
%%%                 p       section of pressures
%%%                 
%%%                 NOTE:   missing values denoted by NaN's
%%%
%%% OUTPUT:         N2      section of N2 values
%%%
%%% UNITS:          salinity                 psu (IPSS-78)
%%%                 temperature              deg. C (ITS-90)
%%%                 pressure                 db
%%%                 N2                       sec^-2
%%%
%%%
%%% AUTHOR:         David Jackett
%%%
%%% CREATED:        April, 2005
%%%
%%% REVISION:       
%%%

g = 9.8;

if nargin ~= 3
  error('invalid  # input arguments in N2_from_ct')
end

[n,ncasts] = size(s);

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