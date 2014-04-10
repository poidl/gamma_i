 function [SAns,CTns,pns,dSAns,dCTns,dpns] = gsw_neutral_surfaces(SA,CT,p,gamma,glevels)
 
% 
%--------------------------------------------------------------------------
% 
% Usage
%  [SAns,CTns,pns,dsns,dtns,dpns) = neutral_surfaces(SA,CT,p,gamma,glevels)   
%
% DESCRIPTION
%  For a section of hydrographic data (SA, CT, p) that has been labelled 
%  with gamma, find the Absolute Salinities, Conservative Temperatures and
%  pressures of neutral density surfaces specified by glevels.
%
% Input        
%  SA          Absolute Salinity of bottles
%  CT          Conservative Temperature of bottles
%  p           pressure of bottles
%  gamma       gamma values of bottles
%  glevels     neutral density surfaces to be outputed
%
% Output 
%  SAns        Absolute Salinity on the neutral density surfaces
%  CTns        Conservative Temperature on the neutral density surfaces
%  pns         pressure of the neutral density surfaces
%  dSAns       Absolute Salinity errors
%  dCTns       Conservative Temperature errors
%  dpns        pressure errors
%
% Note that sns, tns and pns values of nan denotes under or outcropping
%
% Units:          salinity    psu (IPSS-78)
%                 temperature degrees C (IPS-90)
%                 pressure    db
%                 gamma       kg m-3
%
%
% DRJ on 16/06/03
%
%--------------------------------------------------------------------------

[m,n] = size(SA);
ng = length(glevels);

SAns = nan(ng,n); 
CTns = SAns; 
pns = SAns;

dSAns = zeros(ng,n); 
dCTns = dSAns; 
dpns = dSAns;

for k = 1:n
    Idata = find(isfinite(gamma(:,k)));
    if ~isempty(Idata)
      [SAns(:,k),CTns(:,k),pns(:,k),dSAns(:,k),dCTns(:,k),dpns(:,k)] = ...
                        gsw_neutral_surfaces0(SA(Idata,k),CT(Idata,k),p(Idata,k),gamma(Idata,k),glevels);
    end
end

return