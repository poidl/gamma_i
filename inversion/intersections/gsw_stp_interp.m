function [SA0,CT0] = gsw_stp_interp(SA,CT,p,p0)

% gsw_stp_interp                                      pressure interpolated
%==========================================================================
% 
% USAGE:         
%  [SA0,CT0] = stp_interp(SA,CT,p,p0)
%  
% DESCRIPTION:
%  Linearly interpolate Absolute Salinity and Conservative Temperature on a 
%  cast to a specified pressure.
%
% INPUT:          
%  SA             cast Absolute Salinities
%  CT             cast Conservative Temperatures
%  p              cast pressures
%  p0             specified pressure
% 
% OUTPUT:
%  SA0            interpolated Absolute Salinity
%  CT0            interpolated Conservative Temperature
%
%  DRJ on 17/06/03
%
%==========================================================================

%  Find the index of a scalar in a monotonically increasing array
n = length(p);
% k = NaN;
if p(1) < p0 && p0 < p(n)
    I_p = find(p >= p0);
    k = I_p(1) - 1;
elseif p0 <= p(1)
    k = 1;
elseif p0 >= p(n)
    k = n - 1;
else
    SA0 = NaN;
    CT0 = NaN;
    return
end

r = (p0 - p(k))/(p(k+1) - p(k));
SA0 = SA(k) + r*(SA(k+1) - SA(k));
CT0 = CT(k) + r*(CT(k+1) - CT(k));

end