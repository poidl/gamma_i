function [SA0,CT0] = stp_interp(SA,CT,p,p0)

%  Linearly interpolate salinity and conservative temperature on a cast 
%  to a specified pressure
%
%  Usage :         [SA0,CT0] = stp_interp(SA,CT,p,p0)
%
%  Input :         SA             cast Absolute Salinities
%                  CT             cast Conservative Temperatures
%                  p             cast pressures
%                  p0            pressure level
%
%  Output :        s0            interpolated Absolute Salinity
%                  t0            interpolated Conservative Temperature
%
%  Units :         SA      g/kg (TEOS-10)
%                  CT   degrees C (IPS-90)
%                  pressure      db
%                  density       kg m-3
%  DRJ on 17/06/03


k = indx(p,p0);

r = (p0 - p(k))/(p(k+1) - p(k));

SA0 = SA(k) + r*(SA(k+1) - SA(k));

CT0 = CT(k) + r*(CT(k+1) - CT(k));


return