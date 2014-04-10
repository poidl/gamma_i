function ydot = adh_diff(x,y,x0,rmax)


%%%    adh_diff          Anderssen-de Hoog differentiator (multiple points)
%%%
%%%    Usage:          ydot = adh_diff(x,y,x0,rmax)
%%%
%%%    Input:          x - input abscissa (equally spaced)
%%%                    y - input ordinates of data
%%%                    x0 - points at which derivatives are desired
%%%                    rmax - maximum value of r (integral)
%%%
%%%    Output:         ydot - estimates of derivative
%%%
%%%    Author:         David Jackett
%%%
%%%    Date:           7/2/97
%%%


n = length(x); 

if length(find(diff(diff(x))>1.0e-5)) ~= 0,
   error('adh_diff error: x is not regular'); end

rmax = round(rmax); nx0 = length(x0);

ydot = nan*ones(nx0,1);

for i = 1:nx0
  ydot(i) = adh_d0(x(:),y(:),x0(i),rmax);
end



return

