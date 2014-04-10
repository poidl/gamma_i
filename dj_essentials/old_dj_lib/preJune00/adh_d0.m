function ydot = adh_d0(x,y,x0,rmax)

%%%    adh_d0		Anderssen-de Hoog differentiator (single point)
%%%
%%%    Usage:    ydot = adh_d0(x,y,x0,rmax)
%%%
%%%    Description:    Anderssen - de Hoog data differentiator
%%%
%%%    Input:          x - input abscissa (equally spaced)
%%%                    y - input ordinates of data
%%%                    x0 - point at which derivative is desired
%%%                    rmax - maximum value of r (integer)
%%%
%%%    Output:         ydot - estimate of derivative
%%%
%%%    Author:         David Jackett
%%%
%%%    Date:           12/6/96
%%%


n = length(x); 

ind = find(x==x0); h = mean(diff(x));

if ind < 1 | ind > n, error('error in adh_diff.m: not a valid point'); end

if ind==1, ydot =  (y(2)-y(1))/h; return, end

if ind==n, ydot =  (y(n)-y(n-1))/h; return, end

rhi = min([ind-1,n-ind,rmax]);

for r = 1:rhi
  Kr = r*(r+1)*(2*r+1)/6; ylo = 0; yhi = 0;
  for k = 1:r
    wk = k/(2*Kr*h); ylo = ylo+wk*y(ind-k); yhi = yhi+wk*y(ind+k);
  end
  ydot = yhi-ylo;
end


return

