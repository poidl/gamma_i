function ysum = dj_quadi(x,y,xlo,xhi)

%%%    dj_quadi             quadrature estimate of integral
%%%
%%%    Usage:         ysum = dj_quadi(x,y,xlo,xhi)
%%%
%%%    Input:         x - x values (including NaNs)
%%%                   y - y values (including NaNs)
%%%                   xlo,xhi - limits
%%%
%%%    Output:        ysum - integral sum
%%%
%%%    Author:        David Jackett
%%%
%%%    Date:          8/7/98
%%%


x = x(:); y = y(:);

inds = find(~isnan(x+y)); n = length(inds);

if n==0, error('Error 1 in dj_quadi'); end


x1 = x(inds); y1 = y(inds);

if xlo<x1(1) | xhi>x1(n), ...
   [xlo, xhi]
   x1
      error('Error 2 in dj_quadi: out of limits range')
end


ylo = interp1(x1,y1,xlo);
yhi = interp1(x1,y1,xhi);

inds = find(xlo<x1 & x1<xhi);

xx = [xlo;x1(inds);xhi]; yy = [ylo;y1(inds);yhi];


ysum = dj_quad(xx,yy);



return

