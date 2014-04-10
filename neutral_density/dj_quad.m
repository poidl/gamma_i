function ysum = dj_quad(x,y)

%%%    dj_quad 		- trapezoidal quadrature estimate of integral
%%%
%%%    Usage:        ysum = dj_quad(x,y)
%%%
%%%    Input:        x - x values (including NaNs)
%%%                  y - y values (including NaNs)
%%%
%%%
%%%    Output:       ysum - integral sum
%%%
%%%    Author:       David Jackett
%%%
%%%    Date:         8/7/98
%%%


x = x(:); y = y(:);

inds = find(~isnan(x+y)); n = length(inds);

if n<=1, error('error in dj_quad: no valid data'); end

x1 = x(inds); y1 = y(inds);

dx = diff(x1); ymid = (y1(1:n-1)+y1(2:n))/2;

ysum = sum(ymid.*dx);


return

