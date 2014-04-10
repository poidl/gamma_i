function k = indx(x,x0)
	  
%%  Find the index of a scalar in a monotonically increasing array
%%  
%%  Usage :      k = indx(x,x0)
%%
%%  Input :      x        array of increasing values
%%               x0       scalar
%%  
%%  Output :     k        when x(k) <= x0 < x(k+1), or
%%               n-1      when x0 = x(n)
%%               n        when x0 > x(n)
%%               1    	  when x0 <= x(1)

%%  DRJ on 17/06/03
  
[Inn]=find(~isnan(x));
n = length(x(Inn));

if x((Inn(1))<x0&x0<x(Inn(n)))
  inds_x = find(x(Inn)>=x0);
  k = inds_x(1)-1;
elseif x0<=x(Inn(1))
  k = 1;
elseif x0==x(Inn(n))
  k = n-1;
elseif x0>x(Inn(n))
  k = n;
else
  disp('ERROR in indx.m : not possible!')
  x, x0
end


return