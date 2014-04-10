function [lambda,mu] = dj_fitex(x,y)

%%%    Usage:    [lambda,mu] = dj_fitex(x,y)
%%%
%%%    Description:     fit an exponential to y=f(x)
%%%
%%%    Input:           x           - x axis
%%%                     y	        - y axis
%%%
%%%    Output:          lambda      - multiplicative constant
%%%                     mu          - exponential constant
%%%
%%%    Author:          David Jackett
%%%
%%%    Date:            25/3/97
%%%


iplot = 0;

x = x(:); y = y(:);

inds = find(y<=0); y(inds) = eps*ones(length(inds),1);

n = length(x);

lny = log(y);

if iplot ==1, figure; plot(x,lny); end

%%
%%      fit straight line in log space
%%

Y = lny; X = [ones(1,n);x']';

line = X \ Y;

%%
%%      overplot fit
%%

lambda = exp(line(1)); mu = -line(2);

x1 = x(1):0.1:x(n);

fit = lambda*exp(-mu*x1);

if iplot ==1, figure; plot(x,y); hold on; ...
					  plot(x1,fit,'c'); dj_pause(5); end
					  

return

