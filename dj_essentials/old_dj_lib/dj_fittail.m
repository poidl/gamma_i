function [mu,lambda] = fit_tail(x,greensf)

%%%    Usage:    [lambda,mu] = fit_tail(x,greensf)
%%%
%%%    Description:     fit an exponential to the tail of the ocean's
%%%                     Green's function
%%%
%%%    Input:           x           - time axis
%%%                     greensf     - ocean's Green's function
%%%
%%%    Output:          lambda      - multiplicative constant
%%%                     mu          - exponential constant
%%%
%%%    Author:          David Jackett
%%%
%%%    Date:            25/6/96
%%%


x0 = 380; iplot = 0;

n = length(x);
 
if x(n) < 400, error('error in fit_tail.m: Greens function too short'); end

k0 = find(x==x0);

x = x(k0:n); y = greensf(k0:n);

lny = log(y);

if iplot ==1, figure; plot(x,lny); end

%%
%%      fit straight line in log space
%%

Y = lny'; X = [ones(size(x))',x'];

line = X \ Y;

%%
%%      overplot fit
%%

mu = exp(line(1)); lambda = -line(2);

fit = mu*exp(-lambda*x);

if iplot ==1, figure; plot(x,y); hold on; plot(x,fit,'c'); end

return

