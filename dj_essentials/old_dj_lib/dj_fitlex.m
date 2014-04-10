function [alfa,beta] = dj_fitlex(x,y,x0,xeps,iplot)

%%%    Usage:    [alfa,beta] = dj_fitlex(x,y,x0)
%%%
%%%    Description:     fit a local exponential to y=f(x) at x0
%%%
%%%    Input:           x         	- x vector
%%%                     y	       	- y vector
%%%                     x0          - scalar
%%%
%%%    Output:          alfa	      - multiplicative constant
%%%                     beta        - exponential constant
%%%
%%%    Author:          David Jackett
%%%
%%%    Date:            2/4/98
%%%


n = length(x); x = x(2:n); y = y(2:n);

if iplot == 1, h1 = figure; plot(x,y); hold on; grid on;
               h2 = figure; lny = log(y); plot(x,lny); grid on; end

%%			find a local grid

inds = find(x0-xeps<x&x<x0+xeps); 

x = x(inds); y = y(inds);

n = length(x);

%%
%%      fit straight line in log space
%%


Y = log(y); X = [ones(1,n);x']';

line = X \ Y;

%%
%%      overplot fit
%%

alfa = exp(line(1)); beta = -line(2);

x1 = x(1):0.1:x(n);

fit = alfa*exp(-beta*x1);

if iplot ==1, figure(h1); hold on; ...
					  plot(x1,fit,'r'); grid on; zoom on; end
					  
return

