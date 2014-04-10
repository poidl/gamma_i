function conv = dj_conv(g,tg,f,tf)

%%%    dj_conv		convolution of Greens function and data
%%%
%%%    Usage:			conv = dj_conv(g,tg,f,tf)
%%%
%%%    Input:        g - Greens function
%%%                  tg - time ordinates of the g data
%%%                  f - forcing array
%%%                  tf - time ordinates of the f data
%%%
%%%    Output:       conv - convolution integral of g and f
%%%
%%%    Assumptions:  g starts from time t=0
%%%                  g is longer than f in time
%%%                  the length of time elapsed in f appears in g
%%%
%%%    Author:       David Jackett
%%%
%%%    Date:         20/3/98


g = g(:); tg = tg(:); f = f(:); tf = tf(:);

ng = length(g); nf = length(f); dtf = tf(nf)-tf(1); n = find(tg==dtf); 


if tg(1) ~= 0; error('**** NOT A GREENS FUNCTION IN DJ_CONV ****'); end

if tg(ng) < dtf; error('**** GREENS FUNCTION NOT LONG ENOUGH IN DJ_CONV ****'); end

if length(n) ~= 1; error('**** RESOLUTIONS WRONG IN DJ_CONV ****'); end 


tt = tf(1)+tg(1:n); gg = g(1:n); ff = interp1(tf,f,tt);

yy = flipud(gg).*ff;

conv =  trapz(tt,yy);


return
