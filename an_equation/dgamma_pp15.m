function dgg_pp = dgamma_pp15(ss0,ctt0,alpha0,beta0,c)


c1 = c(1:15); c2 = c(16:30); c3 = c(31:45); c4 = c(46:60);

dgg_p1 = dgamma_p15(ss0,ctt0,alpha0,beta0,c1);

dgg_p2 = dgamma_p15(ss0,ctt0,alpha0,beta0,c2);

dgg_p3 = dgamma_p15(ss0,ctt0,alpha0,beta0,c3);

dgg_p4 = dgamma_p15(ss0,ctt0,alpha0,beta0,c4);


dgg_pp = dgg_p1 + dgg_p2 + dgg_p3 + dgg_p4;


return