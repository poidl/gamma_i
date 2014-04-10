function rms = rfunc_7_9_ss(ss,ctt,gg,coeffs)

global s t p g g_rf longs lats ocean n

% global f_called indss ss tt ctt pp gg gg_rf normalizeQ

% if normalizeQ==1, ss = ss/40; ctt = ctt/30; gg = gg/30;       DON'T NEED
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
A(:,1) = ones(length(ss),1);
A(:,2) = ctt;
A(:,3) = ctt.^2;
A(:,4) = ctt.^3;
A(:,5) = ss;
A(:,6) = ss.*ctt;
A(:,7) = ss.*ss;

A(:,8) = ctt;
A(:,9) = ctt.^2;
A(:,10) = ctt.^3;
A(:,11) = ctt.^4;
A(:,12) = ss;
A(:,13) = ss.*ctt;
A(:,14) = ss.*ctt.^3;
A(:,15) = ss.*sqrt(ss);
A(:,16) = ss.*sqrt(ss).*ctt.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gnum = A(:,1:7)*coeffs(1:7);

gden = 1+A(:,8:16)*coeffs(8:16);

gg_rf = gnum./gden; %g_rf = g; g_rf(indss) = gg_rf;

rms = sqrt(mean((gg-gg_rf).*(gg-gg_rf)));

%ss = ss*40; ctt = ctt*30;                              DON'T NEED


return