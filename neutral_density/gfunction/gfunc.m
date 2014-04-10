function g = gfunc(s,ct)

load gamma_coeffs_ct

ss = s(:); ctt = ct(:);

A(:,1) = ones(length(ss),1);
A(:,2) = ctt;
A(:,3) = ctt.^2;
A(:,4) = ctt.^3;
A(:,5) = ss;
A(:,6) = ss.*ctt;
A(:,7) = ss.*ss;

A(:,8) = ctt;
A(:,9) = ctt.^2;
A(:,10) =  ctt.^3;
A(:,11) =  ctt.^4;
A(:,12) = ss;
A(:,13) =  ss.*ctt;
A(:,14) =  ss.*ctt.^3;
A(:,15) = ss.*sqrt(ss);
A(:,16) = ss.*sqrt(ss).*ctt.^2;

gnum = A(:,1:7)*gcoeffs(1:7);

gden = 1+A(:,8:16)*gcoeffs(8:16);

g = gnum./gden; g = g-1000;


return