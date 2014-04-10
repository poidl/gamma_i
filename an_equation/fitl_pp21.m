global s ct t p g longs lats ocean n

[A,b] = Agen_pp21(s,ct,p,g,ocean,longs,lats);

[B1,c1] = Bgen_pp21(1); [B2,c2] = Bgen_pp21(2);

[B3,c3] = Bgen_pp21(3); [B4,c4] = Bgen_pp21(4);

B = [B1; B2; B3; B4]; c = [c1; c2; c3; c4];

A = [A;10^2*B]; b = [b;c];

x3 = pinv(A)*b; reshape(x3,21,4)

save gamma_pp21.dat x3 -ascii -double

check_pp_contours

