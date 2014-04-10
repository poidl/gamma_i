function [B,c] = Bgen_pp21(k)


ss_b1 = 0:0.5:30; ctt_b1 = -5:0.5:40;                 %   [0,30]x[-5,40]
[ss_b2,ctt_b2] = meshgrid(ss_b1,ctt_b1);
ss0 = ss_b2(:); ctt0 = ctt_b2(:); 

ss_b1 = 30:0.5:42; ctt_b1 = 25:0.5:40;              %   [30,42]x[25,40]
[ss_b2,ctt_b2] = meshgrid(ss_b1,ctt_b1);
ss0 = [ss0; ss_b2(:)]; ctt0 = [ctt0; ctt_b2(:)]; 
     
pp0 = zeros(size(ss0));       
[rho,rhoss,rhoctt,rhopp] = eosall_from_ct(ss0,ctt0,pp0);           
alpha0 = -rhoctt./rho; beta0 = rhoss./rho;

no_boundary_eqns = length(ss0)

B = zeros(no_boundary_eqns,84); c = zeros(no_boundary_eqns,1);

kk = (k-1)*21;

B(:, kk+2)  = alpha0;   
B(:, kk+3)  = beta0;

B(:, kk+4)  = 2*alpha0.*ss0; 
B(:, kk+5)  = alpha0.*ctt0 + beta0.*ss0; 
B(:, kk+6)  = 2*beta0.*ctt0;  

B(:, kk+7)  = 1.5*ss0.*B(:, kk+4);   
B(:, kk+8)  = ctt0.*B(:, kk+4) + beta0.*ss0.*ss0; 
B(:, kk+9)  = alpha0.*ctt0.*ctt0 + ss0.*B(:, kk+6);  
B(:, kk+10) = 1.5*ctt0.*B(:, kk+6);   

B(:, kk+11) = (4/3)*ss0.*B(:, kk+7); 
B(:, kk+12) = ss0.*B(:, kk+8) + ctt0.*B(:, kk+7)/3;   
B(:, kk+13) = ctt0.*(B(:, kk+8) + beta0.*ss0.*ss0);   
B(:, kk+14) = ctt0.*B(:, kk+9) + ss0.*B(:, kk+10)/3;
B(:, kk+15) = (4/3)*ctt0.*B(:, kk+10);   

B(:, kk+16) = 5*alpha0.*ss0.^4;
B(:, kk+17) = 4*alpha0.*ctt0.*ss0.^3 + beta0.*ss0.^4;
B(:, kk+18) = 3*alpha0.*ctt0.^2.*ss0.^2 + 2*beta0.*ctt0.*ss0.^3;
B(:, kk+19) = 2*alpha0.*ctt0.^3.*ss0 + 3*beta0.*ctt0.^2.*ss0.^2;
B(:, kk+20) = alpha0.*ctt0.^4 + 4*beta0.*ctt0.^3.*ss0;
B(:, kk+21) = 5*beta0.*ctt0.^4;


return