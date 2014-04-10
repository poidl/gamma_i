function dgg_p = dgamma_p21(ss0,ctt0,alpha0,beta0,coeffs)


nn0 = length(ss0); B = zeros(nn0,21);

B(:,2)  = alpha0;   
B(:,3)  = beta0;
    
B(:,4)  = 2*alpha0.*ss0;                 
B(:,5)  = alpha0.*ctt0 + beta0.*ss0;         
B(:,6)  = 2*beta0.*ctt0;              
    
B(:,7)  = 1.5*ss0.*B(:,4);           
B(:,8)  = ctt0.*B(:,4) + beta0.*ss0.*ss0;         
B(:,9)  = alpha0.*ctt0.*ctt0 + ss0.*B(:,6);  
B(:,10) = 1.5*ctt0.*B(:,6);           
    
B(:,11) = (4/3)*ss0.*B(:,7);         
B(:,12) = ss0.*B(:,8) + ctt0.*B(:,7)/3;           
B(:,13) = ctt0.*(B(:,8) + beta0.*ss0.*ss0);           
B(:,14) = ctt0.*B(:,9) + ss0.*B(:,10)/3;        
B(:,15) = (4/3)*ctt0.*B(:,10);   

B(:,16) = 5*alpha0.*ss0.^4;
B(:,17) = 4*alpha0.*ctt0.*ss0.^3 + beta0.*ss0.^4;
B(:,18) = 3*alpha0.*ctt0.^2.*ss0.^2 + 2*beta0.*ctt0.*ss0.^3;
B(:,19) = 2*alpha0.*ctt0.^3.*ss0 + 3*beta0.*ctt0.^2.*ss0.^2;
B(:,20) = alpha0.*ctt0.^4 + 4*beta0.*ctt0.^3.*ss0;
B(:,21) = 5*beta0.*ctt0.^4;
   

dgg_p = B(:,:)*coeffs(:);


return