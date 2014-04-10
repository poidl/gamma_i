function dgg_p = dgamma_p15(ss0,ctt0,alpha0,beta0,coeffs)


nn0 = length(ss0); B = zeros(nn0,15);

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

dgg_p = B(:,:)*coeffs(:);


return