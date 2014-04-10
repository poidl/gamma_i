function rms = rfunc_15_0_ss(coeffs)

global f_called indss ss ctt pp gg gg_rf 

global ss0 ctt0 alpha0 beta0 h_normalise

nn = length(ss); A = zeros(nn,15); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
A(:,1)  = ones(nn,1);   
A(:,2)  = ss;
A(:,3)  = ctt ;
    
A(:,4)  = ss.*ss;                %   s^2
A(:,5)  = ss.*ctt ;                %   s   t
A(:,6)  = ctt .*ctt ;                %       t^2
    
A(:,7)  = ss.*A(:,4);            %   s^3
A(:,8)  = ss.*A(:,5);            %   s^2 t
A(:,9)  = ss.*A(:,6);            %   s   t^2
A(:,10) = ctt .*A(:,6);            %       t^3
    
A(:,11) = ss.*A(:,7);            %   s^4 
A(:,12) = ss.*A(:,8);            %   s^3 t
A(:,13) = ss.*A(:,9);            %   s^2 t^2
A(:,14) = ss.*A(:,10);           %   s   t^3
A(:,15) = ctt .*A(:,10);           %       t^4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gg_rf = A(:,:)*coeffs(:);

rms_g = sqrt(mean((gg-gg_rf).*(gg-gg_rf)));

dgg_rf = B(:,:)*coeffs(:);

rms_dg = sqrt(mean(dgg_rf.*dgg_rf));

if h_normalise==1, rms_g = 30*rms_g; end

rms_dg = 20*rms_dg; 

rms = (rms_g + rms_dg); 

rmss = [rms_g rms_dg rms]

%[ss(1),ctt(1),gg(1),gg_rf(1)], 	dj_pause(0)



return