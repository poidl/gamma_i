function gg_p = gamma_p15(ss,ctt,coeffs)


nn = length(ss); A = zeros(nn,15); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
A(:,1)  = ones(nn,1);   
A(:,2)  = ss;
A(:,3)  = ctt ;
    
A(:,4)  = ss.*ss;                %   s^2
A(:,5)  = ss.*ctt ;              %   s   t
A(:,6)  = ctt .*ctt ;            %       t^2
    
A(:,7)  = ss.*A(:,4);            %   s^3
A(:,8)  = ss.*A(:,5);            %   s^2 t
A(:,9)  = ss.*A(:,6);            %   s   t^2
A(:,10) = ctt .*A(:,6);          %       t^3
    
A(:,11) = ss.*A(:,7);            %   s^4 
A(:,12) = ss.*A(:,8);            %   s^3 t
A(:,13) = ss.*A(:,9);            %   s^2 t^2
A(:,14) = ss.*A(:,10);           %   s   t^3
A(:,15) = ctt .*A(:,10);         %       t^4


gg_p = A(:,:)*coeffs(:);


return