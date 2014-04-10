global A C b


inversion_method = 7; 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dj_disp('inside fitl_rfunc  -  assembling matrix ...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(handles.denominator,'0')                                          %   a polynomial
    
    indss = find(finite(s+ct+g));
    
    ss = s(indss); ctt = ct(indss); gg = g(indss); wt = 1;
    
    if handles.boundary==1
        
        ss_b1 = 0:3:33; ctt_b1 = -5:3:40;                %   [0,33]x[-5,40]
        [ss_b2,ctt_b2] = meshgrid(ss_b1,ctt_b1);
        ss0 = ss_b2(:); ctt0 = ctt_b2(:); 
%         
%         ss_b1 = 0:3:42; ctt_b1 = 27:3:40;                %   [0,42]x[27,40]
%         [ss_b2,ctt_b2] = meshgrid(ss_b1,ctt_b1);
%         ss0 = [ss0; ss_b2(:)]; ctt0 = [ctt0; ctt_b2(:)]; 
%         
%         ss_b1 = 38:1:42; ctt_b1 = -5:3:40;                %   [38,42]x[-5,40]
%         [ss_b2,ctt_b2] = meshgrid(ss_b1,ctt_b1);
%         ss0 = [ss0; ss_b2(:)]; ctt0 = [ctt0; ctt_b2(:)]; 
        
        pp0 = zeros(size(ss0));       
        [rho,rhoss,rhoctt,rhopp] = eosall_from_ct(ss0,ctt0,pp0);           
        alpha0 = -rhoctt./rho; beta0 = rhoss./rho;
        
        load gamma_bdry
         
        inds_gb = find(finite(ctt_gbdry));
        no_boundary_gammas = length(inds_gb)

        load sig0_bdry

        no_boundary_sig0s = length(ss)

        if no_boundary_sig0s>0
            ss_gbdry = [ss_gbdry; ss_sigbdry];
            ctt_gbdry = [ctt_gbdry; ctt_sigbdry];
            gg_gbdry = [gg_gbdry; gg_sigbdry];
        end    
            
        
    end

    if handles.normalise==1
        ss = ss/40; ctt = ctt/30; gg = gg/30; wt = 30;   %mean(ss),mean(ctt),mean(gg)
        if handles.boundary==1
            ss0 = ss0(:)/40; ctt0 = ctt0(:)/30;   
            alpha0 = alpha0(:)*0.75; beta0 = beta0(:);
            
            ss_gbdry = ss_gbdry(:)/40; ctt_gbdry = ctt_gbdry(:)/30; gg_gbdry = gg_gbdry(:)/30;  

        end
    end
    
    cmd = ['A = Agen', handles.numerator, '(ss,ctt);'], eval(cmd);
    
    b = gg;
    
    awt = 5;
    
    A = awt*A; b = awt*b;
    
    [nA,m] = size(A)
    
    if handles.boundary==1
        
%         cmd = ['B = Bgen', handles.numerator, '(ss0,ctt0,alpha0,beta0);'], eval(cmd);
%         
%         B = 1*B;
%         
%         [nB,m] = size(B)
%         
%         b0 = zeros(nB,1);
%         
%         A = [A;B]; b = [b(:);b0(:)];
        
        cmd = ['C = Agen', handles.numerator, '(ss_gbdry,ctt_gbdry);'], eval(cmd);
        
    	c = gg_gbdry(:);
        
        [nC,m] = size(C)
        
        cwt = 1;
        
        C = cwt*C; c = cwt*c;
        
        A = [A;C]; b = [b; c];
        
    end
    
else                                                                        %   a rational function
    
    cmd = ['A = Agen', handles.numerator, '_', handles.denominator, '(s,ct,g);'], eval(cmd);
    b = g/30;                                    
    
end


[neqs,nvars] = size(A)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dj_disp('inverting ...'), figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%   solution 1

x1= A\b;

if strcmp(handles.denominator,'0')
    b1 = A(1:nA,:)*x1(:); resids1 = wt*(b(1:nA)-b1);       %   resids are in gamma units
    rms1 = sqrt(mean(resids1.*resids1)); max1 = max(abs(resids1));
    subplot(2,2,1), hist(resids1,200), grid on, title('residuals 1')
else
    p1 = 1e8*A(:,:)*x1(:); resids1 = 1e-4*(p-p1);       %   pressure comparison in dbars
    rms1 = sqrt(mean(resids1.*resids1)); max1 = max(abs(resids1));
    subplot(2,2,1), hist(resids1,200), grid on, title('residuals 1')
end

%   solution 2

[Q,R] = qr(A,0); y = Q'*b; x2 = R\y; 

if strcmp(handles.denominator,'0')
    b2 = A(1:nA,:)*x2(:); resids2 = wt*(b(1:nA)-b2);       %   resids are in gamma units
    rms2 = sqrt(mean(resids2.*resids2)); max2 = max(abs(resids2));
    subplot(2,2,2), hist(resids2,200), grid on, title('residuals 2')
else
    p2 = 1e8*A(:,:)*x2(:); resids2 = 1e-4*(p-p2);
    rms2 = sqrt(mean(resids2.*resids2)); max2 = max(abs(resids2));
    subplot(2,2,2), hist(resids2,200), grid on, title('residuals 2')
end


%   solution 3

% x3 = R\(R'\(A'*b1)); r = b1-A*x3; e = R\(R'\(A'*r)); x3 = x3+e;

x3 = pinv(A)*b;

if strcmp(handles.denominator,'0')
    b3 = A(1:nA,:)*x3(:); resids3 = wt*(b(1:nA)-b3);       %   resids are in gamma units
    rms3 = sqrt(mean(resids3.*resids3)); max3 = max(abs(resids3));
    subplot(2,2,3), hist(resids3,200), grid on, title('residuals 3')
else
    p3 = 1e8*A(:,:)*x3(:); resids3 = 1e-4*(p-p3);
    rms3 = sqrt(mean(resids3.*resids3)); max3 = max(abs(resids3));
    subplot(2,2,3), hist(resids3,200), grid on, title('residuals 3'), figure(gcf)
end

 

solutions =  [x1(:),x2(:),x3(:)]

rms = [rms1,rms2,rms3], max_abs = [max1,max2,max3]

x0 = x3;

if eval(handles.denominator)~=0
    cmd = ['save rfunc_', handles.numerator, '_', handles.denominator, '.dat x0 -ascii -double']
else
    cmd = ['save gamma_p', handles.numerator, '.dat x0 -ascii -double']
end

eval(cmd)

