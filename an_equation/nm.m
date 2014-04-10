global f_called cmin cmin1 cmin2 cmax3 cmin4 pmax 

global cost_plot x_plot cost_plot1 cost_plot2 cost_plot3  

global s ct t p g longs lats ocean n

global f_called indss ss ctt gg ss_gbdry ctt_gbdry gg_gbdry

global ss0 ctt0 alpha0 beta0

global ss_overlap ctt_overlap gg_overlap

global h_numerator h_denominator h_normalise h_boundary ave

global xden gmodel npoly


disp(' '), ok = 'inside nm ...'

by_volume = 1; gmodel = 'gamma_pp'; npoly = 21;

wts_lim = 1;

path0 = path;

indss = find(finite(s+ct+p+g)); 

ss = s(indss); ctt = ct(indss); pp = p(indss); gg = g(indss);

[nz,ny,nx] = size(g); lats3 = reshape(ones(nz,1)*lats(:)',nz,ny,nx);


h_numerator = handles.numerator; h_denominator = handles.denominator; h_normalise = handles.normalise;

h_boundary = handles.boundary;

if strcmp(handles.numerator,'7') & strcmp(handles.denominator,'9')
    gg = gg+1000; 
end


%       boundary data

if handles.boundary==1        
        
        ss_b1 = 0:0.1:30; ctt_b1 = -5:0.1:40;                 %   [0,30]x[-5,40]
        [ss_b2,ctt_b2] = meshgrid(ss_b1,ctt_b1);
        ss0 = ss_b2(:); ctt0 = ctt_b2(:); 
        
        ss_b1 = [30:0.1:42]; ctt_b1 = 25:0.1:40;              %   [30,42]x[25,40]
        [ss_b2,ctt_b2] = meshgrid(ss_b1,ctt_b1);
        ss0 = [ss0; ss_b2(:)]; ctt0 = [ctt0; ctt_b2(:)]; 
     
        pp0 = zeros(size(ss0));       
        [rho,rhoss,rhoctt,rhopp] = eosall_from_ct(ss0,ctt0,pp0);           
        alpha0 = -rhoctt./rho; beta0 = rhoss./rho;
        
        boundary_condtions_neqs = length(ss0)

end

%         
%         load gamma_bdry
%          
%         inds_gb = find(finite(ctt_gbdry));
%         no_boundary_gammas = length(inds_gb)
% 
%         load sig0_bdry
% 
%         no_boundary_sig0s = length(ss_sigbdry)
% 
%         if no_boundary_sig0s>0
%             ss_gbdry = [ss_gbdry; ss_sigbdry];
%             ctt_gbdry = [ctt_gbdry; ctt_sigbdry];
%             gg_gbdry = [gg_gbdry; gg_sigbdry];
%             if strcmp(handles.numerator,'7') & strcmp(handles.denominator,'9')
%                 gg_gbdry = gg_gbdry+1000;
%             end
%         end
%         
% end
% 
% if handles.normalise==1
%         ss = ss/40; ctt = ctt/30; gg = gg/30; wt = 30;   
%         if handles.boundary==1
%             ss0 = ss0(:)/40; ctt0 = ctt0(:)/30;   
%             alpha0 = alpha0(:)*0.75; beta0 = beta0(:);
%         end
% end

f_called = 0

% cmd = ['copyfile(''rfunc_', handles.numerator, '_', handles.denominator, 'a.dat'', ''input_data.dat'')']
% 
% eval(cmd)
% 
% load input_data.dat
% 
% xa = input_data(:); 
% 
% cmd = ['copyfile(''rfunc_', handles.numerator, '_', handles.denominator, 'b.dat'', ''input_data.dat'')']
% 
% eval(cmd)
% 
% load input_data.dat
% 
% xb = input_data(:); 

if eval(h_numerator)==7 & eval(h_denominator)==9
    load gamma_prf_7_9_23.dat
    z = gamma_prf_7_9_23;
    xa = [z(1:7); z(15:23)]; xb = z(8:14);    x14 = [xa(1:7);xb(1:7)], xden = xa(8:16)
    epsx = 0.01; x14 = (1+epsx*(2*rand(size(x14))-1)).*x14;
else
    if strcmp(gmodel,'gamma_p')
        cmd = ['load gamma_p', handles.numerator, '.dat']; eval(cmd)    
        cmd = ['x = gamma_p', handles.numerator], eval(cmd)
    elseif strcmp(gmodel,'gamma_pp')
%         x = generate_coefficients(npoly)
        load gamma_pp21.dat
        x = gamma_pp21
    end
    epsx = 0.0; x = (1+epsx*(2*rand(size(x))-1)).*x;
end

options = optimset('Display', 'iter', 'TolX', 1e-4, 'TolFun', 1e-3, 'MaxFunEvals', 100000, 'MaxIter', 1000000);

if eval(h_numerator)==7 & eval(h_denominator)==9
    [y14,fval] = fminsearch(@nm_function,x14,options)
%     y = [y14; xden];
else
    [y,fval] = fminsearch(@nm_function,x,options)
end


return
