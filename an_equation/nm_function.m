function z = nm_function(x)

global cmin cmin1 cmin2 cmin3 cmin4 pmax 

global cost_plot x_plot cost_plot1 cost_plot2 cost_plot3 cost_plot4 cost3_model

global s ct t p g g_rf longs lats ocean n

global f_called ss_gbdry ctt_gbdry gg_gbdry

global indss ss tt ctt pp gg gg_fit

global ss0 ctt0 alpha0 beta0

global ss_overlap ctt_overlap gg_overlap

global gfunc_coeffs xden gmodel npoly

global h_numerator h_denominator h_normalise h_boundary ave

wt1 = 1;  wt2 = 1; wt3 = 10; wt4 = 10^6; f_called = f_called+1;


if eval(h_numerator)==7 & eval(h_denominator)==9
    xx = [x; xden];
    cmd = ['g_fit = gamma_prf_', h_numerator, '_', h_denominator, '(s,ct,g,ocean,lats,xx);']
    eval(cmd)
    gg_fit = g_fit(indss);
else
    xx = x;
    if strcmp(gmodel,'gamma_p')
        cmd = ['gg_fit = gamma_p', h_numerator, '(ss,ctt,xx);']
        eval(cmd)
    elseif strcmp(gmodel,'gamma_pp')
        cmd = ['g_fit = gamma_pp', h_numerator, '(s,ct,g,ocean,lats,xx);']
        eval(cmd)
        gg_fit = g_fit(indss);
    end
end

gfunc_coeffs = xx;

rms = sqrt(mean((gg-gg_fit).*(gg-gg_fit)));

cost1 = rms

% cmd = ['cost1a = rfunc_', h_numerator, '_', h_denominator, '_ss(ss_overlap,ctt_overlap,gg_overlap,xx);']
% 
% eval(cmd)
% 
% cost1 = cost1 + 10*cost1a;

% if h_normalise==1, cost1 = cost1*30; end

if eval(h_numerator)==7 & eval(h_denominator)==9
    cmd = ['save gamma_prf_', h_numerator, '_', h_denominator, '_23.dat xx -ascii -double'];
    eval(cmd)
else
    if strcmp(gmodel,'gamma_p')
        cmd = ['save gamma_p', h_numerator, '.dat xx -ascii -double'];
    elseif strcmp(gmodel,'gamma_pp')
        cmd = ['save gamma_pp', h_numerator, '.dat xx -ascii -double'];
    end    
    eval(cmd)
end



g(indss) = gg_fit; 

 cd errors
         get_grads_nm
 cd ..

cost2 = N_stats(1)

cost3 = N_stats(3)

if strcmp(gmodel,'gamma_p')
    cmd = ['dgg_p = dgamma_p', h_numerator, '(ss0,ctt0,alpha0,beta0,xx);']
elseif strcmp(gmodel,'gamma_pp')
    cmd = ['dgg_p = dgamma_pp', h_numerator, '(ss0,ctt0,alpha0,beta0,xx);']
end
eval(cmd)

cost4 = sqrt(mean(dgg_p.*dgg_p));


costs = [cost1,cost2,cost3,cost4,wt1*cost1+wt2*cost2+wt3*cost3+wt4*cost4]


function_evaluation = [f_called,costs];

%dj_toc, %dj_pause(0)

z = costs(length(costs));

z1 = costs(1); z2 = costs(2); z3 = costs(3); z4 = costs(4);


%%      plot penalty function

nplot = 51; 
if f_called==1
    cmin = z; cmin1 = z1; cmin2 = z2; cmin3 = z3; cmin4 = z4;
    cost_plot = nan*ones(nplot,1); x_plot = nan*ones(nplot,1);
    cost_plot1 = nan*ones(nplot,1);
    cost_plot2 = nan*ones(nplot,1);
    cost_plot3 = nan*ones(nplot,1);
    cost_plot4 = nan*ones(nplot,1);
    cost_plot(f_called) = z; x_plot(f_called) = f_called;
    cost_plot1(f_called) = z1; 
    cost_plot2(f_called) = z2; 
    cost_plot3(f_called) = z3; 
    cost_plot4(f_called) = z4; 
elseif f_called<=nplot
    cost_plot(f_called) = z; x_plot(f_called) = f_called;
    cost_plot1(f_called) = z1; 
    cost_plot2(f_called) = z2; 
    cost_plot3(f_called) = z3; 
    cost_plot4(f_called) = z4; 
else
    cost_plot = [cost_plot(2:nplot); z]; x_plot = [x_plot(2:nplot); f_called];
    cost_plot1 = [cost_plot1(2:nplot); z1]; 
    cost_plot2 = [cost_plot2(2:nplot); z2]; 
    cost_plot3 = [cost_plot3(2:nplot); z3]; 
    cost_plot4 = [cost_plot4(2:nplot); z4]; 
end

cmin = min([cmin;z]); cmin1 = min([cmin1;z1]); cmin2 = min([cmin2;z2]); cmin3 = min([cmin3;z3]); cmin4 = min([cmin4;z4]);

if mod(f_called,10)==0
    
    figure(1)
    
    subplot(3,2,1)
        plot(x_plot,cost_plot1), grid on, hold on
        plot(x_plot,cmin1*ones(size(x_plot)),'r'), hold off
        title([num2str(wt1), '    ', num2str(wt1*z1,3), '    ', num2str(z1,3)])

    subplot(3,2,2)
        plot(x_plot,cost_plot2), grid on, hold on
        plot(x_plot,cmin2*ones(size(x_plot)),'r'), hold off
        title([num2str(wt2), '    ', num2str(wt2*z2,3), '    ', num2str(z2,3)])
        
    subplot(3,2,3)
        plot(x_plot,cost_plot3), grid on, hold on
        plot(x_plot,cmin3*ones(size(x_plot)),'r'), hold off
        title([num2str(wt3), '    ', num2str(wt3*z3,3), '    ', num2str(z3,3)])
        
    subplot(3,2,4)
        plot(x_plot,cost_plot4), grid on, hold on
        plot(x_plot,cmin4*ones(size(x_plot)),'r'), hold off
        title([num2str(wt4), '    ', num2str(wt4*z4,3), '    ', num2str(z4,3)])
        
    subplot(3,2,6)
        plot(x_plot,cost_plot), grid on, hold on
        plot(x_plot,cmin*ones(size(x_plot)),'r'), hold off
        title(num2str(z,3))
        
    figure(gcf), dj_pause(1)
    
end
    
return