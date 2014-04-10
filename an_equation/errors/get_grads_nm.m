global cos_mat dp my_computer dir0 fig_title ntpx ntpy

global s ct t p g g_rf longs lats ocean n

global indss ss tt ctt pp gg gg_fit

global f_called gfunc_coeffs

global ave h_numerator h_denominator h_normalise

global N_stats0


%%      initialise

by_volume = 1; options = zeros(2,1); 

[nz,ny,nx] = size(s);

if f_called==1
    
    dir0 = 'd:/neutrals/ness8/an_equation/data';
    
	compile_qgrads_nm

	s = change(s,'==',nan,-99); g = change(g,'==',nan,-99); 
  
	ocean = change(ocean,'==',nan,-99);
  
	inds = find(s==-99);
  
	if by_volume==1
        cos_mat = cos(pi*lats/180);
        cos_mat = reshape(ones(nz,1)*cos_mat(:)',nz,ny,nx);
        cos_mat(inds) = nan;
        p_mid = (p(1:nz-1,:,:)+p(2:nz,:,:))/2;
        dp(1,:,:) = p_mid(1,:,:)-p(1,:,:); 
        dp(2:nz-1,:,:) = diff(p_mid(:,:,:));
        dp(nz,:,:) = p(nz,:,:)-p_mid(nz-1,:,:);
	end 

  
%%		local tangent plane gradients

    options(1) = 1;
    
    disp('  tangent plane gradients')
    
    which quick_gradients
    
    [ntpx,ntpy] = quick_gradients(s,t,p,g,ocean,n,longs,lats,options);
    
    ntp_data = [ntpx(:),ntpy(:)];
    ntpx = change(ntpx,'==',-99,nan); ntpy = change(ntpy,'==',-99,nan);
    ntpx = reshape(ntpx,nz,ny,nx); ntpy = reshape(ntpy,nz,ny,nx);
    inds = find(finite(ntpx+ntpy)); n_ntp = length(inds);
   
end


%%		gamma function gradients

options(1) = 2; 

	s = change(s,'==',nan,-99); g = change(g,'==',nan,-99); 
  
	ocean = change(ocean,'==',nan,-99);
    
   
[gfuncx,gfuncy] = quick_gradients(s,ct,p,g,ocean,n,longs,lats,options);

	s = change(s,'==',-99,nan); g = change(g,'==',-99,nan); 
  
	ocean = change(ocean,'==',-99,nan);

gfuncx = change(gfuncx,'<=',-99,nan); gfuncy = change(gfuncy,'<=',-99,nan);

gfuncx = reshape(gfuncx,nz,ny,nx); gfuncy = reshape(gfuncy,nz,ny,nx);
    

%%      surface statistic computations


nsx = gfuncx; nsy = gfuncy;


indsss = find(finite(g+s+nsx+ntpx+nsy+ntpy));  n_poss = length(indsss)

if ~isempty(indsss)
    
  xx1 = nsx(indsss); yy1 = ntpx(indsss); xx2 = nsy(indsss); yy2 = ntpy(indsss);
  norm = (xx1-yy1).*(xx1-yy1)+(xx2-yy2).*(xx2-yy2);

  D = 1000*norm; logD0 = log10(D);
  
  logD_srtd = sort(logD0); nD = length(logD_srtd);            %   95th percentile
    
  p50 = interp1((1:nD)/nD, logD_srtd, 0.50);
  
  p95 = interp1((1:nD)/nD, logD_srtd, 0.95);
    
  stats = [p50, p95]
    
  logD = nan*ones(size(s)); logD(indsss) = logD0;
    
  [D_bins,x_bins] = hist(logD0,[-15:0.1:0]);
  nD = sum(D_bins); D_integral = dj_quad(x_bins,D_bins); D_bins = D_bins/D_integral;

  percent = 100*dj_quadi(x_bins,D_bins,-5,0); 

  ave = dj_quad(x_bins,x_bins.*D_bins); % N_stats = [ave,p95,percent]
    
  inds_poss = find(g>=0&p>=200); percent_checked = 100*length(indsss)/length(inds_poss);
  inds_poss1 = find(s>=0&p>=200); percent_checked1 = 100*length(indsss)/length(inds_poss1)
     
else
  N_stats = [0,100]; percent_checked = 0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


indsss = find(finite(g+s+nsx+ntpx+nsy+ntpy)); nn = length(indsss);

epsilon_x = nsx-ntpx; epsilon_y = nsy-ntpy; 

if nn>0
  
  figure(2)
  
  set(gcf,'Name','Veronis errors','NumberTitle','off','Color',[0.961 0.988 0.965])
    
  subplot(2,2,1)
  
	xx1 = nsx(indsss); yy1 = ntpx(indsss); xx2 = nsy(indsss); yy2 = ntpy(indsss);
    
	norm = (xx1-yy1).*(xx1-yy1)+(xx2-yy2).*(xx2-yy2);
    
    sx_diff = ntpx(indsss)-nsx(indsss); sy_diff = ntpy(indsss)-nsy(indsss);

	D = 1000*norm; logD0 = log10(D);
   
    logD = nan*ones(nz,ny,nx);

    logD(indsss) = logD0; 

	[D_bins,x_bins] = hist(logD0,[-15:0.1:0]);
    
	nD = sum(D_bins); D_integral = dj_quad(x_bins,D_bins); D_bins = D_bins/D_integral;

	percent = 100*dj_quadi(x_bins,D_bins,-5,0);

	ave = dj_quad(x_bins,x_bins.*D_bins); N_stats = [ave,percent];
 
    logD_srtd = sort(logD0); nD = length(logD_srtd);            %   95th percentile
    
    p50 = interp1((1:nD)/nD, logD_srtd, 0.50);
  
    p95 = interp1((1:nD)/nD, logD_srtd, 0.95);
    
    stats = [p50, p95]
  
  
	plot(x_bins,D_bins), grid on

	title('log(\it{D_V})'); ylabel('frequency');
    
    Dmax = 1.01*max(D_bins); 
    
%     D_cum = [0,sum(D_bins)]; p95 = interp1(x_bins,D_cum);
    
    set(gca,'ylim',[0,Dmax])

    if by_volume==1
        indss5 = find(logD>=-5); percent1 = 100*length(indss5)/length(indsss);
        i3 = zeros(size(dp)); i3_5 = i3; 
        i3(indsss) = 1; i3_5(indss5) = 1;
        zz = i3.*cos_mat.*dp; zz5 = i3_5.*cos_mat.*dp;
        percent2 = 100*nansum(zz5(:))/nansum(zz(:));  N_stats(3) = percent2;     
    end

    
%%       gamma frequency distributions

  subplot(2,2,2)
  
    ggg = gg;
  
    gmean0 = mean(ggg), gstd0 = std(ggg);
    
    if gmean0>1000
        gmean0 = gmean0-1000; ggg = ggg-1000;
    end
    
	[g_bins,x_bins] = hist(ggg,100); g_bins = g_bins/sum(g_bins(:));
      
	plot(x_bins,g_bins), grid on
    
    xmin = 0.99*min(ggg); xmax = 1.01*max(ggg); gmax = 1.01*max(g_bins);

    set(gca,'xlim',[xmin,xmax],'ylim',[0,gmax])

	title(['\gamma     ', num2str(gmean0,4),'     ',  num2str(gstd0,3)]);
    
    
	ggg = gg_fit;
    
    if nanmean(ggg)>1000, ggg = ggg-1000; end
     
    gmean = mean(ggg); gstd = std(ggg); 
    
	[g_bins,x_bins] = hist(ggg,100); g_bins = g_bins/sum(g_bins(:));
      
	hold on, plot(x_bins,g_bins,'r'), hold off
    
    xmin = min([xmin; 0.99*min(ggg)]); xmax = max([xmax; 1.01*max(ggg)]); gmax = max([gmax; 1.01*max(g_bins)]);

    set(gca,'xlim',[xmin,xmax],'ylim',[0,gmax])

	xlabel(['\gamma     ', num2str(gmean,4),'     ',  num2str(gstd,3)]);
    
    
%%      logD map

    logD2 = reshape(nanmax(logD),ny,nx);
    minz = nanmin(logD2(:)); [maxz0,indsm0] = nanmax(logD2(:));
    minz = max(round(minz),-9); maxz = round(maxz0);
    logD2 = change(logD2,'<=',minz,minz);
    logD2 = change(logD2,'>=',maxz,maxz);
    if maxz>maxz0, logD2(indsm0(1)) = maxz; end
    cmd = ['colormap jet(',  int2str(max(1,2*(maxz-minz))), ')']; eval(cmd)

  subplot(2,2,3)
  
    dj_pltmp(nanmean(longs),nanmean(lats'),logD2), grid on
    
    title([num2str(ave,3), '    ', num2str(p50,3), '    ', num2str(p95,3), '    ', num2str(percent2,3), '%'])

    
        
%   percent of ocean labelled

inds_ntp_good = find(finite(ntpx+ntpy));
inds_ns_good = find(finite(nsx+nsy));
percent_ocean_checked = 100*length(inds_ns_good)/length(inds_ntp_good)

N_stats(4) = percent_ocean_checked;


%%      contour map
    
  subplot(2,2,4)
   
    smin = 25; smax = 42; sby = (smax-smin)/100;

    ctmin = -5; ctmax = 40; ctby = (ctmax-ctmin)/100;

    s0 = smin:sby:smax; ct0 = ctmin:ctby:ctmax; ns = length(s0); nct = length(ct0);

    [sss,cttt] = meshgrid(s0,ct0); sss = sss(:); cttt = cttt(:);

    if h_normalise==1
        sss = sss/40; cttt = cttt/30;
    end


    cd ..
      if eval(h_denominator)==0
        coeffs = gfunc_coeffs(1:eval(h_numerator));
        cmd = ['g0 = gamma_p', h_numerator, '(sss,cttt,coeffs);']; eval(cmd)
      end
    cd errors

    g0 = reshape(g0,nct,ns);

    if h_normalise==1
    	g0 = 30*g0;
    end

    if nanmean(g0)>1000, g0 = g0-1000; end

    fpcolor(s0,ct0,g0), colorbar
    inds = find(finite(s(:))); sss = s(inds); cttt = ct(inds);
    hold on, plot(sss,cttt,'w.')
    contour(s0,ct0,g0,20,'k'), hold off
    
else

    inds = find(finite(g));     ng = length(inds)
    inds = find(finite(s));     ns = length(inds)
    inds = find(finite(ntpx));  nntpx = length(inds)
    inds = find(finite(nsx));   nnsx = length(inds)
    inds = find(finite(ntpy));  nntpy = length(inds)
    inds = find(finite(nsy));   nnsy = length(inds)

    disp('*****        no data')
    
    N_stats = 1.05*N_stats0;
    
end


N_stats0 = N_stats;
dj_toc
nanmean(g(:))
if f_called==1
    dj_pause(0)
else
    dj_pause(1)
end
    