function [ave,percent2,epsilon_x,epsilon_y] = get_gradients(handles)

global s ct t p g ocean n longs lats

global ntpx ntpy nsx nsy logD gmean gstd

global called_ggrads iters dv_orig

global gn limit_panel2 logD2
    
    global xx1 xx2 yy1 yy2 norm

dj_disp('computing Veronis errors ...'), dj_disp('')

dj_tic

figure_no = 0;

if handles.ntp==1
    ntp_surfaces = 1;
else
    ntp_surfaces = 0 ;
end

if handles.ns==1
    ns_surfaces = 1;
    surface_comparison = 2;
else
    ns_surfaces = 0 ;
end

if handles.gfunc==1
    gfunc_surfaces = 1;
    surface_comparison = 3;
else
    gfunc_surfaces = 0 ;
end

if handles.sigp==1
    sigp_surfaces = 1;
    surface_comparison = 4;
else
    sigp_surfaces = 0 ;
end

by_volume = 1; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = zeros(2,1);

% if called_ggrads==0
%     compile_qgrads_t
%     called_ggrads = 1;
% end

[nz,ny,nx] = size(s);

indss = find(isfinite(s));

ss = s(indss); tt = t(indss); pp = p(indss);

s = change_ak(s,'==',nan,-99); t = change_ak(t,'==',nan,-99);
p = change_ak(p,'==',nan,-99); g = change_ak(g,'==',nan,-99); 
ocean = change_ak(ocean,'==',nan,-99); n = change_ak(n,'==',nan,-99);
%longs = longs(:); lats = lats(:); [longs2,lats2] = meshgrid(longs,lats);

inds = find(s==-99); 
    
if by_volume==1
    
    cos_mat = cos(pi*lats/180);
%    cos_mat = cos_mat*ones(1,nx);
    cos_mat = reshape(ones(nz,1)*cos_mat(:)',nz,ny,nx);
    cos_mat(inds) = nan;

    p_mid = (p(1:nz-1,:,:)+p(2:nz,:,:))/2;

    dp(1,:,:) = p_mid(1,:,:)-p(1,:,:); 
    dp(2:nz-1,:,:) = diff(p_mid(:,:,:));
    dp(nz,:,:) = p(nz,:,:)-p_mid(nz-1,:,:);

    i3 = zeros(size(dp));

    inds3 = find(isfinite(dp)); i3(inds3) = 1; 

    zz = i3.*cos_mat.*dp;

%    ocean_volume = nansum(zz(:))*111.2d3*111.2d3
%    inferred_mean_depth = ocean_volume/3.6e14
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp('computing ...')

if ntp_surfaces == 1
%    disp('  neutral tangent plane gradients ...')
%    dj_tic
%    if handles.temp==1
%         options(1) = 1; 
%         [ntpx,ntpy] = quick_gradients(s,t,p,g,ocean,n,longs,lats,options);
%    elseif handles.temp==3
        [ntpx,ntpy] = ntp_gradients(s,ct,p,g,ocean,n,longs,lats);
%    end
%    dj_toc
    ntpx = change_ak(ntpx,'==',-99,nan); ntpy = change_ak(ntpy,'==',-99,nan);
    ntpx = reshape(ntpx,nz,ny,nx); ntpy = reshape(ntpy,nz,ny,nx);
    inds = find(isfinite(ntpx+ntpy));
%    cmd = ['save data/', dataset_str, '/ntp', region_str, temp_str, ' ntpx ntpy'];
%    eval(cmd)
    save ./data/ntp_gradients ntpx ntpy
else
%    cmd = ['load data/', dataset_str, '/ntp', region_str, temp_str];
%    eval(cmd)
    load ./data/ntp_gradients
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ns_surfaces == 1
%    dj_disp('  neutral density gradients ...')
%    dj_tic
%    if handles.temp==1
%         options(1) = 2;
%         inds_s = find(s~=-99); ns = length(inds_s);
%         inds_g = find(g~=-99); ng = length(inds_g); 
%       [nsx,nsy] = quick_gradients(s,t,p,g,ocean,n,longs,lats,options);
%    elseif handles.temp==3
       [nsx,nsy] = ns_gradients(s,ct,p,g,ocean,n,longs,lats);
%    end
%    dj_toc
    nsx = change_ak(nsx,'==',-99,nan); nsy = change_ak(nsy,'==',-99,nan);
    nsx = reshape(nsx,nz,ny,nx); nsy = reshape(nsy,nz,ny,nx);
%    cmd = ['save data/', dataset_str, '/ns', region_str, temp_str, ' nsx nsy'];
%    eval(cmd)
    save ./data/ns_gradients nsx nsy
elseif surface_comparison==2
%    cmd = ['load data/', dataset_str, '/ns', region_str, temp_str];
%    eval(cmd)
    load ./data/ns_gradients
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if gfunc_surfaces == 1
    dj_disp('  gamma_a gradients ...')
%    dj_tic
%    if handles.temp==1
        options(1) = 3;
%         inds_s = find(s~=-99); ns = length(inds_s);
%         inds_g = find(g~=-99); ng = length(inds_g);
%         percent_labeled = 100*ng/ns
        [nsx,nsy] = quick_gradients(s,t,p,g,ocean,n,longs,lats,options);
%    elseif handles.temp==3
%        [nsx,nsy] = ns_gradients(s,ct,p,g,ocean,n,longs,lats);
%    end
%    dj_toc
    nsx = change_ak(nsx,'==',-99,nan); nsy = change_ak(nsy,'==',-99,nan);
    nsx = reshape(nsx,nz,ny,nx); nsy = reshape(nsy,nz,ny,nx);
%    cmd = ['save data/', dataset_str, '/ns', region_str, temp_str, ' nsx nsy'];
%    eval(cmd)
    save data/gammaa_gradients nsx nsy
elseif surface_comparison==3
%    cmd = ['load data/', dataset_str, '/ns', region_str, temp_str];
%    eval(cmd)
    load data/gammaa_gradients
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sigp_surfaces == 1
    disp('  sigma_p gradients ...')
%    dj_tic
%    if handles.temp==1
        options(1) = 4; options(2) = handles.pr0;
        [sigpx,sigpy] = quick_gradients(s,t,p,g,ocean,n,longs,lats,options);
%    elseif handles.temp==3
%        [ntpx,ntpy] = ntp_gradients(s,ct,p,g,ocean,n,longs,lats);
%    end
%    dj_toc
    sigpx = change_ak(sigpx,'==',-99,nan); sigpy = change_ak(sigpy,'==',-99,nan);
    sigpx = reshape(sigpx,nz,ny,nx); sigpy = reshape(sigpy,nz,ny,nx);
    inds = find(isfinite(sigpx+sigpy)); n_sigp = length(inds)
%    cmd = ['save data/', dataset_str, '/sigp', region_str, temp_str, ' sigpx sigpy'];
%    eval(cmd)
    save data/sigp_gradients sigpx sigpy
elseif surface_comparison==4
%    cmd = ['load data/', dataset_str, '/sigp', region_str, temp_str];
%    eval(cmd)
    load data/sigp_gradients
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%							surface statistics computations
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


s = change_ak(s,'==',-99,nan); t = change_ak(t,'==',-99,nan);
p = change_ak(p,'==',-99,nan); g = change_ak(g,'==',-99,nan);
ocean = change_ak(ocean,'==',-99,nan); n = change_ak(n,'==',-99,nan);


%j_inds = find(lats>64); s(:,j_inds,:) = nan;   %%      nan out s for lats > 60


if sigp_surfaces==1
    nsx = sigpx; nsy = sigpy;
end

    
%   percent of ocean labelled

inds_s_good = find(isfinite(s));
inds_g_good = find(isfinite(g));
percent_ocean_labelled = 100*length(inds_g_good)/length(inds_s_good)



indsss = find(isfinite(g+s+nsx+ntpx+nsy+ntpy)); nn = length(indsss);

epsilon_x = nsx-ntpx; epsilon_y = nsy-ntpy; 

if nn>0
    
    figure(1)
    clf
    [hax,hax_noxlabel,hax_noylabel,hax_left,hax_bottom] = nfigaxes([2 2],[0.08 0.08],[0.08 0.92],[0.08 0.92]);
    
    axes(hax(1)) % subplot(2,2,1)
    
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
    
    stats = [p50, p95, percent_ocean_labelled]
    
    if isempty(dv_orig)
        dv_orig(:,1) = x_bins(:);
        dv_orig(:,2) = D_bins(:);
    end
    plot(dv_orig(:,1),dv_orig(:,2),'r'), hold on
    plot(x_bins,D_bins), grid on, hold off
    
    title('log(\it{D_V})'); ylabel('frequency');
    
    Dmax = 1.01*max(D_bins);
    
    %     D_cum = [0,sum(D_bins)]; p95 = interp1(x_bins,D_cum);
    
    set(gcf,'Name','Veronis errors','NumberTitle','off','Color',[0.961 0.988 0.965])
    
    set(gca,'ylim',[0,Dmax])
    
    if by_volume==1
        indss5 = find(logD>=-5); percent1 = 100*length(indss5)/length(indsss);
        i3 = zeros(size(dp)); i3_5 = i3;
        i3(indsss) = 1; i3_5(indss5) = 1;
        zz = i3.*cos_mat.*dp; zz5 = i3_5.*cos_mat.*dp;
        percent2 = 100*nansum(zz5(:))/nansum(zz(:));  N_stats(3) = percent2;
    end
    
    if figure_no==1 | figure_no==2
        xlo = -15; xhi = 0; xrange = xhi-xlo; xx = xlo+0.02*xrange;
        ylo = 0; yhi = 0.375; yrange = yhi-ylo; yy = yhi-0.06*yrange;
        set(gca,'ylim',[ylo yhi])
        text(xx,yy,'(a)')
    elseif figure_no==3 | figure_no==4
        xlo = -15; xhi = 0; xrange = xhi-xlo; xx = xlo+0.02*xrange;
        ylo = 0; yhi = 0.4; yrange = yhi-ylo; yy = yhi-0.06*yrange;
        set(gca,'ylim',[ylo yhi])
        text(xx,yy,'(a)')
    elseif figure_no==7 | figure_no==8
        xlo = -15; xhi = 0; xrange = xhi-xlo; xx = xlo+0.02*xrange;
        ylo = 0; yhi = 0.35; yrange = yhi-ylo; yy = yhi-0.06*yrange;
        set(gca,'ylim',[ylo yhi])
        text(xx,yy,'(a)')
    elseif figure_no==10
        xlo = -15; xhi = 0; xrange = xhi-xlo; xx = xlo+0.02*xrange;
        ylo = 0; yhi = 0.55; yrange = yhi-ylo; yy = yhi-0.06*yrange;
        set(gca,'ylim',[ylo yhi])
        text(xx,yy,'(a)')
    elseif figure_no==14 | figure_no==15
        xlo = -15; xhi = 0; xrange = xhi-xlo; xx = xlo+0.02*xrange;
        ylo = 0; yhi = 0.4; yrange = yhi-ylo; yy = yhi-0.06*yrange;
        set(gca,'ylim',[ylo yhi])
        text(xx,yy,'(a)')
    end
    
    
    %       gamma frequency distribution or gamma residual plot
    %
    %       and gamma frequency distribution stats
    
    axes(hax(2)) % subplot(2,2,2)
    
    inds_g = find(isfinite(g)); gg = g(inds_g);
    
    if mean(gg)>1000, gg = gg-1000; end
    
    [g_bins,x_bins] = hist(gg,100); g_bins = g_bins/sum(g_bins(:));
    
    if sum(size(gn))==0
        whos, plot(x_bins,g_bins), grid on
        xmin = 0.99*min(gg); xmax = 1.01*max(gg); gmax = 1.01*max(g_bins);
        set(gca,'xlim',[xmin,xmax],'ylim',[0,gmax])
    else
        gnn = gn(inds_g); plot(gg,gg-gnn,'.'), grid on
        if sum(size(limit_panel2))~=0, set(gca,'xlim',[27,29],'ylim',[-0.2,0.2]), 
        else set(gca,'ylim',[(min(gg-gnn)-0.0001) (max(gg-gnn)+0.0001)]),end
    end
    %keyboard
    gmean = mean(gg); gstd = std(gg); if gmean>1000, gmean = gmean-1000; end
    title(['Av.{\gamma} ', num2str(gmean,4),',   \sigma ',  num2str(gstd,3)]);
    
    %ylabel('frequency'); %title(zz)
    
    if figure_no==1 | figure_no==2
        xlo = 20; xhi = 30; xrange = xhi-xlo; xx = xlo+0.02*xrange;
        ylo = 0; yhi = 0.25; yrange = yhi-ylo; yy = yhi-0.06*yrange;
        set(gca,'xlim',[xlo xhi],'ylim',[ylo yhi])
        text(xx,yy,'(b)')
    elseif figure_no==3 | figure_no==4
        xlo = 24; xhi = 30; xrange = xhi-xlo; xx = xlo+0.02*xrange;
        ylo = 0; yhi = 0.10; yrange = yhi-ylo; yy = yhi-0.06*yrange;
        set(gca,'xlim',[xlo xhi],'ylim',[ylo yhi])
        text(xx,yy,'(b)')
    elseif figure_no==7 | figure_no==8
        xlo = 20; xhi = 30; xrange = xhi-xlo; xx = xlo+0.02*xrange;
        ylo = 0; yhi = 0.14; yrange = yhi-ylo; yy = yhi-0.06*yrange;
        set(gca,'xlim',[xlo xhi],'ylim',[ylo yhi])
        text(xx,yy,'(b)')
    elseif figure_no==10
        xlo = 22.5; xhi = 29; xrange = xhi-xlo; xx = xlo+0.02*xrange;
        ylo = 0; yhi = 0.075; yrange = yhi-ylo; yy = yhi-0.06*yrange;
        set(gca,'xlim',[xlo xhi],'ylim',[ylo yhi])
        text(xx,yy,'(b)')
    elseif figure_no==14 | figure_no==15
        xlo = 21; xhi = 29; xrange = xhi-xlo; xx = xlo+0.02*xrange;
        ylo = 0; yhi = 0.13; yrange = yhi-ylo; yy = yhi-0.06*yrange;
        set(gca,'xlim',[xlo xhi],'ylim',[ylo yhi])
        text(xx,yy,'(b)')
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%      logD map
    
    logD2 = reshape(nanmax(logD),ny,nx);
    minz = nanmin(logD2(:)); [maxz0,indsm0] = nanmax(logD2(:));
    minz = max(round(minz),-9); maxz = round(maxz0);
    logD2 = change_ak(logD2,'<=',minz,minz);
    logD2 = change_ak(logD2,'>=',maxz,maxz);
    if maxz>maxz0, logD2(indsm0(1)) = maxz; end
    cmd = ['colormap jet(',  int2str(max(1,2*(maxz-minz))), ')']; eval(cmd)
    
    axes(hax(3)) %subplot(2,2,3)
    % whos, dj_pause(0)
    dj_pltmp(nanmean(longs),nanmean(lats'),logD2), grid on
    
    title([' Av.= ',num2str(ave,3),',  50% = ',num2str(p50,3),',  95% = ',num2str(p95,3)])
    
    if figure_no==1 | figure_no==2
        text(-3,20,'(c)')
    elseif figure_no==3 | figure_no==4
        text(-3,-75,'(c)')
    elseif figure_no==7 | figure_no==8
        text(-3,20,'(c)')
    elseif figure_no==10
        xlo = 329.9; xhi = 330.1; xrange = xhi-xlo; xx = xlo+0.02*xrange;
        ylo = -70; yhi = 68; yrange = yhi-ylo; yy = yhi-0.08*yrange;
        text(xx,yy,'(c)')
    elseif figure_no==14 | figure_no==15
        text(-3,20,'(c)')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    axes(hax(4)) % subplot(2,2,4)
    p = reshape(p,size(s)); %pstats = [nanmin(p(:)),nanmax(p(:))]
    pmean = nanmean(p(:,:)');
    k0 = find(pmean>200);
    k1 = change_ak(s,'>=',0,1); k1mean = nansum(k1(k0,:)');
    %    [pmean',[diff(pmean');nan]]
    [nhist,xout] = hist(p(indss5),pmean'); %size(nhist(k0)), size(k1mean),
    k1mean = reshape(k1mean,size(nhist(k0)));
    nhist(k0) = 100*nhist(k0)./k1mean;
    %    barh(xout,nhist)
    barh(xout,nhist,'hist')
    %    plot(nhist,xout)
    
    %    h = findobj(gca,'Type','patch');
    %    set(h,'EdgeColor','w')
    set(gca,'ylim',[xout(1),xout(length(xout))],'ydir','reverse')
    
    perc5 = N_stats(length(N_stats))
    
    if perc5>=10
        title(['5^t^h% = ',num2str(perc5,4),' %']) %, '    ', num2str(p95,3)])
    else
        title(['5^t^h% = ',num2str(perc5,3),' %']) %, '    ', num2str(p95,3)])
    end
    
    grid on

    if figure_no==1
        xlo = 0; xhi = 2.3; xrange = xhi-xlo; xx = xlo+0.9*xrange;
        set(gca,'xlim',[xlo xhi])
        text(xx,700,'(d)')
    elseif figure_no==2
        xlo = 0; xhi = 6.7; xrange = xhi-xlo; xx = xlo+0.9*xrange;
        set(gca,'xlim',[xlo xhi])
        text(xx,350,'(d)')
    elseif figure_no==3
        xlo = 0; xhi = 8; xrange = xhi-xlo; xx = xlo+0.9*xrange;
        set(gca,'xlim',[xlo xhi])
        text(xx,350,'(d)')
    elseif figure_no==4
        xlo = 0; xhi = 55; xrange = xhi-xlo; xx = xlo+0.9*xrange;
        set(gca,'xlim',[xlo xhi])
        text(xx,350,'(d)')
    elseif figure_no==7
        xlo = 0; xhi = 0.3; xrange = xhi-xlo; xx = xlo+0.9*xrange;
        set(gca,'xlim',[xlo xhi])
        text(xx,350,'(d)')
    elseif figure_no==8
        xlo = 0; xhi = 30; xrange = xhi-xlo; xx = xlo+0.9*xrange;
        set(gca,'xlim',[xlo xhi])
        text(xx,350,'(d)')
    elseif figure_no==10
        xlo = 0; xhi = 1.2; xrange = xhi-xlo; xx = xlo+0.9*xrange;
        set(gca,'xlim',[xlo xhi])
        text(xx,350,'(d)')
    elseif figure_no==14
        xlo = 0; xhi = 1.2; xrange = xhi-xlo; xx = xlo+0.9*xrange;
        set(gca,'xlim',[xlo xhi])
        text(xx,600,'(d)')
    elseif figure_no==15
        xlo = 0; xhi = 2; xrange = xhi-xlo; xx = xlo+0.9*xrange;
        set(gca,'xlim',[xlo xhi])
        text(xx,600,'(d)')
    end
    
else
    
    whos
    
    inds = find(isfinite(g)); ng = length(inds)
    inds = find(isfinite(s)); ns = length(inds)
    inds = find(isfinite(ntpx)); nntpx = length(inds)
    inds = find(isfinite(nsx)); nnsx = length(inds)
    inds = find(isfinite(ntpy)); nntpy = length(inds)
    inds = find(isfinite(nsy)); nnsy = length(inds)
    
    error('*****        no data')
    
end

dj_toc

% inds_poss = find(s>=0&p>=200); percent_checked = 100*length(indsss)/length(inds_poss)
% 
 dj_pause(1)
    

return
