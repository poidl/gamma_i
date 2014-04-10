global sns ctns tns pns

clc, dj_tic

dj_disp('extending density surfaces ...'), dj_disp('')

global h_glevels

glevels = eval(h_glevels); ng = length(glevels);

ss_gbdry = []; ctt_gbdry = []; gg_gbdry = []; figure(1)

for level = 1:ng

    [snss,inds_srtd] = sort(sns(level,:));
  
    indss = find(finite(snss(:))); nn = length(indss);
    
    if nn>0

      snss = snss(indss);
      ctnss = ctns(level,inds_srtd); ctnss = ctnss(indss);
      pnss = pns(level,inds_srtd); pnss = pnss(indss);

      for kk = 1:1
        
        deltas = -0.01;
        [sss,cttt] = add_sctp_data(snss(kk),ctnss(kk),pnss(kk),deltas);
        ggg = glevels(level)*ones(size(sss));
        
        ss_gbdry = [ss_gbdry; sss];
        ctt_gbdry = [ctt_gbdry; cttt];
        gg_gbdry = [gg_gbdry; ggg];
        
        hold on; plot(sss,cttt,'m.'), hold off
        dj_pause(1)

        
        deltas = 0.05;
        [sss,cttt] = add_sctp_data(snss(nn-kk+1),ctnss(nn-kk+1),pnss(nn-kk+1),deltas);
        ggg = glevels(level)*ones(size(sss));
        
        ss_gbdry = [ss_gbdry; sss];
        ctt_gbdry = [ctt_gbdry; cttt];
        gg_gbdry = [gg_gbdry; ggg];
        
        hold on; plot(sss,cttt,'m.'), hold off    
        dj_pause(1)
        
      end
    
    end
  
end

save ../an_equation/gamma_bdry ss_gbdry ctt_gbdry gg_gbdry

dj_toc