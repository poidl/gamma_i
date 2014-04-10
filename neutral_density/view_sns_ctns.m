global sns ctns pns longs lats g

%%      latitudinal spectrum

cmap = colormap(jet(64)); nloop = 0; longs1 = longs; lats1 = lats;

inds = find(~finite(g(1,:))); longs1(inds) = nan; lats1(inds) = nan;

longss = nanmean(longs); latss = nanmean(lats');

[ng,ny,nx] = size(sns+ctns); figure, dj_tic

for j = 1:ny   
for i = 1:nx

  ss = sns(:,j,i); ctt = ctns(:,j,i); 
    
  inds = find(finite(ss));    % [k_inds,i_inds] = ind2sub(size(ss),inds);

  if length(inds)>0
       
      nloop = nloop+1; 
      
      ind_map = 1+round(((longs(j,i)-longss(1))/(longss(nx)-longss(1)))*63);
         
      subplot(2,2,1), plot(ss,ctt,'.','color',cmap(ind_map,:))

      if nloop==1, grid on, hold on, end
      
      
      ind_map = 1+round(((lats(j,i)-latss(1))/(latss(ny)-latss(1)))*63);
    
      subplot(2,2,2), plot(ss,ctt,'.','color',cmap(ind_map,:))     
      
      if nloop==1, grid on, hold on, end
          
      if nloop==1
          subplot(2,2,3), dj_pltmp(longss,latss,longs1), grid on
          subplot(2,2,4), dj_pltmp(longss,latss,lats1), grid on
      end

      
    end
    
  end
    
  if mod(j,1)==0 | j==ny
        set(gca,'zdir','reverse')
        figure(gcf), dj_pause(1)
  end
    
end

hold off

dj_toc

return


inds = find(finite(g)&p>plo&p<phi); nn = length(inds);

ss = s(inds); tt = t(inds); pp = p(inds);

gg = g(inds); gmin = min(gg), gmax = max(gg)

for kk = 1:nn
    
    ind_map = 1+round(((gg(kk)-gmin)/(gmax-gmin))*63); 

    figure(3), set(gcf,'Position',[709   446   706   489])
   
    plot3(ss(kk),tt(kk),pp(kk),'.','color',cmap(ind_map,:))
    
    if kk==1, grid on, hold on, end
    
    if mod(kk,10000)==0 | kk==nn
        set(gca,'zdir','reverse')  
        figure(gcf), dj_pause(1)
    end
  
end


