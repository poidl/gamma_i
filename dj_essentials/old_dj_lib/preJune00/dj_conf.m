function dj_conf(x,y,z)

%%%    dj_conf           modified filled contour function
%%%
%%%    Usage:       dj_conf(x,y,z)
%%%
%%%    Author:      David Jackett
%%%
%%%    Date:        29/4/98
%%%


colormap('jet');

discrete = 1; cc = 12:2:30; 

discrete = 1; cc = -1:0.5:5; 

n = length(get(gcf,'CurrentAxes')); if n>0, figure; end

cmin = nanmin(nanmin(z)), cmax = nanmax(nanmax(z))

%[C,h] = contour(x,y,z,'k-'); hold on;

[C,h1,CF] = contourf(x,y,z,[cmin,cc,cmax]);

caxis([cmin,cmax]); set(gca,'DataAspectRatio',[1,1,1])

if discrete==1
	cmap = colormap;
	inds = round(1+63*(cc-cmin)/(cmax-cmin)); inds = [1,inds];
	n1 = 256;
	inds1 = floor(1+(n1-1)*(cc-cmin)/(cmax-cmin)); inds1 = [1,inds1,n1+1];
	nn = length(inds1)-1; cmap1 = zeros(n1,3);
	for ii = 1:nn
		cmap1(inds1(ii):inds1(ii+1)-1,:) = ones(inds1(ii+1)-inds1(ii),1)*cmap(inds(ii),:);
	end
	colormap(cmap1)
	h2 = colorbar('horizontal');
	set(h2,'TickLength',[0,0])
else
	h2 = colorbar('horizontal');
end


hold on; coast('k'); dj_mpbdr(x,y); hold off

clabel(C,h1)



return
