function dj_pltmp(long,lat,z)

%    dj_pltmp      plot a map of an array
%
%    Usage:        dj_pltmp(longs,lats,z)
%
%    Input:        long 
%                  lats   
%                  z              - array (ny by nx) to plot
%
%    Output:       a plot
%
%    Author:       David Jackett
%
%    Date:         20/11/96          initial write
%                  14/10/99          fixed documentation
%


pcolor(long,lat,z)
shading interp
colormap('jet')
colorbar('horizontal')
hold on
coast('k')

% nx = length(longs); 
% ny = length(lats);
% dlo = 0.022*(longs(nx)-longs(1));
% dla = 0.034*(lats(ny)-lats(1));
% set(gca,'xlim',[longs(1)-dlo,longs(nx)+dlo]);
% set(gca,'ylim',[lats(1)-dla,lats(ny)+dla]);
% xx = [longs(1)-dlo,longs(nx)+dlo,longs(nx)+dlo,longs(1)-dlo,longs(1)-dlo];
% yy = [lats(1)-dla,lats(1)-dla,lats(ny)+dla,lats(ny)+dla,lats(1)-dla];
% plot(xx,yy,'k.')

hold off
axis([min(long(:)) max(long(:)) min(lat(:)) max(lat(:))])
figure(gcf)

return
