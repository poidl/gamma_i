function dj_mpbdr(longs,lats)

%%%    dj_mpbdr          plot a border around a map 
%%%
%%%    Usage:    			 dj_mpbdr(longs,lats)
%%%
%%%    Input:            longs	- longitudes of data range
%%%                      lats	- latitudes of data range
%%%
%%%    Output:           a border
%%%
%%%    Author:           David Jackett
%%%
%%%    Date:             4/11/96
%%%



nx = length(longs); ny = length(lats);

dlo = 0.022*(longs(nx)-longs(1));
dla = 0.034*(lats(ny)-lats(1));
		
set(gca,'xlim',[longs(1)-dlo,longs(nx)+dlo]);
set(gca,'ylim',[lats(1)-dla,lats(ny)+dla]);

xx = [longs(1)-dlo,longs(nx)+dlo,longs(nx)+dlo,longs(1)-dlo,longs(1)-dlo];
yy = [lats(1)-dla,lats(1)-dla,lats(ny)+dla,lats(ny)+dla,lats(1)-dla];

plot(xx,yy,'k')


return
