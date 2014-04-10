function dj_pltmpv(longs,lats,z,ifig)

%%%    dj_pltmp          plot a map of an array
%%%
%%%    Usage:        dj_pltmp(longs,lats,z,ifig)
%%%
%%%    Input:        longs & lats   - obvious (nx and ny)
%%%                  z              - array (ny by nx) to plot
%%%                  ifig           - 0, no new figure
%%%                                   1, new figure
%%%
%%%    Output:       a plot
%%%
%%%    Author:       David Jackett
%%%
%%%    Date:         20/11/96          initial write
%%%                  14/10/99          fixed documentation
%%%


if nargin ~=4, ifig = 1; end

%if	ifig==2
%	zz = get(gcf,'CurrentAxes')
%	n = length(get(gcf,'CurrentAxes'));
%	if n>0, figure; end
%end

fpcolor(longs,lats,z)

colormap('jet');  colorbar('vertical')

hold on; 

axis([longs(1) longs(length(longs)) lats(1) lats(length(lats))])      

coast('k')

dj_mpbdr(longs,lats)

set(gca,'DataAspectRatio',[1,1,1])

%area_ave = dj_arave(longs,lats,z);

%title(['area average : ',num2str(area_ave,4)]);

hold off

figure(gcf)


return
