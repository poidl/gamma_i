function dj_pltz(longs,lats,z)

%%%    dj_pltmp          plot a map of an array
%%%
%%%    Usage:        dj_pltmp(longs,lats,z)
%%%
%%%    Input:        longs & lats   - obvious (nx and ny)
%%%                  z              - array (ny by nx) to plot
%%%                  ipar           - 0, no new figure
%%%                                   1, new figure
%%%
%%%    Output:       a plot
%%%
%%%    Author:       David Jackett
%%%
%%%    Date:         20/11/96
%%%


figure

fpcolor(longs,lats,z)

colormap('jet'); colorbar('horizontal')

%set(gca,'DataAspectRatio',[1,1,1])

area_ave = nanmean(nanmean(z));

title(['average value: ',num2str(area_ave,4)]);

figure(gcf)


return
