function ocean = dj_ocean0(long,lat)

%%%    dj_ocean0        ocean of single long/lat observation
%%%
%%%    Usage:     		ocean = 	dj_ocean0(long,lat)
%%%
%%%    Input:     		long  - 	longitude [0,360]
%%%							lat   - 	latitude [-90,90] 
%%%
%%%    Output:    		ocean - 	0  land,
%%%                       			1-6  main oceans,
%%%                       			>6   Arctic & marginal seas

%%%
%%%    Author:    David Jackett
%%%
%%%    Date:      7/7/98
%%%


ifig = 0;


po_long = [100, 140, 240, 262, 268, 278.5, 280.6, 296 ,290 ,290 ,314 ,...
      	  296, 290, 146, 146, 130.73, 126.5, 115.5, ...
           106.75, 103, 101.7, 98.81, 100];
     
po_lat = [20, 66, 66, 18, 16.7, 8, 9.4, -4, -50, -52, -61.5, -65.5, ...
          -90, -90, -43, -12.38, -8.5, -8.5, -7, -4, ...
          3.13, 10 , 20];
    
    
if ifig==1
  plot(po_long,po_lat); hold on; coast('k'); hold off
  dj_mpbdr([0,360],[-90,90]);
  set(gca,'DataAspectRatio',[1,1,1]); figure(gcf)
end


i_pacific = inpolygon(long,lat,po_long,po_lat);


%% pacific ocean

if i_pacific == 1
	  if lat<=0
	    ocean = 1;
	  else
	    ocean = 2;
	  end

%% indian ocean

elseif 20<=long & long<=150 & -90<=lat & lat<=30
   
	  if lat<=0
	    ocean = 3;
	  else
	    ocean = 4;
     end
    
%% atlantic ocean
     
elseif (0<=long & long<=20 & -90<=lat & lat<=64) | ...
       (20<=long & long<=40 & 28<=lat & lat<=44) | ...
       (260<=long & long<=360 & -90<=lat & lat<=64)
    
	if lat<=0
	    ocean = 5;
	else
	    ocean = 6;
	end

%% arctic ocean

elseif 0<=long & long<=360 & 64<lat & lat<=90
   
   ocean = 7;

   
else
   
   ocean = 8;
   
%   plot(po_long,po_lat); hold on; coast('k'); 
%   plot(long,lat,'rx'); hold off
%   dj_mpbdr([0,360],[-90,90]);set(gca,'DataAspectRatio',[1,1,1]);
%   grid on; zoom on; figure(gcf)
%   [long,lat]
%   error('impossible alternative')
   
end


%% red sea

if 31.25<=long & long<=43.25 & 13<=lat & lat<=30.1
   
	ocean = 10;
      
%% persian gulf
     
elseif 40<=long & long<=56 & 22<=lat & lat<=32
   
	ocean = 10;

%% baltic sea

elseif 10<=long & long<=60 & 50<=lat & lat<=60
	ocean = 9;
elseif 16<=long & long<=50 & 60<=lat & lat<=66
	ocean = 9;

%% black/caspian seas

elseif 40<=long & long<=112 & 30<=lat & lat<=60
   ocean = 8;
elseif 26<=long & long<=40 & 40<=lat & lat<=50
	ocean = 8;

% african lakes
elseif 20<=long & long<=38 & -16<=lat & lat<=20
	ocean = 10;

% hudson bay
elseif 250<=long & long<=284 & 40<=lat & lat<=70
	ocean = 7;

% mediterranean sea
elseif 0<=long & long<=40 & 28<=lat & lat<=47
	ocean = 8;
elseif 355<=long & long<=360 & 34<=lat & lat<=42
   ocean = 8;
   
end


return
