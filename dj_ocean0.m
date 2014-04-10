function ocean = dj_ocean0(long,lat)

%   dj_ocean0        ocean of single long/lat observation
%
%    Usage:     		ocean = 	dj_ocean0(long,lat)
%
%    Input:     		long  - 	longitude [0,360]
%					    lat   - 	latitude [-90,90] 
%
%    Output:    		ocean - 	0  land,
%                       			1-6  main oceans,
%                                   7    Arctic
%                                   8    Med
%
%    Author:    David Jackett
%
%    Date:      7/7/98
%

global po_long po_lat

ifig = 0;


%po_long = [100, 140, 240, 262, 268, 278.5, 280.6, 296 ,290 ,290 ,314 ,...
%      	  296, 290, 146, 146, 130.73, 126.5, 115.5, ...
%           106.75, 103, 101.7, 98.81, 100];
     
%po_lat = [20, 66, 66, 18, 16.7, 8, 9.4, -4, -50, -52, -61.5, -65.5, ...
%          -90, -90, -43, -12.38, -8.5, -8.5, -7, -4, ...
%          3.13, 10 , 20];

po_long = [100, 140, 240, 260, 272.59, 276.5, 278.65, 280.73, 295.217 ,290 , ...
      	  300, 294, 290, 146, 146, 133.9, 126.94, 123.62, 120.92, 117.42, 114.11, ...
           107.79, 102.57, 102.57, 98.79, 100];
     
po_lat = [20, 66, 66, 19.55, 13.97, 9.6, 8.1, 9.33, 0, -52, -64.5, -67.5, -90, -90, -41,...
          -12.48, -8.58, -8.39, -8.7, -8.82, -8.02, -7.04, -3.784 , 2.9, 10, 20];
    
    
na_long = [260, 272, 285, 310, 341.5, 353.25, 356, 360, 354.5, 360, ...
           360, 295.217, 280.73, 278.65, 276.5,  ...
		   272.59, 260, 260];
     
na_lat = [60, 72, 82, 81.75, 65, 62.1, 56.5, 44, 36, 20, ...
                                    0, 0, 9.33, 8.1, 9.6, 13.97, 19.55, 60];
    

if ifig==1
  plot(po_long,po_lat); hold on; coast('k')
  plot(na_long,na_lat,'g');
  dj_mpbdr([0,360],[-90,90]);
  set(gca,'DataAspectRatio',[1,1,1]); figure(gcf)
  hold off,dj_pause(0)
end


i_pacific = inpolygon(long,lat,po_long,po_lat);


% pacific ocean

if i_pacific == 1

  if lat<=0
	ocean = 2;
  else
	ocean = 1;
  end

% indian ocean

elseif 20<=long && long<=150 && -90<=lat && lat<=30
  
  if lat<=0
    ocean = 4;
  else
	ocean = 3;
  end
    
% atlantic ocean
     
elseif (0<=long && long<=20 && -90<=lat && lat<=90) || ...
       (20<=long && long<=40 && 28<=lat && lat<=44) || ...
       (260<=long && long<=360 && -90<=lat && lat<=90)
    
	i_natlantic = inpolygon(long,lat,na_long,na_lat);

	if lat<=0
	    ocean = 6;
	elseif (i_natlantic==1) || (0<=long && long<=15 && 0<lat && lat<=10)
	    ocean = 5;
    else
	    ocean = 7;
	end

% arctic ocean

elseif 0<=long && long<=360 && 64<lat && lat<=90
   
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


% red sea

if 31.25<=long && long<=43.25 && 13<=lat && lat<=30.1
   
	ocean = 9;
      
% persian gulf
     
elseif 40<=long && long<=56 && 22<=lat && lat<=32
   
	ocean = 9;

% baltic sea

elseif 12<=long && long<=33.5 && 50<=lat && lat<=60
	ocean = 9;
elseif 16<=long && long<=33.5 && 60<=lat && lat<=66
	ocean = 9;

% black/caspian seas

elseif 40<=long && long<=112 && 30<=lat && lat<=60
	ocean = 10;
elseif 27.5<=long && long<=40 && 40<=lat && lat<=50
	ocean = 10;

% african lakes
elseif 29<=long && long<=36 && -16<=lat && lat<=4
	ocean = 10;

% hudson bay
elseif 266<=long && long<=284 && 40<=lat && lat<=70
	ocean = 5;

% mediterranean sea
elseif 0<=long && long<=40 && 28<=lat && lat<=47
	ocean = 8;
elseif 355<=long && long<=360 && 34<=lat && lat<=42
	ocean = 8;  
end


return
