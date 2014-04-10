clear all, close all, clc, dj_tic

global resolution symbol

symbol = 'g'; resolution = 4;

load ../data/global_ocean

filename = ['data/intersections_east_', symbol, int2str(resolution), '.mat'];
    if exist(filename)~=0
	    cmd = ['load ',filename]; eval(cmd)
    else
	    ok = 'computing easterly intersections ...'
        [k_east,r_east] = intersections_east(s,eta,p,ocean,longs,lats);
    end
    k_east = change(k_east,'<=',0,1);
    k_east = change(k_east,'==',nan,1);
    dj_toc,    dj_pause(0)
    

filename = ['data/intersections_north_', symbol, int2str(resolution), '.mat'];
    if exist(filename)~=0
	    cmd = ['load ',filename]; eval(cmd)
    else
	    [k_north,r_north] = intersections_north(s,eta,p,ocean,longs,lats);
    end
    k_north = change(k_north,'<=',0,1);
    k_north = change(k_north,'==',nan,1);
    dj_pause(1)
    
    
filename = ['data/intersections_west_', symbol, int2str(resolution), '.mat'];
    if exist(filename)~=0
	    cmd = ['load ',filename]; eval(cmd)
    else
	    [k_west,r_west] = intersections_west(s,eta,p,ocean,longs,lats);
    end
    k_west = change(k_west,'<=',0,1);
    k_west = change(k_west,'==',nan,1);
    dj_pause(0)
    

filename = ['data/intersections_south_', symbol, int2str(resolution), '.mat'];
    if exist(filename)~=0
	    cmd = ['load ',filename]; eval(cmd)
    else
	    [k_south,r_south] = intersections_south(s,eta,p,ocean,longs,lats);
    end
    k_south = change(k_south,'<=',0,1);
    k_south = change(k_south,'==',nan,1);    
    

