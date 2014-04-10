function dj_2Dbathy(x,z,data)

%%%    dj_2Dbath		2D bathymetric plot of section
%%%
%%%    Usage:	dj_bathy(x,z,data)
%%%
%%%    Description:	plot bathymetric profile of a section of data
%%%
%%%    Input:        x - horizontal coordinate of section
%%%			         z - depth of section
%%%
%%%	 Output:       plot of bathymetry
%%%
%%%    Author:       David Jackett
%%%
%%%    Date:         23/5/96
%%%


[dval,inds] = pm_lgood(data);

nx = length(x); xlims = [x(1),x(nx)];
nz = length(z); zlims = [z(1),z(nz)];

x = [x(1), x, x(nx)];

z = z([nz, inds, nz]);

%set(gca,'Color',[.5,0,1]);

set(gca,'XLim',xlims); set(gca,'YLim',zlims)

hf = fill(x,z,[1,0,1]);

set(hf,'EdgeColor',[0,1,1])

z'

xl = [xlims(2),xlims(2)]
zl = [zlims(2),z(length(z)-1)]

plot(xl,zl,'c')



return
