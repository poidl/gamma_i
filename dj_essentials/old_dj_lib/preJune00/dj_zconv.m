function zout = dj_zconv(z)

%%%    dj_zconv		converts (nz,nx,ny) array to (nz,ny,nx) array
%%%
%%%    Usage:        zout = dj_zconv(z)
%%%
%%%    Input:        z - 3D array (nz,nx,ny) to be converted
%%%
%%%    Output:       zout - (nz,ny,nx) array
%%%
%%%    Assumptions:  all arrays are 3D
%%%
%%%    Author:       David Jackett
%%%
%%%    Date:         24/11/98



if length(size(z)) ~= 3, error('array not 3D in dj_zconv'), end 

zout = permute(z,[1 3 2]);


return
