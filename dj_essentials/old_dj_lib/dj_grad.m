function [fz,fx,fy] = dj_grad(f,z,x,y)

%%%    dj_grad		three dimensional gradient (uneven grid spacing)
%%%
%%%    Usage:		[fz,fx,fy] = dj_grad(f,z,x,y)
%%%
%%%    Input:     f - three dimensional function
%%%			      z - vector of depths
%%%               longs - vector of longitudes
%%%               lats - vector of latitudes
%%%
%%%	 Output:    fz - z component of gradient
%%%               fx - x component of gradient
%%%               fy - y component of gradient
%%%
%%%    Author:    David Jackett
%%%
%%%    Date:      27/8/98
%%%

global dx

iwrite = 0;

factor = x(2)-x(1);


nx = length(x); ny = length(y); nz = length(z);




[fx,fz,fy] = gradient(f);




if iwrite == 1
  disp(' '); disp('  adjusting gradient ...'); dj_pause(1)
end


%     x-direction

dx0 = cos(pi*y/180); % size(dx0)

dx0 = dx0*ones(1,nx); % size(dx0)

dx = zeros(size(f));

for kk = 1:nz, dx(kk,:,:) = dx0'; end; % size(dx)

dx = factor*111.2e3*dx;


fx = fx./dx;


%     y-direction

dy = factor*111.2e3;

fy = fy/dy;


%     z-direction

R = 1; nxy = nx*ny; k = 1:nz;

yy = adh_diff(k,z,k,R)';

dzdk = yy(1:nz)*ones(1,nxy); dzdk = reshape(dzdk,nz,nx,ny);


fz = fz./dzdk;




if iwrite == 1, figure(1);
              subplot(2,1,1); plot(k,z); grid on; ylabel('z');
              subplot(2,1,2); plot(k,y); grid on; xlabel('k');
                                                  ylabel('dz/dk');

%       check the derivative

              figure(2);
              for x = z'

                if    x == z(1),  yint = z(1);
                else              nn = length(yint)+1;
                              yint = [yint;dj_quad(k(1:nn),y(1:nn))];    
                end;
              end;
              subplot(2,1,1); plot(z,yint); grid on;
              subplot(2,1,2); plot(z,yint-z); grid on;
end

