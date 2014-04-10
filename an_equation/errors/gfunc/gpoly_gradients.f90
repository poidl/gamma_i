subroutine gpoly_gradients(s,t,p,longs,lats,iocean,nx,ny,nz,gpolyx,gpolyy)

implicit real*8(a-h,o-z)

integer, dimension(ny,nx) :: iocean

real*8, dimension(ny,nx) ::  longs
real*8, dimension(ny,nx) ::  lats
real*8, dimension(nz,ny,nx) ::  p

real*8, dimension(nz,ny,nx) :: s,t,gpolyx,gpolyy

real*8, dimension(:,:,:), allocatable :: gpoly

real*8, dimension(:), allocatable :: s_east, t_east, p_east, gpoly_east
real*8, dimension(:), allocatable :: s_west, t_west, p_west, gpoly_west
real*8, dimension(:), allocatable :: s_north, t_north, p_north, gpoly_north
real*8, dimension(:), allocatable :: s_south, t_south, p_south, gpoly_south

real*8, dimension(:), allocatable :: pns_e,pns_w,pns_n,pns_s

save called
	
data called/0/

data pi/3.14159265358979d0/, plevel/200/



!    allocate array sizes
      
if(called.eq.0) then
  allocate(gpoly(1:nz,1:ny,1:nx))
  allocate(s_east(1:nz));allocate(t_east(1:nz));allocate(p_east(1:nz));allocate(gpoly_east(1:nz))
  allocate(s_west(1:nz));allocate(t_west(1:nz));allocate(p_west(1:nz));allocate(gpoly_west(1:nz))
  allocate(s_north(1:nz));allocate(t_north(1:nz));allocate(p_north(1:nz));allocate(gpoly_north(1:nz))
  allocate(s_south(1:nz));allocate(t_south(1:nz));allocate(p_south(1:nz));allocate(gpoly_south(1:nz))
  allocate(pns_e(1:nz));allocate(pns_w(1:nz));allocate(pns_n(1:nz));allocate(pns_s(1:nz))
  called = 1
endif


!						initialize variables


gpolyx = -99.d0; gpolyy = -99.d0; ng = 1;

nsum = 0; gsum = 0.d0;

do j0 = 1,ny
do i0 = 1,nx
do k0 = 1,nz
    if(s(k0,j0,i0).gt.0.0) then
!		th0 = theta(s(k0,j0,i0),t(k0,j0,i0),p(k0,j0,i0),pr0)
		gpoly(k0,j0,i0) = gamma(s(k0,j0,i0),t(k0,j0,i0))
		nsum = nsum+1; gsum = gsum+gpoly(k0,j0,i0)
	end if
end do
end do
end do

gmean = gsum/nsum

open(11,file='gmean.dat',status='unknown')
write(11,*) gmean
close(11)


!						main loop
		

do j0 = 2,ny-1
!print *, j0
do i0 = 2,nx-1

	denx = 111.2d3*cos(pi*lats(j0,i0)/180)*(longs(j0,i0+1)-longs(j0,i0-1))

	deny = 111.2d3*(lats(j0+1,i0)-lats(j0-1,i0))

	if(iocean(j0,i0).ge.1.and.iocean(j0,i0).le.8.and.                   &
!	    77.5d0.le.longs(j0,i0).and.longs(j0,i0).lt.150.d0.and.				&
!	     -90.d0.le.lats(j0,i0).and.lats(j0,i0).le.-60.d0.and.					&
	       s(1,j0,i0).ne.-99.d0.and.s(1,j0,i0+1).ne.-99.d0.and.		    &
   		    s(1,j0+1,i0).ne.-99.d0.and.s(1,j0,i0-1).ne.-99.d0.and.	   	&
								    s(1,j0-1,i0).ne.-99.d0) then

		n_east = 0
		do k1 = 1,nz
			if(s(k1,j0,i0+1).gt.0.d0) then
				n_east = n_east+1 
				s_east(n_east) = s(k1,j0,i0+1);
				t_east(n_east) = t(k1,j0,i0+1);
				p_east(n_east) = p(k1,j0,i0+1);
				gpoly_east(n_east) = gpoly(k1,j0,i0+1);
			end if
		end do

		n_west = 0
		do k1 = 1,nz
			if(s(k1,j0,i0-1).gt.0.d0) then
				n_west = n_west+1 
				s_west(n_west) = s(k1,j0,i0-1);
				t_west(n_west) = t(k1,j0,i0-1);
				p_west(n_west) = p(k1,j0,i0-1);
				gpoly_west(n_west) = gpoly(k1,j0,i0-1);
			end if
		end do

		n_north = 0
		do k1 = 1,nz
		if(s(k1,j0+1,i0).gt.0.d0) then
				n_north = n_north+1 
				s_north(n_north) = s(k1,j0+1,i0);
				t_north(n_north) = t(k1,j0+1,i0);
				p_north(n_north) = p(k1,j0+1,i0);
				gpoly_north(n_north) = gpoly(k1,j0+1,i0);
			end if
		end do

		n_south = 0
		do k1 = 1,nz
			if(s(k1,j0-1,i0).gt.0.d0) then
				n_south = n_south+1 
				s_south(n_south) = s(k1,j0-1,i0);
				t_south(n_south) = t(k1,j0-1,i0);
				p_south(n_south) = p(k1,j0-1,i0);
				gpoly_south(n_south) = gpoly(k1,j0-1,i0);
			end if
		end do


!										now the work

		do k0 = 1,nz

		  if(p(k0,j0,i0).ge.plevel.and.s(k0,j0,i0).gt.0.d0) then

		    gpoly0 = gpoly(k0,j0,i0)

			call gpoly_surfaces(s_east,t_east,p_east,gpoly_east,n_east,gpoly0,pns_e(k0))

			call gpoly_surfaces(s_west,t_west,p_west,gpoly_west,n_west,gpoly0,pns_w(k0))

			call gpoly_surfaces(s_north,t_north,p_north,gpoly_north,n_north,gpoly0,pns_n(k0))

			call gpoly_surfaces(s_south,t_south,p_south,gpoly_south,n_south,gpoly0,pns_s(k0))
		  else
			pns_e(k0) = -99.d0
			pns_w(k0) = -99.d0
			pns_n(k0) = -99.d0
			pns_s(k0) = -99.d0
          end if

	  	end do

!						the x-gradient

		do k0 = 1,nz
			if(p(k0,j0,i0).ge.plevel.and.pns_e(k0).gt.0.d0.and.pns_w(k0).gt.0.d0) then
				gpolyx(k0,j0,i0) =  -(pns_e(k0)-pns_w(k0))/denx
			else
				gpolyx(k0,j0,i0) = -99.d0
			end if
		end do

!						the y-gradient

		do k0 = 1,nz
			if(p(k0,j0,i0).ge.plevel.and.pns_n(k0).gt.0.d0.and.pns_s(k0).gt.0.d0) then
				gpolyy(k0,j0,i0) =  -(pns_n(k0)-pns_s(k0))/deny
			else
				gpolyy(k0,j0,i0) = -99.d0
			end if
		end do

	end if

end do
end do

!gpolyx = gpoly

	
return
end
