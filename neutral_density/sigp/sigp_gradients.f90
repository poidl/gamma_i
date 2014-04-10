subroutine sigp_gradients(s,t,p,longs,lats,pr,nx,ny,nz,sigpx,sigpy)

implicit real*8(a-h,o-z)

real*8, dimension(ny,nx) ::  longs, lats

real*8, dimension(nz,ny,nx) :: s,t,p,sigpx,sigpy

real*8, dimension(:,:,:), allocatable :: sigp

real*8, dimension(:), allocatable :: s_east, t_east, p_east, sigp_east
real*8, dimension(:), allocatable :: s_west, t_west, p_west, sigp_west
real*8, dimension(:), allocatable :: s_north, t_north, p_north, sigp_north
real*8, dimension(:), allocatable :: s_south, t_south, p_south, sigp_south

real*8, dimension(:), allocatable :: pns_e,pns_w,pns_n,pns_s

save called
	
data called/0/

data pi/3.14159265358979d0/, plevel/200.d0/




!    allocate array sizes
      
if(called.eq.0) then
  allocate(sigp(1:nz,1:ny,1:nx))
  allocate(s_east(1:nz));allocate(t_east(1:nz));allocate(p_east(1:nz));allocate(sigp_east(1:nz))
  allocate(s_west(1:nz));allocate(t_west(1:nz));allocate(p_west(1:nz));allocate(sigp_west(1:nz))
  allocate(s_north(1:nz));allocate(t_north(1:nz));allocate(p_north(1:nz));allocate(sigp_north(1:nz))
  allocate(s_south(1:nz));allocate(t_south(1:nz));allocate(p_south(1:nz));allocate(sigp_south(1:nz))
  allocate(pns_e(1:nz));allocate(pns_w(1:nz));allocate(pns_n(1:nz));allocate(pns_s(1:nz))
  called = 1
endif


!						initialize variables


sigpx = -99.d0; sigpy = -99.d0; nsigpy = 0; ncasts = 0


do j0 = 1,ny
do i0 = 1,nx
do k0 = 1,nz
	if(s(k0,j0,i0).gt.0.d0) then
		ct = ct_from_t(s(k0,j0,i0),t(k0,j0,i0),p(k0,j0,i0))
		sigp(k0,j0,i0) = rho_from_ct(s(k0,j0,i0),ct,pr)-1000
	end if
end do
end do
end do




!						main loop
		

do j0 = 2,ny-1

do i0 = 2,nx-1

	denx = 111.2d3*cos(pi*lats(j0,i0)/180)*(longs(j0,i0+1)-longs(j0,i0-1))

	deny = 111.2d3*(lats(j0+1,i0)-lats(j0-1,i0))


	if(s(1,j0,i0).ne.-99.d0.and.s(1,j0,i0+1).ne.-99.d0.and.			    &
   		s(1,j0+1,i0).ne.-99.d0.and.s(1,j0,i0-1).ne.-99.d0.and.			&
								    s(1,j0-1,i0).ne.-99.d0) then
		ncasts = ncasts+1
		n_east = 0
		do k1 = 1,nz
			if(s(k1,j0,i0+1).gt.0.d0) then
				n_east = n_east+1 
				s_east(n_east) = s(k1,j0,i0+1);
				t_east(n_east) = t(k1,j0,i0+1);
				p_east(n_east) = p(k1,j0,i0+1);
				sigp_east(n_east) = sigp(k1,j0,i0+1);
			end if
		end do

		n_west = 0
		do k1 = 1,nz
			if(s(k1,j0,i0-1).gt.0.d0) then
				n_west = n_west+1 
				s_west(n_west) = s(k1,j0,i0-1);
				t_west(n_west) = t(k1,j0,i0-1);
				p_west(n_west) = p(k1,j0,i0-1);
				sigp_west(n_west) = sigp(k1,j0,i0-1);
			end if
		end do

		n_north = 0
		do k1 = 1,nz
			if(s(k1,j0+1,i0).gt.0.d0) then
				n_north = n_north+1 
				s_north(n_north) = s(k1,j0+1,i0);
				t_north(n_north) = t(k1,j0+1,i0);
				p_north(n_north) = p(k1,j0+1,i0);
				sigp_north(n_north) = sigp(k1,j0+1,i0);
			end if
		end do

		n_south = 0
		do k1 = 1,nz
			if(s(k1,j0-1,i0).gt.0.d0) then
				n_south = n_south+1 
				s_south(n_south) = s(k1,j0-1,i0);
				t_south(n_south) = t(k1,j0-1,i0);
				p_south(n_south) = p(k1,j0-1,i0);
				sigp_south(n_south) = sigp(k1,j0-1,i0);
			end if
		end do

!		if(j0.eq.102.and.i0.eq.99) then 
!			open(21,file='test.dat',status='unknown')
!			do kk = 1,nz
!				write(21,*) sigp_east(kk),sigp(kk,j0,i0),sigp_east(kk)
!			end do
!			close(21)
!		end if


!										now the work

		do k0 = 1,nz


		  if(p(k0,j0,i0).ge.plevel.and.s(k0,j0,i0).gt.0.d0) then

		    sigp0 = sigp(k0,j0,i0)

			call sigp_surface(s_east,t_east,p_east,sigp_east,n_east,pr,sigp0,sns,tns,pns_e(k0))

			call sigp_surface(s_west,t_west,p_west,sigp_west,n_west,pr,sigp0,sns,tns,pns_w(k0))

			call sigp_surface(s_north,t_north,p_north,sigp_north,n_north,pr,sigp0,sns,tns,pns_n(k0))

			call sigp_surface(s_south,t_south,p_south,sigp_south,n_south,pr,sigp0,sns,tns,pns_s(k0))

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
				sigpx(k0,j0,i0) =  -(pns_e(k0)-pns_w(k0))/denx
			else
				sigpx(k0,j0,i0) = -99.d0
			end if
		end do

!						the y-gradient

		do k0 = 1,nz
			if(p(k0,j0,i0).ge.plevel.and.pns_n(k0).gt.0.d0.and.pns_s(k0).gt.0.d0) then
				sigpy(k0,j0,i0) =  -(pns_n(k0)-pns_s(k0))/deny
				nsigpy = nsigpy+1
			else
				sigpy(k0,j0,i0) = -99.d0
			end if
		end do

	end if

end do
end do


	
return
end
