subroutine sigp_surfaces(s,t,p,longs,lats,pr,sig_0,nx,ny,nz,s_sfce,t_sfce,p_sfce)

implicit real*8(a-h,o-z)

real*8, dimension(ny,nx) ::  longs, lats

real*8, dimension(nz,ny,nx) :: s,t,p,s_sfce,t_sfce,p_sfce

real*8, dimension(:,:,:), allocatable :: sigp


save called
	
data called/0/



!    allocate array sizes
      
if(called.eq.0) then
  allocate(sigp(1:nz,1:ny,1:nx))
  called = 1
endif


!						initialize variables


s_sfce = -99.d0; t_sfce = -99.d0; p_sfce = -99.d0;

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
		

do j0 = 1,ny
do i0 = 1,nx

	if(s(1,j0,i0).ne.-99.d0) then

		do k0 = 1,nz
			if(s(k0,j0,i0).ne.-99.d0) n_cast = k0
		end do
		
		call sigp_surface(s(1,j0,i0),t(1,j0,i0),p(1,j0,i0),sigp(1,j0,i0),n_cast,pr,sig_0,s_sfce,t_sfce,p_sfce)

	end if

end do
end do

				

	
return
end
