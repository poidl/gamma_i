subroutine ns_gradients(s,t,p,g,ocean,longs,lats,nx,ny,nz,nsx,nsy)

implicit real*8(a-h,o-z)

real*8 longs(ny,nx),lats(ny,nx),ocean(ny,nx)
real*8 s(nz,ny,nx),t(nz,ny,nx),p(nz,ny,nx),g(nz,ny,nx)

real*8, dimension(:), allocatable :: s_east, t_east, p_east, g_east
real*8, dimension(:), allocatable :: s_west, t_west, p_west, g_west
real*8, dimension(:), allocatable :: s_north,t_north,p_north,g_north
real*8, dimension(:), allocatable :: s_south,t_south,p_south,g_south

real*8, dimension(:), allocatable :: pns_e,pns_w,pns_n,pns_s

real*8 nsx(nz,ny,nx),nsy(nz,ny,nx)


data pi/3.14159265358979d0/, plevel/200.d0/

      
allocate(s_east(1:nz)); allocate(t_east(1:nz))
allocate(p_east(1:nz)); allocate(g_east(1:nz))
allocate(s_west(1:nz)); allocate(t_west(1:nz))
allocate(p_west(1:nz)); allocate(g_west(1:nz))
allocate(s_north(1:nz)); allocate(t_north(1:nz))
allocate(p_north(1:nz)); allocate(g_north(1:nz))
allocate(s_south(1:nz)); allocate(t_south(1:nz))
allocate(p_south(1:nz)); allocate(g_south(1:nz))
allocate(pns_e(1:nz)); allocate(pns_w(1:nz))
allocate(pns_n(1:nz)); allocate(pns_s(1:nz))


nsx = -99.d0; nsy = -99.d0; ng = 1

do j0 = 2,ny-1
do i0 = 2,nx-1
  deny = 111.2d3*(lats(j0+1,i0)-lats(j0-1,i0))
  denx = 111.2d3*cos(pi*lats(j0,i0)/180)*(longs(j0,i0+1)-longs(j0,i0-1))

  if((ocean(j0,i0).ge.1.d0.and.ocean(j0,i0).le.6.d0).and.  &
       s(1,j0,i0).ge.0.d0.and.s(1,j0,i0+1).ge.0.d0.and. &
        s(1,j0+1,i0).ge.0.d0.and.s(1,j0,i0-1).ge.0.d0.and.     &
         s(1,j0-1,i0).ge.0.d0)    then

    n_east = 0
    do k1 = 1,nz
      if(g(k1,j0,i0+1).ge.0.d0) then
        n_east = n_east+1 
        s_east(n_east) = s(k1,j0,i0+1)
        t_east(n_east) = t(k1,j0,i0+1)
        p_east(n_east) = p(k1,j0,i0+1)
        g_east(n_east) = g(k1,j0,i0+1)
      end if
    end do

    n_west = 0
    do k1 = 1,nz
      if(g(k1,j0,i0-1).ge.0.d0) then
        n_west = n_west+1 
        s_west(n_west) = s(k1,j0,i0-1)
        t_west(n_west) = t(k1,j0,i0-1)
        p_west(n_west) = p(k1,j0,i0-1)
        g_west(n_west) = g(k1,j0,i0-1)
      end if
    end do

    n_north = 0
    do k1 = 1,nz
      if(g(k1,j0+1,i0).ge.0.d0) then
        n_north = n_north+1 
        s_north(n_north) = s(k1,j0+1,i0)
        t_north(n_north) = t(k1,j0+1,i0)
        p_north(n_north) = p(k1,j0+1,i0)
        g_north(n_north) = g(k1,j0+1,i0)
      end if
    end do

    n_south = 0
    do k1 = 1,nz
      if(g(k1,j0-1,i0).ge.0.d0) then
        n_south = n_south+1 
        s_south(n_south) = s(k1,j0-1,i0)
        t_south(n_south) = t(k1,j0-1,i0)
        p_south(n_south) = p(k1,j0-1,i0)
        g_south(n_south) = g(k1,j0-1,i0)
      end if
    end do


!        				now the work

    do k0 = 1,nz
 
      if(p(k0,j0,i0).ge.plevel.and.g(k0,j0,i0).ge.0.d0) then
  
        call neutral_surfaces(s_east,t_east,p_east,g_east,n_east, &
                        g(k0,j0,i0),ng,sns,tns,pns_e(k0),dsns,dtns,dpns)
	    if(dpns.gt.0.0) pns_e(k0) = -99.d0        


	    call neutral_surfaces(s_west,t_west,p_west,g_west,n_west, &	
  			            g(k0,j0,i0),ng,sns,tns,pns_w(k0),dsns,dtns,dpns)
        if(dpns.gt.0.0) pns_w(k0) = -99.d0

 
	    call neutral_surfaces(s_north,t_north,p_north,g_north,n_north, &	
  			            g(k0,j0,i0),ng,sns,tns,pns_n(k0),dsns,dtns,dpns)
	    if(dpns.gt.0.0) pns_n(k0) = -99.d0

        call neutral_surfaces(s_south,t_south,p_south,g_south,n_south, &	
   			            g(k0,j0,i0),ng,sns,tns,pns_s(k0),dsns,dtns,dpns)
        if(dpns.gt.0.0) pns_s(k0) = -99.d0
  
      else
        pns_e(k0) = -99.d0
        pns_w(k0) = -99.d0
        pns_n(k0) = -99.d0
        pns_s(k0) = -99.d0
      end if
   
    end do


!        	the x-gradient

    do k0 = 1,nz
      if(p(k0,j0,i0).ge.plevel.and.pns_e(k0).gt.0.and.pns_w(k0).gt.0.0) then
        nsx(k0,j0,i0) = -(pns_e(k0)-pns_w(k0))/denx
      else
        nsx(k0,j0,i0) = -99.d0
       end if
      end do


!        	the y-gradient

    do k0 = 1,nz
      if(p(k0,j0,i0).ge.plevel.and.pns_n(k0).gt.0.and.pns_s(k0).gt.0.0) then
        nsy(k0,j0,i0) =  -(pns_n(k0)-pns_s(k0))/deny
      else
        nsy(k0,j0,i0) = -99.d0
      end if
    end do
 
  end if

end do
end do


!		deallocate arrays

deallocate(s_east, t_east, p_east, g_east, stat=ierr)
deallocate(s_west, t_west, p_west, g_west, stat=ierr)
deallocate(s_north,t_north,p_north,g_north, stat=ierr)
deallocate(s_south,t_south,p_south,g_south, stat=ierr)

deallocate(pns_e,pns_w,pns_n,pns_s, stat=ierr)

	
return
end
