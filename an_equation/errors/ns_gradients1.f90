subroutine ns_gradients(s,t,p,g,ocean,longs,lats,opts,nx,ny,nz,nsx,nsy)

implicit real*8(a-h,o-z)
 
parameter(nlev=50)


real*8 longs(ny,nx),lats(ny,nx),p(nz,ny,nx),ocean(ny,nx)
real*8 s(nz,ny,nx),t(nz,ny,nx),g(nz,ny,nx)

real*8 s_east(nlev), t_east(nlev), p_east(nlev), g_east(nlev)
real*8 s_west(nlev), t_west(nlev), p_west(nlev), g_west(nlev)
real*8 s_north(nlev),t_north(nlev),p_north(nlev),g_north(nlev)
real*8 s_south(nlev),t_south(nlev),p_south(nlev),g_south(nlev)

real*8 pns_e(nlev),pns_w(nlev),pns_n(nlev),pns_s(nlev)

real*8 nsx(nz,ny,nx),nsy(nz,ny,nx)

real*8 opts(2)

data pi/3.14159265358979d0/, plevel/200.d0/



nsx = -99.d0; nsy = -99.d0; ng = 1

open(11,file='ns_ok.dat',status='unknown')
write(11,*) opts
close(11)


do j0 = 2,ny-1

do i0 = 2,nx-1

  denx = 111.2d3*cos(pi*lats(j0,i0)/180)*(longs(j0,i0+1)-longs(j0,i0-1))
  deny = 111.2d3*(lats(j0+1,i0)-lats(j0-1,i0))

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

 !open(21,file='junk.dat',status='unknown',access='append')
 !write(21,*) i0,j0,k0,1
 !do kk = 1,n_south
 !  write(21,*) p_south(kk),g_south(kk)
 !end do
 !write(21,*) g(k0,j0,i0),ng
 !close(21)

        call neutral_surfaces(s_south,t_south,p_south,g_south,n_south, &	
   			g(k0,j0,i0),ng,sns,tns,pns_s(k0),dsns,dtns,dpns)
        if(dpns.gt.0.0) pns_s(k0) = -99.d0

 !open(21,file='junk.dat',status='unknown',access='append')
 !write(21,*) i0,j0,k0,2
 !write(21,*) k0,pns_s(k0)
 !close(21)
  
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

 !open(21,file='junk.dat',status='unknown',access='append')
 !write(21,*) i0,j0,' done 3'
 !close(21)
 
! print *, i0,j0,' done 3'
! pause
 
  end if

 
end do
end do

	
return
end
