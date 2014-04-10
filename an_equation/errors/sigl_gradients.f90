subroutine sigl_gradients(s,t,p,ocean,longs,lats,nx,ny,nz,siglx,sigly)

implicit real*8(a-h,o-z)

parameter(nlev=50)


real*8 longs(ny,nx),lats(ny,nx),p(nz,ny,nx),ocean(ny,nx)
real*8 s(nz,ny,nx),t(nz,ny,nx)

real*8 s_east(nlev),t_east(nlev), p_east(nlev), g_east(nlev)
real*8 s_west(nlev),t_west(nlev), p_west(nlev), g_west(nlev)
real*8 s_north(nlev),t_north(nlev),p_north(nlev),g_north(nlev)
real*8 s_south(nlev),t_south(nlev),p_south(nlev),g_south(nlev)

real*8 sns_e(nlev),sns_w(nlev),sns_n(nlev),sns_s(nlev)
real*8 tns_e(nlev),tns_w(nlev),tns_n(nlev),tns_s(nlev)
real*8 pns_e(nlev),pns_w(nlev),pns_n(nlev),pns_s(nlev)

real*8 siglx(nz,ny,nx),sigly(nz,ny,nx)

data pi/3.14159265358979d0/, plevel/200.d0/, write_gradient_data/0/



open(11,file='sigl_ok.dat',status='unknown')
write(11,*) 'ok ok'
close(11)

print *, 'here sigl fortran ok'


if(write_gradient_data.eq.1) then
  open(21,file='gradient_data.dat',status='unknown')
end if

siglx = -99.0d0; sigly = -99.0d0


do j0 = 2,ny-1

do i0 = 2,nx-1

  denx = 111.2d3*cos(pi*lats(j0,i0)/180)*(longs(j0,i0+1)-longs(j0,i0-1))
  deny = 111.2d3*(lats(j0+1,i0)-lats(j0-1,i0))


  if(ocean(j0,i0).ge.1.d0.and.ocean(j0,i0).le.7.d0.and.      &
       s(1,j0,i0).ge.0.d0.and.s(1,j0,i0+1).ge.0.d0.and.      &	
 	     s(1,j0+1,i0).ge.0.d0.and.s(1,j0,i0-1).ge.0.d0.and.  &	
   	       s(1,j0-1,i0).ge.0.d0) then

    n_east = 0
    do k1 = 1,nz
      if(s(k1,j0,i0+1).ge.0.d0) then
        n_east = n_east+1 
        s_east(n_east) = s(k1,j0,i0+1)
        t_east(n_east) = t(k1,j0,i0+1)
        p_east(n_east) = p(k1,j0,i0+1)
      end if
    end do

    n_west = 0
    do k1 = 1,nz
      if(s(k1,j0,i0-1).ge.0.d0) then
        n_west = n_west+1 
        s_west(n_west) = s(k1,j0,i0-1)
        t_west(n_west) = t(k1,j0,i0-1)
        p_west(n_west) = p(k1,j0,i0-1)
      end if
    end do

    n_north = 0
    do k1 = 1,nz
      if(s(k1,j0+1,i0).ge.0.d0) then
        n_north = n_north+1 
        s_north(n_north) = s(k1,j0+1,i0)
        t_north(n_north) = t(k1,j0+1,i0)
        p_north(n_north) = p(k1,j0+1,i0)
      end if
    end do

    n_south = 0
    do k1 = 1,nz
      if(s(k1,j0-1,i0).ge.0.d0) then
        n_south = n_south+1 
        s_south(n_south) = s(k1,j0-1,i0)
        t_south(n_south) = t(k1,j0-1,i0)
        p_south(n_south) = p(k1,j0-1,i0)
      end if
    end do


!		now the work

    sns_e = -99.d0; tns_e = -99.d0; pns_e = -99.d0
    sns_n = -99.d0; tns_n = -99.d0; pns_n = -99.d0
    sns_w = -99.d0; tns_w = -99.d0; pns_w = -99.d0
    sns_s = -99.d0; tns_s = -99.d0; pns_s = -99.d0
    
    do k0 = 1,nz
      if(p(k0,j0,i0).ge.plevel.and.s(k0,j0,i0).ge.0.d0) then
        call depth_ns(s_east,t_east,p_east,n_east,s(k0,j0,i0),          &   
         	    t(k0,j0,i0),p(k0,j0,i0),sns_e(k0),tns_e(k0),pns_e(k0))
        call depth_ns(s_west,t_west,p_west,n_west,s(k0,j0,i0),          &   
         	    t(k0,j0,i0),p(k0,j0,i0),sns_w(k0),tns_w(k0),pns_w(k0))
        call depth_ns(s_north,t_north,p_north,n_north,s(k0,j0,i0),      &
             	t(k0,j0,i0),p(k0,j0,i0),sns_n(k0),tns_n(k0),pns_n(k0))
        call depth_ns(s_south,t_south,p_south,n_south,s(k0,j0,i0),      &
            	t(k0,j0,i0),p(k0,j0,i0),sns_s(k0),tns_s(k0),pns_s(k0))
      else
        pns_e(k0) = -99.d0
        pns_w(k0) = -99.d0
        pns_n(k0) = -99.d0
        pns_s(k0) = -99.d0
      end if
    end do

!          the x-gradient

    do k0 = 1,nz
      if(p(k0,j0,i0).ge.plevel.and.pns_e(k0).ge.0.d0.and.pns_w(k0).ge.0.d0) then
        siglx(k0,j0,i0) =  -(pns_e(k0)-pns_w(k0))/denx
        if(write_gradient_data.eq.1.and.k0.lt.nz.and.s(k0+1,j0,i0).ne.-99.d0) then
          write(21,'(5f15.1)') ocean(j0,i0),lats(j0,i0),p(k0-1,j0,i0),p(k0,j0,i0),p(k0+1,j0,i0)
          write(21,'(5e15.6)') sns_e(k0),sns_w(k0),s(k0-1,j0,i0),s(k0,j0,i0),s(k0+1,j0,i0)
          write(21,'(5e15.6)') tns_e(k0),tns_w(k0),t(k0-1,j0,i0),t(k0,j0,i0),t(k0+1,j0,i0)
        end if
      else
        siglx(k0,j0,i0) = -99.d0
      end if
    end do

!          the y-gradient

    do k0 = 1,nz
      if(p(k0,j0,i0).ge.plevel.and.pns_n(k0).ge.0.d0.and.pns_s(k0).ge.0.d0) then
        sigly(k0,j0,i0) =  -(pns_n(k0)-pns_s(k0))/deny
        if(write_gradient_data.eq.1.and.k0.lt.nz.and.s(k0+1,j0,i0).ne.-99.d0) then
          write(21,'(5f15.1)') ocean(j0,i0),lats(j0,i0),p(k0-1,j0,i0),p(k0,j0,i0),p(k0+1,j0,i0)
          write(21,'(5e15.6)') sns_n(k0),sns_s(k0),s(k0-1,j0,i0),s(k0,j0,i0),s(k0+1,j0,i0)
          write(21,'(5e15.6)') tns_n(k0),tns_s(k0),t(k0-1,j0,i0),t(k0,j0,i0),t(k0+1,j0,i0)
        end if
      else
        sigly(k0,j0,i0) = -99.d0
      end if
    end do

  end if

end do
end do

if(write_gradient_data.eq.1) then
  open(21,file='gradient_data.dat',status='unknown')
end if

	
return
end
