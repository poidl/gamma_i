	real function gamma(s,t,lat)

	parameter(n=19)

	implicit real*8(a-h,o-z)

	common /ifd/ icalled_from_driver

	real*8 a(n),s,t,lat

	save a


	if(icalled_from_driver.eq.0) then
		open(10,file='d:/neutrals/ns04/an_equation/e19.dat',status='old')
		read(10,*) a
		close(10)
		icalled_from_driver = 1
	end if

	s = s/40; t = t/30; lat = lat/90

	gamma = a(1)+a(2)*s+a(3)*t+a(4)*lat +                           &
			a(5)*s*s+a(6)*s*t+a(7)*s*lat +                          &
			a(8)*t*t+a(9)*t*lat+a(10)*lat*lat +                     &
			a(11)*s**3.d0+a(12)*s**2.d0*t+a(13)*s**2.d0*lat	+       &
			a(14)*t**3.d0+a(15)*t**2.d0*s+a(16)*t**2.d0*lat	+       &
 			a(17)*lat**3.d0+a(18)*lat**2.d0*s+a(19)*lat**2.d0*t	

	s = s*40; t = t*30; lat = lat*90

		
	return

	end
