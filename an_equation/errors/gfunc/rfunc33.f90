	real function gamma(s,t)

	parameter(n=33)

	implicit real*8(a-h,o-z)

	common /ifd/ icalled_from_driver

	real*8 a(n)

	save a


	if(icalled_from_driver.eq.0) then
		open(10,file='d:/neutrals/ns04/an_equation/e33.dat',status='old')
		read(10,*) a
		close(10)
		icalled_from_driver = 1
	end if

	s = s/40; t = t/30

	gamma = a(1)+a(2)*s+a(3)*t+   &
			a(4)*s*s+a(5)*s*t+a(6)*t*t+   &
			a(7)*s**3.d0+a(8)*(s**2.d0)*t+a(9)*s*(t**2.d0)+a(10)*t**3.d0 +    &
			a(11)*s**4.d0+a(12)*(s**3.d0)*t+a(13)*(s**2.d0)*(t**2.d0) +    &
			a(14)*s*(t**3.d0)+a(15)*t**4.d0 +    &
			a(16)*s**5.d0+a(17)*(s**4.d0)*t+a(18)*(s**3.d0)*(t**2.d0) +    &
			a(19)*(s**2.d0)*(t**3.d0)+a(20)*s*(t**4.d0)+a(21)*t**5.d0 +    &
			a(22)*s**6.d0+a(23)*(s**5.d0)*t+a(24)*(s**4.d0)*(t**2.d0) +    &
			a(25)*(s**3.d0)*(t**3.d0)+a(26)*(s**2.d0)*(t**4.d0) +    &
			a(27)*s*(t**5.d0)+a(28)*t**6.d0
			
	gamma = gamma/(1+a(29)*s+a(30)*t+a(31)*s*s+a(32)*s*t+a(33)*t*t) 

	s = s*40; t = t*30

		
	return

	end
