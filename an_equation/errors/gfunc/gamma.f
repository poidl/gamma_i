	real function gamma(s,t,p)

	parameter(n=16)

	implicit real*8(a-h,o-z)

	common /ifd/ icalled_from_driver

	real*8 a(n)

	save a



	if(icalled_from_driver.eq.0) then
		open(10,file='d:/neutrals/ns04/an_equation/gpoly.dat',
     &												status='old')
		read(10,*) a
		close(10)
	    print *, a
		icalled_from_driver = 1
	end if

	ct = ct_from_t(s,t,p)


	gnum = a(1)+ct*(a(2)+ct*(a(3)+a(4)*ct))+s*(a(5)+a(6)*ct+s*a(7))

	gden = 1+a(8)*ct+a(9)*ct**2.0+a(10)*ct**3.0+a(11)*ct**4.0

	gden = gden+a(12)*s+a(13)*s*ct+a(14)*s*ct**3.0
	
      gden = gden+a(15)*s**1.5+a(16)*s**1.5*ct**2.0

	gamma = gnum/gden
		
		
	  return

	  end
