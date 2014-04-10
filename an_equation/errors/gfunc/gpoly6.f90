	real function gamma(s,t)

	parameter(n=6)

	implicit real*8(a-h,o-z)

	common /ifd/ icalled_from_driver

	real*8 a(n),s,t

	save a


	if(icalled_from_driver.eq.0) then
		open(10,file='d:/neutrals/ns04/an_equation/e6.dat',status='old')
		read(10,*) a
		close(10)
		icalled_from_driver = 1
	end if

	ss = s/40; tt = t/30;

	gamma = a(1) + a(2)*ss + a(3)*tt + a(4)*ss*ss + a(5)*ss*tt + a(6)*tt*tt


		
	return

	end
