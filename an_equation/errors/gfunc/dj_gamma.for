	real function dj_gamma(s,th)

	real p15(15),p5(5),A(15)

	character*80 filename

	data pr0/0.0/, dj_gamma_called/0/

	save dj_gamma_called,p15,p5,p15_glo,p15_ghi,p5_hi


c					read data

	if(dj_gamma_called.eq.0) then
		filename = './gpoly/gamma_polynomials.dat'
		open(10,file=filename,status='old')
		read(10,*) (p15(k),k=1,15),(p5(k),k=1,5),p15_glo,p15_ghi,p5_hi
		close(10)
		dj_gamma_called = 1
!		print *, 'read dj_gamma data'
	end if


c	print *, p15
c	print *, p5
c	print *, p15_glo, p15_ghi, p5_hi


c					the 15-term polynomial


	A(1) = 1

	A(2) = th
	A(3) = A(2)*th
	A(4) = A(3)*th
	A(5) = A(4)*th
	A(6) = A(5)*th

	A(7) = s
	A(8) = A(7)*th
	A(9) = A(8)*th
	A(10) = A(9)*th
	A(11) = A(10)*th

	A(12) = s**1.5
	A(13) = A(12)*th
	A(14) = A(13)*th

	A(15) = s*s


	dj_gamma = DOT_PRODUCT(A,p15)

c	print *, 'g1',dj_gamma
	
      
	return
	
      gamma_15 = dj_gamma


c					+ the 5 term polynomial					

	if(p15_glo.lt.gamma_15.and.gamma_15.le.p15_ghi) then
		xx = gamma_15-p15_glo 
		yy = p5(1)
		do kk = 2,5
			yy = p5(kk)+xx*yy
		end do
		dj_gamma = dj_gamma+yy
	
c					+ the constant term			

	elseif(gamma_15.gt.p15_ghi)	then

		dj_gamma = dj_gamma+p5_hi

	end if

c	print *, 'g2',dj_gamma


c					and blend into sigma0

	band_lo = 23.0; band_hi = 24.0

	sdummy = svan(s,th,pr0,sig0)

c	print *, 'g3',dj_gamma,sig0

	if(dj_gamma.le.band_lo)	then
		dj_gamma = sig0
	elseif(band_lo.lt.dj_gamma.and.dj_gamma.le.band_hi)	then
		r = (dj_gamma-band_lo)/(band_hi-band_lo)
		dj_gamma = (1-r) * sig0 + r * dj_gamma
	end if


c	print *, 'g4',dj_gamma,sig0


	return

	end
