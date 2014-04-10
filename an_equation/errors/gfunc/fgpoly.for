	real function  fgpoly(p)

 	common /fgp/  ss(2),tt(2),pp(2),gpp


	data n2/2/


	if(p.eq.pp(1)) then
		s = ss(1)
		t = tt(1)
	elseif(p.eq.pp(2)) then
		s = ss(2)
		t = tt(2)
	else
		call stp_interp(ss,tt,pp,n2,s,t,p)
	end if


!	th = theta(s,t,p,pr0)

	gpoly = gamma(s,t)


	fgpoly = gpoly-gpp


	
	return
	end