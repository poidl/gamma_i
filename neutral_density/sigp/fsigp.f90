	real*8 function fsigp(p)

	implicit real*8(a-h,o-z)

 	common /fsig/  ss(2),tt(2),pp(2),sigg,prr

	data n2/2/


	if(p.eq.pp(1)) then
		s = ss(1); t = tt(1); ct = ct_from_t(s,t,p)
	elseif(p.eq.pp(2)) then
		s = ss(2); t = tt(2); ct = ct_from_t(s,t,p)
	else
		call sctp_interp(ss,tt,pp,n2,s,t,p)
	end if

 	sigp = rho_from_ct(s,ct,prr)-1000.d0


	fsigp = sigp-sigg


	
	return
	end