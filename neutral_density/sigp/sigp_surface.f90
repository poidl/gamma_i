	subroutine  sigp_surface(s,t,p,sig,n,pr,sigp,sns,tns,pns)

	implicit real*8(a-h,o-z)

	common /fsig/  ss(2),tt(2),pp(2),sigg,prr

	real*8 s(n),t(n),p(n),sig(n),pr,sigp,sns,tns,pns

	external dzbren, fsigp


	data errabs/5.d-5/, errrel/1.d-5/, n2/2/
	

!
!							find all crossings
!


		
	ncr = 0

	do k = 1,n-1
		sig_lo = min(sig(k),sig(k+1)); sig_hi = max(sig(k),sig(k+1))
		if(sig_lo.le.sigp.and.sigp.lt.sig_hi) then
			ncr = ncr+1
			k0 = k
		end if
	end do

	if(sigp.eq.sig(n)) then
		ncr = ncr+1
		k0 = k
		if(ncr.eq.1) then
			pns = p(k0)
			return
		end if
	end if


!
!							find single crossing
!
	
	if(ncr.eq.1) then
	  if(sigp.eq.sig(k0)) then
	    sns = s(k0); tns = t(k0); pns = p(k0)
	  elseif(sigp.eq.sig(k0+1)) then
		sns = s(k0+1); tns = t(k0+1); pns = p(k0+1)
	  else
		ss(1) = s(k0); tt(1) = t(k0); pp(1) = p(k0)
		ss(2) = s(k0+1); tt(2) = t(k0+1); pp(2) = p(k0+1)
		sigg = sigp; prr = pr
		a = pp(1); b = pp(2); maxfn = 100
!	print *, sig(k0),sig(k0+1),sigp

		call dzbren(fsigp,errabs,errrel,a,b,maxfn)
		pns = b
		call sctp_interp(ss,tt,pp,n2,sns,tns,pns)
	  end if
	else
		sns = -99.d0; tns = -99.d0; pns = -99.d0
	end if



	return
	end
