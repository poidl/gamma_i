	subroutine  gpoly_surfaces(s,t,p,gpoly,n,gpoly0,pns)

	common /fgp/  ss(2),tt(2),pp(2),gpp


	real s(n),t(n),p(n),gpoly(n)


	external zbren, fgpoly


	data errabs/0.00005/, errrel/0.00001/


c
c							find all crossings
c

	ncr = 0

	do k = 1,n-1
		gpoly_lo = min(gpoly(k),gpoly(k+1))
		gpoly_hi = max(gpoly(k),gpoly(k+1))
		if(gpoly_lo.le.gpoly0.and.gpoly0.lt.gpoly_hi) then
			ncr = ncr+1
			k0 = k
		end if
	end do

	if(gpoly0.eq.gpoly(n)) then
		ncr = ncr+1
		if(ncr.eq.1) then
			pns = p(n)
			return
		end if
	end if


c
c							find single crossing
c
	
	if(ncr.eq.1) then
	  if(gpoly0.eq.gpoly(k0)) then
	    pns = p(k0)
	  elseif(gpoly0.eq.gpoly(k0+1)) then
		pns = p(k0+1)
	  else
		ss(1) = s(k0); tt(1) = t(k0); pp(1) = p(k0)
		ss(2) = s(k0+1); tt(2) = t(k0+1); pp(2) = p(k0+1)
		gpp = gpoly0
		a = pp(1); b = pp(2); maxfn = 100
		fa = fgpoly(a); fb = fgpoly(b);
		if(fa*fb.le.0) then
			call dzbren(fgpoly,errabs,errrel,a,b,maxfn)
			pns = b
		else
			pns = -99.0
!			print *, 'oh no',fa,fb
		end if
	  end if
	else
		pns = -99.0
	end if



	return
	end
