	subroutine sctp_interp(s,t,p,n,s0,t0,p0)
ccc
ccc
ccc
ccc	DESCRIPTION:	Interpolate salinity and in situ temperature
ccc					on a cast by linearly interpolating salinity
ccc					and conservative temperature
ccc
ccc	PRECISION:		Real*8
ccc
ccc	INPUT:			s(n)	array of cast salinities
ccc					t(n)	array of cast in situ temperatures
ccc					p(n)	array of cast pressures
ccc					n		length of cast
ccc					p0		pressure for which salinity and
ccc							in situ temperature are required
ccc
ccc	OUTPUT:			s0		interpolated value of salinity
ccc					t0		interpolated value of situ temperature
ccc
ccc	UNITS:			salinities		psu (IPSS-78)
ccc					temperatures	degrees C (ITS-90)
ccc					pressures		db
ccc
ccc
ccc	AUTHOR:			David Jackett
ccc
ccc	CREATED:		March 2007
ccc
ccc	REVISIONS:		
ccc
ccc
ccc
	implicit real*8(a-h,o-z)
	
	dimension s(n),t(n),p(n)

	external indx



	call indx(p,n,p0,k)

	r = (p0-p(k))/(p(k+1)-p(k))

	s0 = s(k) + r*(s(k+1)-s(k))

	ctk = ct_from_t(s(k),t(k),p(k))

	ct0 = ctk + r*(ct_from_t(s(k+1),t(k+1),p(k+1))-ctk)	

	t0 = t_from_ct(s0,ct0,p0)



	return
	end
