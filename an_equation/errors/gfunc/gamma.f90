real function gamma(s,t,p)

parameter(n=15)

implicit real*8(a-h,o-z)


real*8 a(n)

data icalled/0/

save a, icalled



if(icalled.eq.0) then
	open(10,file='d:/neutrals/ns03/an_equation/e15.out',status='old')
	read(10,*) a
	close(10)
	icalled = 1
end if



s1on2 = s**0.5; s3on2 = s1on2*s; st = s*t;

t2 = t*t; t3 = t2*t; t4 = t*t3; t5 = t2*t3

st2 = s*t2; st3 = s*t3; st4 = s*t4; s2 = s*s

p2 = p*p; p3 = p2*p;


gamma = a(1)+a(2)*t+a(3)*t2+a(4)*t3+a(5)*t4+a(6)*t5+ &
			a(7)*s+a(8)*st+a(9)*st2+a(10)*st3+a(11)*st4+ &
				s3on2*(a(12)+a(13)*t+a(14)*t2)+a(15)*s2

		
		
return

end
