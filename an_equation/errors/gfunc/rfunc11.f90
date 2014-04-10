real*8 function gamma(s,t)

parameter(n=11)

implicit real*8(a-h,o-z)

common /ifd/ icalled_from_driver

real*8 a(n)

save a


if(icalled_from_driver.eq.0) then
  open(10,file='d:/neutrals/ns04/an_equation/e11.dat',status='old')
  read(10,*) a
  close(10)
  icalled_from_driver = 1
end if

														   
gnum = a(1)+a(2)*t+a(3)*t*t+a(4)*s+a(5)*s*t+a(6)*s**2.d0*t**1.d0

gden = 1+a(7)*t+a(8)*t*t+a(9)*s+a(10)*s*t+a(11)*s*sqrt(s)

gamma = gnum/gden
		
return

end
