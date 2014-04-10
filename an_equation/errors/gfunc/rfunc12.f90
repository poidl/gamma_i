real*8 function gamma(s,t)

parameter(n=12)

implicit real*8(a-h,o-z)

common /ifd/ icalled_from_driver

real*8 a(n)

save a


if(icalled_from_driver.eq.0) then
  open(10,file='d:/neutrals/ns04/an_equation/e12.dat',status='old')
  read(10,*) a
  close(10)
  icalled_from_driver = 1
end if

														   
ss = s/40; tt = t/30;

gnum = a(1) + a(2)*ss + a(3)*tt + a(4)*ss*ss + a(5)*ss*tt + a(6)*tt*tt

gnum = gnum+a(7)*ss**3.d0+a(8)*ss**2.d0*tt+a(9)*ss*tt**2.d0+a(10)*tt**3.d0

gden = 1.d0 + a(11)*ss + a(12)*tt

gamma = gnum/gden


		
return

end
