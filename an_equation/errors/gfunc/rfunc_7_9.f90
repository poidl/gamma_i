real*8 function gamma(s,t)

parameter(n=16)

implicit real*8(a-h,o-z)

common /ifd/ icalled_from_driver

real*8 a(n)

save a


if(icalled_from_driver.eq.0) then
  open(10,file='d:/neutrals/ness8/an_equation/rfunc_7_9.dat',status='old')
  read(10,*) a
  close(10)
  icalled_from_driver = 1
end if

														   
!ss = s/40; tt = t/30;

ss = s; tt = t;

gnum = a(1) + tt*(a(2)+tt*(a(3)+a(4)*tt))+ss*(a(5)+a(6)*tt+a(7)*ss)

gden = 1.d0 + tt*(a(8)+tt*(a(9)+tt*(a(10)+tt*a(11))))

gden = gden + ss*(a(12)+tt*(a(13)+tt*tt*a(14)) + sqrt(ss)*(a(15)+tt*tt*a(16)))

gamma = gnum/gden

		
return

end
