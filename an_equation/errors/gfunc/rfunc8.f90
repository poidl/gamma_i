real*8 function gamma(s,t)

parameter(n=8)

implicit real*8(a-h,o-z)

common /ifd/ icalled_from_driver

real*8 a(n)

save a


if(icalled_from_driver.eq.0) then
  open(10,file='d:/neutrals/ns06/fitting_gui/rfunc_th.dat',status='old')
  read(10,*) a
  close(10)
  icalled_from_driver = 1
end if

														   
ss = s/40; tt = t/30;

gnum = a(1) + a(2)*ss + a(3)*tt + a(4)*ss*ss + a(5)*ss*tt + a(6)*tt*tt

gden = 1.d0 + a(7)*ss + a(8)*tt

gamma = gnum/gden

!s = s*40; t = t*30
		
return

end
