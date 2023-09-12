function afield_pump(t)
use modmain
use modconstants
implicit none
real(8) gs,t,t1,t2,t3
real(8) afield_pump(3)
  t1=t+tshift

  t2=cos(w_pump*t1)
  t3=cos(0.5*pi*t1/hd_pump)
  t2=cos(w_pump*t1)
  gs=t2*(t3**2)
  
  if ((abs(t1).ge.hd_pump)) gs=0.d0
  if (abs(gs).lt.1.d-20) gs=0.d0
  afield_pump(:)=afield_vec_pump(:)*gs
  
end function
