function afield_probe(t)
use modmain
use modconstants
implicit none
real(8) gs,t,t1,t2,t3
real(8) afield_probe(3)
  t1=t-delta_t+tshift

  t2=cos(w_probe*t1)
  t3=cos(0.5*pi*t1/hd_probe)
  t2=cos(w_probe*t1)
  gs=t2*(t3**2)
  
  if ((abs(t1).ge.hd_probe)) gs=0.d0
  if (abs(gs).lt.1.d-20) gs=0.d0
  afield_probe(:)=afield_vec_probe(:)*gs
end function
