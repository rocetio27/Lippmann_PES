subroutine genafield_probe(probetmp)
use modmain
use modconstants
use modpulse
implicit none
integer it
real(8) gs,t,t1,t2,t3
real(8), intent(inout) :: probetmp(exitmax,3)
do it=1,exitmax
  t=(it-1)*(dt/2.d0)
  t1=t-delta_t+tshift

  t2=cos(w_probe*t1)
  t3=cos(0.5*pi*t1/hd_probe)
  
  if(circular_probe) then
    gs=t3**2
  else
    gs=t2*(t3**2)
  endif
  
  if ((abs(t1).ge.hd_probe)) gs=0.d0
  if (abs(gs).lt.1.d-20) gs=0.d0
  
  if(circular_probe) then
    probetmp(it,:)=afield_intensity_cprobe*(/cos(w_cprobe*t1),sin(w_cprobe*t1),0.d0/)*gs
  else
    probetmp(it,:)=afield_vec_probe(:)*gs
  endif
enddo
end subroutine
