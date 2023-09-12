subroutine genafield_pump(pumptmp)
use modmain
use modconstants
use modpulse
implicit none
integer it
real(8) gs,t,t1,t2,t3
real(8), intent(inout) :: pumptmp(exitmax,3)
do it=1,exitmax

  t=(it-1)*(dt/2.d0)
  t1=t+tshift

  t2=cos(w_pump*t1)
  t3=cos(0.5*pi*t1/hd_pump)

  if(circular_pump) then
    gs=t3**2
  else
    gs=t2*(t3**2)
  endif
  
  if ((abs(t1).ge.hd_pump)) gs=0.d0
  if (abs(gs).lt.1.d-20) gs=0.d0

  if(circular_pump) then
    pumptmp(it,:)=afield_intensity_cpump*(/cos(w_cpump*t1),sin(w_cpump*t1),0.d0/)*gs
  else
    pumptmp(it,:)=afield_vec_pump(:)*gs
  endif

enddo
end subroutine
