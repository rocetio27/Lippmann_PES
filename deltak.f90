function deltak(k_vec)
use modmain
use modconstants
implicit none
real(8) k_vec(3)
complex(8) deltak
real(8) phase1
integer i
deltak=0.d0

do i=1,3
  phase1=nnv(i,1)*k_vec(1)+nnv(i,2)*k_vec(2)
  deltak=deltak-(dcos(phase1)+ci*dsin(phase1))
enddo
end function
