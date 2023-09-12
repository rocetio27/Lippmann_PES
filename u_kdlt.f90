function u_kdlt(k_vec,add_vec,r_vec)
use modmain
use modconstants
implicit none
real(8) k_vec(3), add_vec(3), r_vec(3)
real(8) d_vec(3)
complex(8) u_kdlt
integer n1,n2

u_kdlt=c0

do n1=-10,10
do n2=-10,10
  d_vec=n1*a_1+n2*a_2+d_0*(/0d0,0d0,1d0/)+add_vec+dlt_vec(r_vec)
  u_kdlt=u_kdlt-hpp(d_vec)*pw3(-k_vec,d_vec)
enddo
enddo
end function
