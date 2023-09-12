function  phir(r_vec,j)
use modmain
implicit none
real(8) r_vec(3)
real(8) x,y,z,r
real(8) phir
integer j
x=r_vec(1)
y=r_vec(2)
z=r_vec(3)
r=dsqrt(x**2+y**2+z**2)

phir=1d0/4d0/dsqrt(2d0*pi)*dsqrt(atnum(j))**5*dexp(-atnum(j)*r/2d0)*z
end function phir

