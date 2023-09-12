function pes_rintegrand_vec(iz,iy,ix,params,p_vec,j)
use modmain
implicit none
integer :: iz,iy,ix
real(8) :: params(6)
real(8) :: p_vec(3)
complex(8) :: pes_rintegrand_vec(3)
real(8) :: phase
integer :: j
real(8) :: r_vec(3)
real(8) :: lx, ly, lz, dx, dy, dz
lx=params(1)
ly=params(2)
lz=params(3)
dx=params(4)
dy=params(5)
dz=params(6)
r_vec=(/lx+dx*ix,ly+dy*iy,lz+dz*iz/)
pes_rintegrand_vec=phir(r_vec,j)*fnlstt(r_vec,p_vec,j)
end function pes_rintegrand_vec
