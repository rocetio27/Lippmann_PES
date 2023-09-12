subroutine intgrl_crtsn(p_vec,j,pes_rintegral_vec)
use modmain
use modconstants
implicit none
real(8) :: p_vec(3)
integer :: j,nx,ny,nz,ix,iy,iz
real(8) :: dx,dy,dz,lx,ly,lz,rx,ry,rz,params(6)
complex(8) :: pes_rintegral_vec(3)
complex(8) :: temp1(3), temp2(3), temp3(3), temp4(3)
nx=15
ny=15
nz=15
lx=-5
ly=-5
lz=-5
rx=-lx
ry=-ly
rz=-lz
dx=(rx-lx)/nx
dy=(ry-ly)/ny
dz=(rz-lz)/nz

params=(/lx,ly,lz,dx,dy,dz/)

temp1=c0
temp2=c0
temp3=c0
temp4=c0
pes_rintegral_vec=c0
do iz=0,nz-1
  do iy=0,ny-1
    do ix=0,nx-1
      temp1=temp1+(pes_rintegrand_vec(iz,iy,ix,params,p_vec,j)&
                  +pes_rintegrand_vec(iz,iy,ix+1,params,p_vec,j))*dx/2d0 
      
      temp2=temp2+(pes_rintegrand_vec(iz,iy+1,ix,params,p_vec,j)&
                  +pes_rintegrand_vec(iz,iy+1,ix+1,params,p_vec,j))*dx/2d0
    enddo
    temp3=temp3+(temp1+temp2)*dy/2d0
    temp1=c0
    temp2=c0   
  
    do ix=0,nx-1
      temp1=temp1+(pes_rintegrand_vec(iz+1,iy,ix,params,p_vec,j)&
                  +pes_rintegrand_vec(iz+1,iy,ix+1,params,p_vec,j))*dx/2d0 
      
      temp2=temp2+(pes_rintegrand_vec(iz+1,iy+1,ix,params,p_vec,j)&
                  +pes_rintegrand_vec(iz+1,iy+1,ix+1,params,p_vec,j))*dx/2d0
    enddo
    temp4=temp4+(temp1+temp2)*dy/2d0
    temp1=c0
    temp2=c0   
  enddo
  pes_rintegral_vec=pes_rintegral_vec+(temp3+temp4)*dz/2d0
  temp3=c0
  temp4=c0
enddo
contains
include 'pes_rintegrand_vec.f90'
include 'fnlstt.f90'
include 'phir.f90'
include 'pw.f90'
include 'pw3.f90'
include 'cdot.f90'
endsubroutine
