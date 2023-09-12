function phik(k_vec)
implicit none
real(8) :: k_vec(3)
real(8) :: pi
integer :: nr,ntheta,nphi,ir,itheta,iphi
real(8) :: dr,dtheta,dphi,rng_r,rng_theta,rng_phi
complex(8) :: phase_factor
complex(8) :: temp1,temp2,temp3,temp4,phik
pi=4.d0*ATAN(1.d0)
!!!!!!!!!! Quality 
rng_r=10.d0
rng_theta=pi
rng_phi=2*pi

nr=50
ntheta=50
nphi=50
!!!!!!!!!! Quality
dr=rng_r/nr
dtheta=rng_theta/ntheta
dphi=rng_phi/nphi

temp1=0.d0
temp2=0.d0
temp3=0.d0
temp4=0.d0

phik=0.d0
do ir=0,nr-1
  do itheta=0,ntheta-1
    do iphi=0,nphi-1
      temp1=temp1+(integrand(ir,itheta,iphi,dr,dtheta,dphi,k_vec)&
                  +integrand(ir,itheta,iphi+1,dr,dtheta,dphi,k_vec))*dphi/2.d0 
      
      temp2=temp2+(integrand(ir,itheta+1,iphi,dr,dtheta,dphi,k_vec)&
                  +integrand(ir,itheta+1,iphi+1,dr,dtheta,dphi,k_vec))*dphi/2.d0
    enddo
    temp3=temp3+(temp1+temp2)*dtheta/2.d0
    temp1=0.d0
    temp2=0.d0    
  
    do iphi=0,nphi-1
      temp1=temp1+(integrand(ir+1,itheta,iphi,dr,dtheta,dphi,k_vec)&
                  +integrand(ir+1,itheta,iphi+1,dr,dtheta,dphi,k_vec))*dphi/2.d0
      
      temp2=temp2+(integrand(ir+1,itheta+1,iphi,dr,dtheta,dphi,k_vec)&
                  +integrand(ir+1,itheta+1,iphi+1,dr,dtheta,dphi,k_vec))*dphi/2.d0
    enddo
    temp4=temp4+(temp1+temp2)*dtheta/2.d0
    temp1=0.d0
    temp2=0.d0    
  enddo
  phik=phik+(temp3+temp4)*dr/2.d0
  temp3=0.d0
  temp4=0.d0
enddo
end function phik

function integrand(ir,itheta,iphi,dr,dtheta,dphi,k_vec)
implicit none
real(8) :: k_vec(3)
integer :: ir,itheta,iphi
real(8) :: dr,dtheta,dphi
real(8) :: r,theta,phi

real(8) :: r_vec(3)
complex(8) :: phase_factor
complex(8) :: integrand
r=dr*ir
theta=dtheta*itheta
phi=dphi*iphi
r_vec(:)=(/ r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta)/)
phase_factor=cmplx(cos(-1.d0*(k_vec(1)*r_vec(1)+k_vec(2)*r_vec(2)+k_vec(3)*r_vec(3))),&
                   sin(-1.d0*(k_vec(1)*r_vec(1)+k_vec(2)*r_vec(2)+k_vec(3)*r_vec(3))),8)
integrand=pzorb(r,theta,phi)*phase_factor*(r**2)*sin(theta)
!integrand=pzorb(r,theta,phi)*phase_factor
end function integrand
function  pzorb(r,theta,phi)
implicit none
integer :: Z
real(8) :: r,theta,phi
real(8) :: pzorb
real(8) :: rfunc,yfunc
real(8) :: r_vec(3)
real(8) :: pi
pi=4.d0*ATAN(1.d0)

Z=6
rfunc=r*exp(-Z*r/2d0)*((sqrt(Z/2d0))**3)*Z/sqrt(3d0)
yfunc=cos(theta)*sqrt(3d0/4d0/pi)

pzorb=Rfunc*Yfunc
!pzorb=2.d0*sqrt(2.d0)/pi/(r**2+1.d0)**2/(sqrt(2.d0*pi))**3
end function pzorb

