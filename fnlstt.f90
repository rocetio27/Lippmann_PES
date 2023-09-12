function  fnlstt(r_vec,p_vec,j)
use modconstants
use modmain
implicit none
complex(8) fnlstt(3)
real(8) r_vec(3), z, rho_vec(3)
real(8) p_vec(3), pz, pp_vec(3)
real(8) g_vec(3), G, Gm
real(8) pzG, pzG_vr(3)
complex(8) pzG_vc(3),temp_vec1(3), temp_vec2(3)
real(8) tau_i_vec(3), tau_iz
real(8) tau_j_vec(3), tau_jz
real(8) tau_d_vec(3), tau_dz
complex(8) vc(3), vzhc(3), phasef1, phasef2
real(8) vr(3), vzhr(3), phase1, phase2

integer :: i,j,ig
real(8) :: a_unit
complex(8) :: coeff
a_unit=(2.46d0*astr)**2/2d0
coeff=4d0*pi*ci/a_unit
z=r_vec(3)
pz=p_vec(3)
rho_vec=(/r_vec(1), r_vec(2), 0d0/)
pp_vec=(/p_vec(1), p_vec(2), 0d0/)
tau_j_vec=orb_pos(j,:) !===
tau_jz=orb_pos(j,3)    !===
! 1. Plane wave
fnlstt=pw3(-p_vec,r_vec)*p_vec
vc=(/c0,c0,c0/)
vzhc=(/c0,c0,c1/)
vzhr=(/0d0,0d0,1d0/)
! 2. p_z^G wave
if(1.eq.1) then
do i=1,norb ! The is is index of the final state

tau_i_vec=orb_pos(i,:) !===
tau_iz=orb_pos(i,3)    !===

tau_d_vec=tau_i_vec-tau_j_vec
tau_dz=tau_iz-tau_jz

do ig=1,ngtot
g_vec=g_vecs(ig,:)
vr=pp_vec+g_vec
vc=pp_vec+g_vec
G=norm2(g_vec)
Gm=dsqrt(G**2+mu**2)
!if (ig.eq.((ngtot-1)/2+1)) then 
  ! NOTHING
!else

  if((norm2(p_vec)**2-norm2(pp_vec+g_vec)**2).ge.0d0) then
  pzG=dsqrt(norm2(p_vec)**2-norm2(pp_vec+g_vec)**2)
  pzG_vr=pzG*vzhr
  pzG_vc=pzG*vzhc
    if(z.ge.tau_dz)then
    phase1=cdot(-vr,rho_vec)+cdot(g_vec-(/0d0,0d0,pz/),tau_d_vec)-pzG*(z-tau_dz)
    phase2=cdot(-vr,rho_vec)+cdot(g_vec-(/0d0,0d0,pz/),tau_d_vec) -pz*(z-tau_dz)
    phasef1=dcos(phase1)+ci*dsin(phase1)
    phasef2=dcos(phase2)+ci*dsin(phase2)
    fnlstt=fnlstt+coeff*atnum(i)&
                   *(&
                   -phasef1*(vr+pzG_vr)/(pzG*(pz-pzG+ci*Gm)*(pz-pzG-ci*Gm))&
                   +phasef2*dexp(-Gm*(z-tau_dz))*(vc+(pz-ci*Gm)*vzhc)/(ci*Gm*(pz-pzG-ci*Gm)*(pz+pzG-ci*Gm))&
                    )
    else
    phase1=cdot(-vr,rho_vec)+cdot(g_vec-(/0d0,0d0,pz/),tau_d_vec)+pzG*(z-tau_dz)
    phase2=cdot(-vr,rho_vec)+cdot(g_vec-(/0d0,0d0,pz/),tau_d_vec) -pz*(z-tau_dz)
    phasef1=dcos(phase1)+ci*dsin(phase1)
    phasef2=dcos(phase2)+ci*dsin(phase2)
    fnlstt=fnlstt+coeff*atnum(i)&
                   *(&
                   -phasef1*(vr-pzG_vr)/(pzG*(pz+pzG+ci*Gm)*(pz+pzG-ci*Gm))&
                   +phasef2*dexp(Gm*(z-tau_dz))*(vc+(pz+ci*Gm)*vzhc)/(ci*Gm*(pz-pzG+ci*Gm)*(pz+pzG+ci*Gm))&
                    )

    endif
  else
  pzG=dsqrt(abs(norm2(p_vec)**2-norm2(pp_vec+g_vec)**2))
  pzG_vr=pzG*vzhr
  pzG_vc=pzG*vzhc

    if(z.ge.tau_dz)then
    phase1=cdot(-vr,rho_vec)+cdot(g_vec-(/0d0,0d0,pz/),tau_d_vec)
    phase2=cdot(-vr,rho_vec)+cdot(g_vec-(/0d0,0d0,pz/),tau_d_vec)-pz*(z-tau_dz)
    phasef1=dcos(phase1)+ci*dsin(phase1)
    phasef2=dcos(phase2)+ci*dsin(phase2)
    fnlstt=fnlstt+coeff*atnum(i)&
                   *(&
                    phasef1*dexp(-pzG*(z-tau_dz))*(vc-ci*pzG_vc)/(ci*pzG*(pz+ci*pzG-ci*Gm)*(pz+ci*pzG+ci*Gm))&
                   +phasef2*dexp(-Gm*(z-tau_dz))*(vc+(pz-ci*Gm)*vzhc)/(ci*Gm*(pz+ci*pzG-ci*Gm)*(pz-ci*pzG-ci*Gm))&
                    )
   else
   phase1=cdot(-vr,rho_vec)+cdot(g_vec-(/0d0,0d0,pz/),tau_d_vec)
   phase2=cdot(-vr,rho_vec)+cdot(g_vec-(/0d0,0d0,pz/),tau_d_vec)-pz*(z-tau_dz)
   phasef1=dcos(phase1)+ci*dsin(phase1)
   phasef2=dcos(phase2)+ci*dsin(phase2)
   fnlstt=fnlstt+coeff*atnum(i)&
                   *(&
                    phasef1*dexp(pzG*(z-tau_dz))*(vc+ci*pzG_vc)/(ci*pzG*(pz-ci*pzG-ci*Gm)*(pz-ci*pzG+ci*Gm))&
                   +phasef2*dexp(Gm*(z-tau_dz))*(vc+(pz+ci*Gm)*vzhc)/(ci*Gm*(pz+ci*pzG+ci*Gm)*(pz-ci*pzG+ci*Gm))&
                    )
   endif
 endif
!endif
enddo
enddo
endif
end function fnlstt
