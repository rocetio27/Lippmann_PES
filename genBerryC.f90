subroutine genBerryC(eigvecs_temp,eigvals_temp,temp_grpvel_mat,berryC)
use modmain
implicit none
complex(8), intent(in) :: eigvecs_temp(norb,norb)
real(8), intent(in) :: eigvals_temp(norb)
complex(8), intent(in) :: temp_grpvel_mat(norb,norb,3)
real(8), intent(inout) :: berryC(norb)
complex(8) :: t_l,t_r
real(8) :: b
integer :: i,j,i1,i2
real(8) :: norm
berryC(:)=0.d0
do i=1,norb
do j=1,norb
  if(j.ne.i) then
  t_l=0.d0
  t_r=0.d0
  do i1=1,norb
    do i2=1,norb
      t_l=t_l+conjg(eigvecs_temp(i1,i))*temp_grpvel_mat(i1,i2,1)*eigvecs_temp(i2,j)
      t_r=t_r+conjg(eigvecs_temp(i1,j))*temp_grpvel_mat(i1,i2,2)*eigvecs_temp(i2,i)
    enddo
  enddo
  b=eigvals_temp(i)-eigvals_temp(j)
  b=b**2
  berryC(i)=berryC(i)-2*aimag(t_l*t_r/b)
  endif
enddo
enddo
end subroutine genBerryC
