subroutine gengrpvel(k_vec,grpvl_mat)
use modmain
use modconstants
implicit none
real(8), intent(in) :: k_vec(3)
real(8) :: dk, dkx_vec(3), dky_vec(3), kx_hat(3), ky_hat(3)
complex(8) :: partialxf_mat(norb,norb,3), partialyf_mat(norb,norb,3)
complex(8) :: temp1_mat(norb,norb), temp2_mat(norb,norb)
complex(8), intent(out) :: grpvl_mat(norb,norb,3)
integer :: i,j
dk=0.000001d0
dkx_vec=(/dk,0d0,0d0/)
dky_vec=(/0d0,dk,0d0/)
kx_hat=(/1d0,0d0,0d0/)
ky_hat=(/0d0,1d0,0d0/)
grpvl_mat(:,:,:)=0d0

call genhk(k_vec,temp1_mat)
call genhk(k_vec+dkx_vec,temp2_mat)
do i=1,norb
do j=1,norb
partialxf_mat(i,j,:)=(temp2_mat(i,j)-temp1_mat(i,j))/dk*kx_hat
enddo
enddo

call genhk(k_vec,temp1_mat)
call genhk(k_vec+dky_vec,temp2_mat)
do i=1,norb
do j=1,norb
partialyf_mat(i,j,:)=(temp2_mat(i,j)-temp1_mat(i,j))/dk*ky_hat
enddo
enddo

do i=1,norb
do j=1,norb
grpvl_mat(i,j,:)=partialxf_mat(i,j,:)+partialyf_mat(i,j,:)
enddo
enddo
contains
include 'deltak.f90'
end subroutine
