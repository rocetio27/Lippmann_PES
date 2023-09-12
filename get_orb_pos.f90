subroutine get_orb_pos(tmp_orb_pos)
use modconstants
use modmain
implicit none
real(8), intent(inout) :: tmp_orb_pos(norb,3)

tmp_orb_pos(1,:)=tau_1_vec
tmp_orb_pos(2,:)=0d0

!tmp_orb_pos(3,:)=tau_1_vec-(/0d0,0d0,d_0/)
!tmp_orb_pos(4,:)=-(/0d0,0d0,d_0/)


end subroutine
