subroutine fcn(n,hk_t,x,dxdt)
   use modpulse
   implicit none
   integer :: i,j
   integer, intent(in) :: n
   complex(8), intent(in) :: hk_t(n,n)
   complex(8),intent(in) :: x(n)
   complex(8),intent(out) :: dxdt(n)
   dxdt=0.d0
   do i=1,n
     do j=1,n
       dxdt(i)=dxdt(i)+hk_t(i,j)*x(j)*cmplx(0.d0,-1.d0)
     enddo
   enddo
   contains
end subroutine
