subroutine r3cross(x,y,z)
implicit none
real(8), intent(in) :: x(3),y(3)
real(8), intent(out) :: z(3)
z(1)=x(2)*y(3)-x(3)*y(2)
z(2)=x(3)*y(1)-x(1)*y(3)
z(3)=x(1)*y(2)-x(2)*y(1)
return
end subroutine
