function cdot(vec1,vec2)
implicit none
real(8) cdot
real(8) vec1(3), vec2(3), phase
cdot=vec1(1)*vec2(1)+vec1(2)*vec2(2)+vec1(3)*vec2(3)
end function
