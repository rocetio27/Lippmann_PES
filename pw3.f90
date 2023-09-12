function pw3(k_vec,r_vec)
implicit none
complex(8) pw3
real(8) k_vec(3), r_vec(3), phase
phase=k_vec(1)*r_vec(1)+k_vec(2)*r_vec(2)+k_vec(3)*r_vec(3)
pw3=cmplx(dcos(phase),dsin(phase),8)
end function
