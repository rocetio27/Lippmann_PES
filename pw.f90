function pw(k,r)
implicit none
complex(8) pw
real(8) k, r, phase
phase=k*r
pw=cmplx(dcos(phase),dsin(phase),8)
end function
