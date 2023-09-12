subroutine rk4n(n,hk_t_1,hk_t_2,hk_t_3,dt,xi)

implicit none
integer, intent(in) :: n
complex(8),intent(in) :: hk_t_1(n,n)
complex(8),intent(in) :: hk_t_2(n,n)
complex(8),intent(in) :: hk_t_3(n,n)
real(8), intent(in) :: dt
complex(8), intent(inout) :: xi(n)

integer j
real(8) h ! t : start time t+h : evolved time
complex(8) x(n), dxdt(n)
complex(8) xf(n)
complex(8) k1(n),k2(n),k3(n),k4(n)

h = dt

!* evaluate k1
call fcn(n, hk_t_1, xi, dxdt)
do j=1,n
   k1(j) = h*dxdt(j)
   x(j)  = xi(j) + k1(j)/2.d0  
end do      

!* evaluate k2
call fcn(n, hk_t_2, x, dxdt)
do j=1,n
   k2(j) = h*dxdt(j) 
   x(j)  = xi(j) + k2(j)/2.d0  
end do

!* evaluate k3
call fcn(n, hk_t_2, x, dxdt)
do j=1,n
   k3(j) = h*dxdt(j) 
   x(j)  = xi(j) + k3(j)   
end do     

!* evaluate k4 and the result      
call fcn(n, hk_t_3, x, dxdt) 
do j=1,n
   k4(j) = h*dxdt(j) 
   xf(j) = xi(j) + k1(j)/6.d0+k2(j)/3.d0+k3(j)/3.d0+k4(j)/6.d0
end do     
xi=xf
end subroutine rk4n
