


module m
implicit none
real(kind = 8), parameter :: dx = 0.01, pi = acos(-1.0), dt=0.0001
complex(kind = 8), parameter :: iota = (0.0, 1.0)
end module m


program q2
use m
implicit none
real(kind = 8), dimension(1000) :: x
complex(kind = 8), dimension(1000)  :: s1, s2, s3, s4, temp
complex(kind = 8) :: f
complex(kind = 8), dimension(1000, 300) :: psi
integer :: i, j


do i = 1, 1000
 x(i) = dx*(i-1)
end do

do j = 1, 1000
 psi(j, 1) = ((40.0/pi)**(0.25)) * exp(-1.0* ((x(j)-5.0)**2) / (2.0* (0.1**2)) ) * exp(iota*20.0*x(j))
end do


do i=2, 300
  do j = 2, 999
   temp(j) = psi(j, i-1)
  end do
  
  do j = 2, 999
   s1(j) = f( j, temp)
  end do
  
  do j = 2, 999
   temp(j) = psi(j, i-1) + (dt*s1(j)/2.0)
  end do
  
  do j = 2, 999
   s2(j) = f( j, temp)
  end do
  
  do j = 2, 999
   temp(j) = psi(j, i-1) + (dt*s2(j)/2.0)
  end do
  
  do j = 2, 999
   s3(j) = f( j, temp)
  end do
  
  do j = 2, 999
   temp(j) = psi(j, i-1) + (dt*s3(j))
  end do
  
  do j = 2, 999
   s4(j) = f(j, temp )
  end do
   
  do j = 2, 999
  psi(j, i) = psi( j, i-1) + dt*(s1(j) + 2.0*s2(j) + 2.0*s3(j) + s4(j))/6.0
  end do
 
end do

do j = 1, 1000
  write(10,*) x(j), abs(psi(j, 1))**2
  write(20,*) x(j), abs(psi(j, 150))**2
  write(30,*) x(j), abs(psi(j, 300))**2
end do


end program

complex(kind = 8) function f(j,temp)
 use m
 implicit none
 integer :: j
 complex(kind = 8), dimension(1000) :: temp
 f = iota*( (temp(j+1) + temp(j-1) - 2.0*(temp(j)) )/ (2.0*(dx**2)) )
end function f



