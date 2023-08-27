module params
real(kind=8)::alpha=20.0, x0=-0.5, p0=20.0, m=14500.0, xmin=-2.0, dx=0.02, dt=0.1, h=1.0, k0=20.0
integer(kind=8)::N=257, L, t=5000
complex(kind=8)::iota=(0, 1)
real(kind=8)::pi=acos(-1.0)
real(kind=8), dimension(257)::arr_x, arr_v, arr_k

end module params

program p
use params
implicit none
complex(kind=8), dimension(257)::psi_t0, psi_t1, psi_t2, phi, rphi
real(kind=8), dimension(50, 257)::psi_square
integer::i, j

L=N*dx
open(unit=10, file='100')

!storing x values
do i=1, N
arr_x(i)=xmin+(i-1)*dx
enddo


!storing v values
do i=1, N
if(arr_x(i) .lt. 0) then
arr_v(i)=0.0
else
arr_v(i)=1.0
endif
enddo

!calculating psi at t=0
do i=1, N
psi_t0(i)=(((2.0*alpha)/(pi))**0.25)*(exp(iota*k0*(arr_x(i)-x0)))*(exp(-alpha*((arr_x(i)-x0)**2)))
enddo


!storing k's
do i=0, N-1
arr_k(i+1)=(2.0*pi*(i-(N-1)/2))/(L)
enddo


!calculating psi at t=1
   
   !FT at t=0
   call ft(psi_t0, phi)
   
   !multiplying by -k**2
   do i=1, N
   phi(i)=(-(arr_k(i))**2)*phi(i)
   enddo
   
   !RFT at t=0
   call rft(phi, rphi)
   
   !Euler for t=1
   do i=1, N
   psi_t1(i)=psi_t0(i)+((iota*dt)*(((rphi(i))/(2*m))-((1/h)*(arr_v(i))*(psi_t0(i)))))      
   enddo
   
   
    do i=1, N
        ! write(10, *)  arr_x(i),  (abs(psi_t1(i)))**2
      enddo
   


!calculating the remaining psi's
   
   do j=2, t
   
      !FT at prev time step
      call ft(psi_t1, phi)
   
      !multiplying by -k**2
      do i=1, N
      phi(i)=(-(arr_k(i))**2)*phi(i)
      enddo
   
      !RFT 
      call rft(phi, rphi)
   
      !final formula
      do i=1, N
      psi_t2(i)=psi_t0(i)-(((2*iota*dt)/(h))*((rphi(i)/(-2.0*m))+(psi_t1(i)*arr_v(i))))       
      enddo
      
      if(j .eq. 4700) then
      do i=1, N
        write(10, *) arr_x(i),  (abs(psi_t2(i)))**2
      enddo
      endif
      
      do i=1, N
      psi_t0(i)=psi_t1(i)
      psi_t1(i)=psi_t2(i)
      enddo
      
    enddo  
   
   
   
   
   
   
!subroutine for fourier transformation   
   contains
   subroutine ft(psi, phi)
   use params
   complex(kind=8), dimension(N)::psi, phi
   integer::i, j

      do i=1, N
      phi(i)=(0,0)
         do j=1, N
         phi(i)=phi(i)+((psi(j))*(exp(-iota*arr_k(i)*arr_x(j))))
         enddo
      phi(i)=phi(i)/(sqrt(257.0))
      enddo

   end subroutine



!subroutine for reverse fourier transformation
   subroutine rft(phi, rphi)
   use params
   complex(kind=8), dimension(N)::phi, rphi
   integer::i, j
   
      do i=1, N
      rphi(i)=(0,0)
         do j=1, N
         rphi(i)=rphi(i)+((phi(j))*(exp(iota*arr_k(j)*arr_x(i))))
         enddo
         rphi(i)=rphi(i)/(sqrt(257.0))
      enddo
   end subroutine 
   

end program p
