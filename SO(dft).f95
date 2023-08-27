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
complex(kind=8), dimension(257)::psi_t0, psi_t1,phi, rphi
integer::i, j

L=N*dx
open(unit=10, file='200')

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
do i=0,(N/2)-1
arr_k(i+1) = (2*pi*i)/(N*dx)
enddo
do i=N/2,N-1
arr_k(i+1) = (2*pi*(i-N))/(N*dx)
enddo



   


!calculating the remaining psi's
   
   do j=1, t
   
      !V/2
      do i=1, N
      psi_t1(i)=(exp((-iota*arr_v(i)*dt)/(2.0*h)))*psi_t0(i)
      enddo
   
      !FT at prev time step
      call ft(psi_t1, phi)
   
      do i=1, n
      phi(i)=phi(i)*(exp((-iota*dt*arr_k(i)*arr_k(i))/(2.0*m)))
      enddo
   
      !RFT 
      call rft(phi, rphi)
   
      !V/2
      do i=1, N
      psi_t1(i)=(exp((-iota*arr_v(i)*dt)/(2.0*h)))*rphi(i)
      enddo
      
      if(j .eq. 5000) then
      do i=1, N
        write(10, *) arr_x(i),  (abs(psi_t1(i)))**2
      enddo
      endif
      
      do i=1, N
      psi_t0(i)=psi_t1(i)
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

