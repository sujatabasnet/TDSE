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
complex(kind=8), dimension(257)::psi_t0, psi_t1, psi_t2, dd
real(kind=8), dimension(50, 257)::psi_square
integer::i, j

L=N*dx
!open(unit=10, file='100')

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



!calculating psi at t=1
   
   !calculating derivative at t=0
   call derv(psi_t0, dd)
   
   !Euler for t=1
   do i=1, N
   psi_t1(i)=psi_t0(i)+((iota*dt)*(((dd(i))/(2*m))-((1/h)*(arr_v(i))*(psi_t0(i)))))      
   enddo
   


!calculating the remaining psi's
   
   do j=3, t
   
      call derv(psi_t1, dd)
      
      !final formula
      do i=1, N
      psi_t2(i)=psi_t0(i)-(((2*iota*dt)/(h))*((dd(i)/(-2.0*m))+(psi_t1(i)*arr_v(i))))       
      enddo
      
      if(mod(j, 100) .eq. 0) then
      do i=1, N
        write(j, *) arr_x(i),  (abs(psi_t2(i)))**2
      enddo
      endif
      
      do i=1, N
      psi_t0(i)=psi_t1(i)
      psi_t1(i)=psi_t2(i)
      enddo
      
    enddo  
   
   
   
   
   
   
!subroutine for double derivative using fd   
   contains
   subroutine derv(psi, dd)
   use params
   complex(kind=8), dimension(N)::psi, dd
   integer::i, j

      dd(1)=(0,0)
      dd(N)=(0,0)
      do i=2, N-1       
        dd(i)=(1/(dx*dx))*(psi(i+1)-2*psi(i)+psi(i-1))
      enddo

   end subroutine
   
end program p   




