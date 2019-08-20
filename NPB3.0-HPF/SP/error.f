
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine error_norm(rms)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c this function computes the norm of the difference between the
c computed solution and the exact solution
c---------------------------------------------------------------------

       include 'header.h'

       integer i, j, k, m, d
       double precision xi, eta, zeta, u_exact(5), rms(5), add

       double precision rmsz(5,0:KMAX-1)
!hpf$ distribute(*,block) :: rmsz

       interface
         extrinsic (hpf_local) pure subroutine
     >                       exact_solution(xi,eta,zeta,dtemp)
           double precision, intent(in)::  xi, eta, zeta
           double precision, dimension (:), intent(out) :: dtemp
         end subroutine exact_solution
       end interface

!hpf$ independent, new(add, u_exact)
       do   k = 0, grid_points(3)-1
          do    m = 1, 5
             rmsz(m,k) = 0.0d0
          end do

          zeta = dble(k) * dnzm1
          do   j = 0, grid_points(2)-1
             eta = dble(j) * dnym1
             do   i = 0, grid_points(1)-1
                xi = dble(i) * dnxm1
                call exact_solution(xi, eta, zeta, u_exact(:))
                do   m = 1, 5
                   add = u(m,i,j,k)-u_exact(m)
                   rmsz(m,k) = rmsz(m,k) + add*add
                end do
             end do
          end do
       end do

       do    m = 1, 5
          rms(m) = sum( rmsz(m,0:grid_points(3)-1) )
          do    d = 1, 3
             rms(m) = rms(m) / dble(grid_points(d)-2)
          end do
          rms(m) = dsqrt(rms(m))
       end do

       return
       end



       subroutine rhs_norm(rms)

       include 'header.h'

       integer i, j, k, d, m
       double precision rms(5), add

       double precision rmsz(5,0:KMAX-1)
!hpf$ distribute(*,block) :: rmsz

!hpf$ independent, new(add)
       do k = 1, nz2
          do   m = 1, 5
             rmsz(m,k) = 0.0d0
          end do

          do j = 1, ny2
             do i = 1, nx2
       	       do m = 1, 5
                  add = rhs(m,i,j,k)
                  rmsz(m,k) = rmsz(m,k) + add*add
               end do 
             end do 
          end do 
       end do 

       do   m = 1, 5
          rms(m) = sum( rmsz(m,1:nz2) )
          do   d = 1, 3
             rms(m) = rms(m) / dble(grid_points(d)-2)
          end do
          rms(m) = dsqrt(rms(m))
       end do

       return
       end


