c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine error_norm(rms)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     this function computes the norm of the difference between the
c     computed solution and the exact solution
c---------------------------------------------------------------------

      include 'header.h'
       interface
         pure extrinsic (hpf_local)
     >     subroutine exact_solution(xi,eta,zeta,dtemp)
           include 'cnst.h'
           double precision, intent(in) ::  xi, eta, zeta
           double precision, intent(out), dimension(:) :: dtemp
         end subroutine
       end interface
      double precision rmsa( 0:grid_points(3)-1, 5)
      integer i, j, k, m, d
      double precision xi, eta, zeta, u_exact(5), rms(5), add

      do m = 1, 5
         rms(m) = 0.0d0
      enddo

!HPF$ independent new(add,zeta,eta,xi)
      do k = 0, grid_points(3)-1
         do m = 1, 5
            rmsa(k, m) = 0.0d0
         enddo
         zeta = dble(k) * dnzm1
         do j = 0, grid_points(2)-1
            eta = dble(j) * dnym1
            do i = 0, grid_points(1)-1
               xi = dble(i) * dnxm1
               call exact_solution(xi, eta, zeta, u_exact)

               do m = 1, 5
                  add = u(m,i,j,k)-u_exact(m)
                  rmsa(k, m) = rmsa(k,m) + add*add
               enddo
            enddo
          enddo
       enddo

      do m = 1, 5
         rms(m) = sum( rmsa(:,m))
         do d = 1, 3
            rms(m) = rms(m) / dble(grid_points(d)-2)
         enddo
         rms(m) = dsqrt(rms(m))
      enddo

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine rhs_norm(rms)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      include 'header.h'

      integer d, m
      double precision rms(5), add

      do m = 1, 5
         rms(m) = sum( rhs(m,1:grid_points(1)-2,
     >                 1:grid_points(2)-2,1:grid_points(3)-2)**2)
      enddo

      do m = 1, 5
         do d = 1, 3
            rms(m) = rms(m) / dble(grid_points(d)-2)
         enddo
         rms(m) = dsqrt(rms(m))
      enddo

      return
      end

