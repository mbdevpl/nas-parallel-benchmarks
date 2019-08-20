c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine  initialize
       include 'header.h'
       interface
         pure extrinsic (hpf_local)
     >     subroutine exact_solution(xi,eta,zeta,dtemp)
           include 'cnst.h'
           double precision, intent(in) ::  xi, eta, zeta
           double precision, intent(out), dimension(:) :: dtemp
         end subroutine
       end interface

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     This subroutine initializes the field variable u using
c     tri-linear transfinite interpolation of the boundary values
c---------------------------------------------------------------------

      integer i, j, k, m, ix, iy, iz
      double precision  xi, eta, zeta, Pface(5,3,2), Pxi, Peta,
     >     Pzeta, temp(5)

c---------------------------------------------------------------------
c  Later (in compute_rhs) we compute 1/u for every element. A few of
c  the corner elements are not used, but it convenient (and faster)
c  to compute the whole thing with a simple loop. Make sure those
c  values are nonzero by initializing the whole thing here.
c---------------------------------------------------------------------
!HPF$ independent
      do k = 0, grid_points(3)-1
         do j = 0, grid_points(2)-1
            do i = 0, grid_points(1)-1
      	       do m = 1, 5
                  u(m,i,j,k) = 1.0
               end do
            end do
         end do
      end do
c---------------------------------------------------------------------



c---------------------------------------------------------------------
c     first store the "interpolated" values everywhere on the grid
c---------------------------------------------------------------------
!HPF$ independent  new(zeta, eta, xi, Pface, Pxi, Peta, Pzeta)
      do k = 0, grid_points(3)-1
         zeta = dble(k) * dnzm1
         do j = 0, grid_points(2)-1
            eta = dble(j) * dnym1
            do i = 0, grid_points(1)-1
               xi = dble(i) * dnxm1

               do ix = 1, 2
                  call exact_solution(dble(ix-1), eta, zeta,
     >                    Pface(:,1,ix))
               enddo

               do iy = 1, 2
                  call exact_solution(xi, dble(iy-1) , zeta,
     >                    Pface(:,2,iy))
               enddo

               do iz = 1, 2
                  call exact_solution(xi, eta, dble(iz-1),
     >                    Pface(:,3,iz))
               enddo

               do m = 1, 5
                  Pxi   = xi   * Pface(m,1,2) +
     >                    (1.0d0-xi)   * Pface(m,1,1)
                  Peta  = eta  * Pface(m,2,2) +
     >                    (1.0d0-eta)  * Pface(m,2,1)
                  Pzeta = zeta * Pface(m,3,2) +
     >                    (1.0d0-zeta) * Pface(m,3,1)

                  u(m,i,j,k) = Pxi + Peta + Pzeta -
     >                    Pxi*Peta - Pxi*Pzeta - Peta*Pzeta +
     >                    Pxi*Peta*Pzeta

               enddo
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     now store the exact values on the boundaries
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     west face
c---------------------------------------------------------------------
      i = 0
      xi = 0.0d0
      do k = 0, grid_points(3)-1
         zeta = dble(k) * dnzm1
         do j = 0, grid_points(2)-1
            eta = dble(j) * dnym1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     east face
c---------------------------------------------------------------------

      i = grid_points(1)-1
      xi = 1.0d0
      do k = 0, grid_points(3)-1
         zeta = dble(k) * dnzm1
         do j = 0, grid_points(2)-1
            eta = dble(j) * dnym1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     south face
c---------------------------------------------------------------------
      j = 0
      eta = 0.0d0
      do k = 0, grid_points(3)-1
         zeta = dble(k) * dnzm1
         do i = 0, grid_points(1)-1
            xi = dble(i) * dnxm1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo


c---------------------------------------------------------------------
c     north face
c---------------------------------------------------------------------
      j = grid_points(2)-1
      eta = 1.0d0
      do k = 0, grid_points(3)-1
         zeta = dble(k) * dnzm1
         do i = 0, grid_points(1)-1
            xi = dble(i) * dnxm1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     bottom face
c---------------------------------------------------------------------
      k = 0
      zeta = 0.0d0
      do i =0, grid_points(1)-1
         xi = dble(i) *dnxm1
         do j = 0, grid_points(2)-1
            eta = dble(j) * dnym1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     top face
c---------------------------------------------------------------------
      k = grid_points(3)-1
      zeta = 1.0d0
      do i =0, grid_points(1)-1
         xi = dble(i) * dnxm1
         do j = 0, grid_points(2)-1
            eta = dble(j) * dnym1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine lhsinit(lhs, size)
      implicit none
      integer size
      double precision lhs(5,5,3,0:size)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      integer i, m, n

c---------------------------------------------------------------------
c     zero the whole left hand side for starters
c---------------------------------------------------------------------
      do i = 0, size, size
         do m = 1, 5
            do n = 1, 5
               lhs(m,n,1,i) = 0.0d0
               lhs(m,n,2,i) = 0.0d0
               lhs(m,n,3,i) = 0.0d0
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     next, set all diagonal values to 1. This is overkill, but convenient
c---------------------------------------------------------------------
      do i = 0, size, size
         do m = 1, 5
            lhs(m,m,2,i) = 1.0d0
         enddo
      enddo

      return
      end



