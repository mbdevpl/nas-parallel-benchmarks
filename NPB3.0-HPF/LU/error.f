c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine error

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c
c   compute the solution error
c
c---------------------------------------------------------------------

      implicit none

      include 'applu.incl'
      interface
            pure extrinsic (hpf_local) subroutine exact(i,j,k,u000ijk)
              implicit none
              include 'exact.h'
            end subroutine exact
      end interface

c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------
      integer i, j, k, m
      double precision  tmp
      double precision  u000ijk(5)



      do m = 1, 5
         errnm(m) = 0.0d+00
      end do

!HPF$ INDEPENDENT, new(u000ijk)
         do j = jst, jend
      do k = 2, nz-1
            do i = ist, iend
               call exact( i, j, k, u000ijk )
               do m = 1, 5
                  rsd(m,i,j,k) = u(m,i,j,k) - u000ijk(m)
               end do
            end do
         end do
      end do
      do m = 1, 5
         errnm(m) = sum(rsd(m,ist:iend,jst:jend,2:nz-1)**2)
      end do

      do m = 1, 5
         errnm(m) = sqrt ( errnm(m) / ( (nx0-2)*(ny0-2)*(nz0-2) ) )
      end do

c        write (*,1002) ( errnm(m), m = 1, 5 )

 1002 format (1x/1x,'RMS-norm of error in soln. to ',
     > 'first pde  = ',1pe12.5/,
     > 1x,'RMS-norm of error in soln. to ',
     > 'second pde = ',1pe12.5/,
     > 1x,'RMS-norm of error in soln. to ',
     > 'third pde  = ',1pe12.5/,
     > 1x,'RMS-norm of error in soln. to ',
     > 'fourth pde = ',1pe12.5/,
     > 1x,'RMS-norm of error in soln. to ',
     > 'fifth pde  = ',1pe12.5)

      return
      end
