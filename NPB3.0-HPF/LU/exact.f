c---------------------------------------------------------------------
c---------------------------------------------------------------------
       pure extrinsic (hpf_local) subroutine exact(i,j,k,u000ijk)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c
c   compute the exact solution at (i,j,k)
c
c---------------------------------------------------------------------

      implicit none
c---------------------------------------------------------------------
c  input parameters
c---------------------------------------------------------------------
      include 'exact.h'

c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------
      double precision ce(5,13)
      common/cexact/ ce
      integer nx, ny, nz
      integer nx0, ny0, nz0
      integer ist, iend
      integer jst, jend
      integer ii1, ii2
      integer ji1, ji2
      integer ki1, ki2
      double precision  dxi, deta, dzeta
      double precision  tx1, tx2, tx3
      double precision  ty1, ty2, ty3
      double precision  tz1, tz2, tz3

      common/cgcon/ dxi, deta, dzeta,
     >              tx1, tx2, tx3,
     >              ty1, ty2, ty3,
     >              tz1, tz2, tz3,
     >              nx, ny, nz,
     >              nx0, ny0, nz0,
     >              ist, iend,
     >              jst, jend,
     >              ii1, ii2,
     >              ji1, ji2,
     >              ki1, ki2
      integer m
      double precision xi, eta, zeta

      xi  = ( dble ( i - 1 ) ) / ( nx0 - 1 )
      eta  = ( dble ( j - 1 ) ) / ( ny0 - 1 )
      zeta = ( dble ( k - 1 ) ) / ( nz - 1 )


      do m = 1, 5
         u000ijk(m) =  ce(m,1)
     >        + (ce(m,2)
     >        + (ce(m,5)
     >        + (ce(m,8)
     >        +  ce(m,11) * xi) * xi) * xi) * xi
     >        + (ce(m,3)
     >        + (ce(m,6)
     >        + (ce(m,9)
     >        +  ce(m,12) * eta) * eta) * eta) * eta
     >        + (ce(m,4)
     >        + (ce(m,7)
     >        + (ce(m,10)
     >        +  ce(m,13) * zeta) * zeta) * zeta) * zeta
      end do

      return
      end
