      integer, intent(in) :: ldx, ldy, ldz
      integer, intent(in) :: nx0, ny0, nz0
      integer, intent(in) :: ist, iend
      integer, intent(in) :: jst, jend
      double precision, intent(in) :: v(5,ldx/2*2+1,ldy/2*2+1,*)
      double precision, intent(out) :: summ(5)
