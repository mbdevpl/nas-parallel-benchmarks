c---------------------------------------------------------------------
c compute the roots-of-unity array that will be used for subsequent FFTs.
c---------------------------------------------------------------------
      PURE extrinsic (hpf_local) subroutine CompExp (n, exponent)

      implicit none
      integer, intent(in) :: n
      double complex, dimension(:), intent(inout) :: exponent

      interface
        pure extrinsic (hpf_local) integer function ilog2 (n)
          integer, intent(in) :: n
        end function ilog2
      end interface
      integer m,nu,ku,i,j,ln
      double precision t, ti, pi
      data pi /3.141592653589793238d0/

      nu = n
      m = ilog2(n)
      exponent(1) = m
      ku = 2
      ln = 1
      do j = 1, m
         t = pi / ln
         do i = 0, ln - 1
            ti = i * t
            exponent(i+ku) = dcmplx(cos(ti),sin(ti))
         enddo
         ku = ku + ln
         ln = 2 * ln
      enddo

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      pure extrinsic (hpf_local) integer function ilog2(n)
      implicit none
      integer, intent( in) :: n


c---------------------------------------------------------------------
c---------------------------------------------------------------------


      integer nn, lg
      if (n .eq. 1) then
         ilog2=0
         return
      endif
      lg = 1
      nn = 2
      do while (nn .lt. n)
         nn = nn*2
         lg = lg+1
      end do
      ilog2 = lg
      return
      end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
       extrinsic (HPF)
     * subroutine compute_initial_conditions(u0,d1,d2,d3)

      double complex, dimension(:,:,:), intent(inout) :: u0
      integer, intent(in) :: d1,d2,d3
      double precision x0, start, an, dummy
      double precision, dimension(d3) :: RanStarts
!hpf$    template grid(d3)
!hpf$    distribute(block) :: grid
!hpf$    align(*,*,:) with grid :: u0
!hpf$    align(:) with grid :: RanStarts
      interface
        pure extrinsic(HPF) SUBROUTINE VRanComp (N, X, A, Y,m)
          implicit none
          integer, intent(in) :: n,m
          double precision, intent(inout) :: x, a
          double complex, dimension(N), intent(out) :: y
	end SUBROUTINE VRanComp
      end interface
      integer i,j,k
      double precision seed, a
      parameter (seed = 314159265.d0, a = 1220703125.d0)
      external randlc
      double precision randlc

      start = seed
c---------------------------------------------------------------------
c Jump to the starting element for our first plane.
c---------------------------------------------------------------------
      call ipow46(a, 0, an)
      dummy = randlc(start, an)
      call ipow46(a, 2*d1*d2, an)
c---------------------------------------------------------------------
c Go through by z planes filling in one square at a time.
c---------------------------------------------------------------------
      RanStarts(1) = start
      do k = 2, d3
         dummy = randlc(start, an)
         RanStarts(k) = start
      end do

!HPF$ independent
      do k = 1, d3
         do j = 1, d1
           call VRanComp(d1,RanStarts(k),a,u0(j,:,k),d2)
         end do
         if (k .ne. d3) dummy = randlc(start, an)
      end do
      return
      end
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      pure extrinsic(HPF) SUBROUTINE VRanComp(N, X, A, Y,m)
      implicit none

      integer, intent(in) :: n,m
      double precision, intent(inout) :: x, a
      double complex, dimension(m), intent(out) :: y

      integer*8 i246m1, Lx, La, i
      double precision d2m46, Re, Im

      parameter(d2m46=0.5d0**46)

      data i246m1 /0/

      if ( i246m1.EQ.0 ) then
           i246m1 = '00003FFFFFFFFFFF'x
      endif

      Lx = X
      La = A
      do i = 1, m
         Lx   = iand(Lx*La,i246m1)
         Re = d2m46*dble(Lx)
         Lx   = iand(Lx*La,i246m1)
         Im = d2m46*dble(Lx)
         y(i) = dcmplx ( Re, Im)
      end do
      x = dble(Lx)

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine ipow46(a, exponent, result)
c---------------------------------------------------------------------
c compute a^exponent mod 2^46
c---------------------------------------------------------------------

      implicit none
      double precision a, result, dummy, q, r
      integer exponent, n, n2
      external randlc
      double precision randlc
c---------------------------------------------------------------------
c Use
c   a^n = a^(n/2)*a^(n/2) if n even else
c   a^n = a*a^(n-1)       if n odd
c---------------------------------------------------------------------
      result = 1
      if (exponent .eq. 0) return
      q = a
      r = 1
      n = exponent

      do while (n .gt. 1)
         n2 = n/2
         if (n2 * 2 .eq. n) then
            dummy = randlc(q, q)
            n = n2
         else
            dummy = randlc(r, q)
            n = n-1
         endif
      end do
      dummy = randlc(r, q)
      result = r
      return
      end

c---------------------------------------------------------------------
      extrinsic (HPF)
     >subroutine CalculateChecksum(csum,iterN,u,d1,d2,d3)
        implicit none
        integer, intent (in) :: iterN
        double complex, dimension(:,:,:), intent(in) :: u
        integer, intent (in) :: d1, d2, d3
        double complex, intent (out) :: csum
!hpf$   template grid(d3)
!hpf$   distribute(block) :: grid
!hpf$   align(*,*,:) with grid :: u
      integer i, i1, ii, ji, ki
      integer dl1,dl2,dl3
      dl1=d1-1
      dl2=d2
      dl3=d3
      csum = dcmplx (0.0, 0.0)
      do i = 1, 1024
        i1 = i - 1
        ii = mod (i1, dl2) + 1
        ji = mod (3 * i1, dl1) + 1
        ki = mod (5 * i1, dl3) + 1
        csum = csum + u(ji,ii,ki)
      end do
      csum = csum/dble(dl1*dl2*dl3)
      print*,'T =',iterN,' checksum =',csum
      return
      end

