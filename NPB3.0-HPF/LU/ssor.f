c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine ssor

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c   to perform pseudo-time stepping SSOR iterations
c   for five nonlinear pde's.
c---------------------------------------------------------------------

      implicit none
      include 'applu.incl'
      interface       
            pure extrinsic (hpf) 
     >        subroutine l2norm ( ldx, ldy, ldz, 
     >                    nx0, ny0, nz0,
     >                    ist, iend, 
     >                    jst, jend,
     >                    v, summ )

              implicit none
              include 'l2norm.h'
            end subroutine l2norm
                
            extrinsic (hpf) 
     >        subroutine blts ( ldmx, ldmy, ldmz,
     >                  l,
     >                  omega,
     >                  v,
     >                  ldz, ldy, ldx, d)
              implicit none
              include 'blts.h'
            end subroutine blts
            
            extrinsic (hpf) 
     >        subroutine buts ( ldmx, ldmy, ldmz,
     >                  l,
     >                  omega,
     >                  v,
     >                  d, udz, udy, udx)     
              implicit none
              include 'buts.h'
            end subroutine buts                       
      end interface
c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------
      integer i, j, k, m, n
      integer istep
      integer l, lst, lend, np(isiz1+isiz2+isiz3-3)
      integer indxp(isiz1*isiz2*3/4,isiz1+isiz2+isiz3-3),
     >        jndxp(isiz1*isiz2*3/4,isiz1+isiz2+isiz3-3)
      double precision  tmp
      double precision  delunm(5), tv(5,isiz1/2*2+1,isiz2)

      external timer_read
      double precision timer_read


 
c---------------------------------------------------------------------
c   begin pseudo-time stepping iterations
c---------------------------------------------------------------------
      tmp = 1.0d+00 / ( omega * ( 2.0d+00 - omega ) ) 
      lst = ist + jst + 2
      lend = iend + jend + nz - 1
c---------------------------------------------------------------------
c   initialize a,b,c,d to zero (guarantees that page tables have been
c   formed, if applicable on given architecture, before timestepping).
c---------------------------------------------------------------------
      a = 0.0
      b = 0.0
      c = 0.0
      d = 0.0
      do i = 1, t_last
      	 call timer_clear(i)
      end do

c---------------------------------------------------------------------
c   compute the steady-state residuals
c---------------------------------------------------------------------
      call rhs
 
c---------------------------------------------------------------------
c   compute the L2 norms of newton iteration residuals
c---------------------------------------------------------------------
      call l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0,
     >             ist, iend, jst, jend,
     >             rsd, rsdnm )


c      if ( ipr .eq. 1 ) then
c         write (*,*) '          Initial residual norms'
c         write (*,*)
c         write (*,1007) ( rsdnm(m), m = 1, 5 )
c	 write (*,'(/a)') 'Iteration RMS-residual of 5th PDE'
c      end if
 
 
      do i = 1, t_last
      	 call timer_clear(i)
      end do
      call timer_start(t_total)
      call InitSectionsN()
      call SetLoopBounds()
c      print *, 'Actual Iterations', itmax
 
c---------------------------------------------------------------------
c   the timestep loop
c---------------------------------------------------------------------
      do istep = 1, itmax
         
c         if ( ( mod ( istep, inorm ) .eq. 0 ) .and.
c     >          ipr .eq. 1 ) then
c             write ( *, 1001 ) istep
c         end if
         if (mod ( istep, 20) .eq. 0 .or.
     >         istep .eq. itmax .or.
     >         istep .eq. 1) then
            write( *, 200) istep
 200        format(' Time step ', i4)
         endif
 
c---------------------------------------------------------------------
c   perform SSOR iteration
c---------------------------------------------------------------------
	 if (timeron) call timer_start(t_rhs)
         rsd = dt*rsd
	 if (timeron) call timer_stop(t_rhs)
	 do l = lst, lend
c---------------------------------------------------------------------
c   form the lower triangular part of the jacobian matrix
c---------------------------------------------------------------------
	    if (timeron) call timer_start(t_jacld)
            call jacld(l) 
 	    if (timeron) call timer_stop(t_jacld)
 
c---------------------------------------------------------------------
c   perform the lower triangular solution
c---------------------------------------------------------------------
	    if (timeron) call timer_start(t_blts)
            call blts( isiz1, isiz2, isiz3,
     >                  l,
     >                  omega,
     >                  rsd,
     >                  a, b, c, d)
	    if (timeron) call timer_stop(t_blts)
	 end do
 
	 do l = lend, lst, -1
c---------------------------------------------------------------------
c   form the strictly upper triangular part of the jacobian matrix
c---------------------------------------------------------------------
	    if (timeron) call timer_start(t_jacu)
            call jacu(l)
	    if (timeron) call timer_stop(t_jacu)

c---------------------------------------------------------------------
c   perform the upper triangular solution
c---------------------------------------------------------------------
	    if (timeron) call timer_start(t_buts)
            call buts( isiz1, isiz2, isiz3,
     >                  l,
     >                  omega,
     >                  rsd, 
     >                  d, a, b, c)
	    if (timeron) call timer_stop(t_buts)
	 end do
 
c---------------------------------------------------------------------
c   update the variables
c---------------------------------------------------------------------

	 if (timeron) call timer_start(t_add)
!HPF$ INDEPENDENT new(tmp)
            do j = jst, jend
         do k = 2, nz-1
               do i = ist, iend
                  do m = 1, 5
                     u( m, i, j, k ) = u( m, i, j, k )
     >                    + tmp * rsd( m, i, j, k )
                  end do
               end do
            end do
         end do
	 if (timeron) call timer_stop(t_add)
 
c---------------------------------------------------------------------
c   compute the max-norms of newton iteration corrections
c---------------------------------------------------------------------
         if ( mod ( istep, inorm ) .eq. 0 ) then
	    if (timeron) call timer_start(t_l2norm)
            call l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0,
     >                   ist, iend, jst, jend,
     >                   rsd, delunm )
	    if (timeron) call timer_stop(t_l2norm)
c            if ( ipr .eq. 1 ) then
c                write (*,1006) ( delunm(m), m = 1, 5 )
c            else if ( ipr .eq. 2 ) then
c                write (*,'(i5,f15.6)') istep,delunm(5)
c            end if
         end if
 
c---------------------------------------------------------------------
c   compute the steady-state residuals
c---------------------------------------------------------------------
	 if (timeron) call timer_start(t_rhs)
         call rhs
	 if (timeron) call timer_stop(t_rhs)
 
c---------------------------------------------------------------------
c   compute the max-norms of newton iteration residuals
c---------------------------------------------------------------------
         if ( ( mod ( istep, inorm ) .eq. 0 ) .or.
     >        ( istep .eq. itmax ) ) then
	    if (timeron) call timer_start(t_l2norm)
            call l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0,
     >                   ist, iend, jst, jend,
     >                   rsd, rsdnm )
	    if (timeron) call timer_stop(t_l2norm)
c            if ( ipr .eq. 1 ) then
c                write (*,1007) ( rsdnm(m), m = 1, 5 )
c            end if
         end if

c---------------------------------------------------------------------
c   check the newton-iteration residuals against the tolerance levels
c---------------------------------------------------------------------
         if ( ( rsdnm(1) .lt. tolrsd(1) ) .and.
     >        ( rsdnm(2) .lt. tolrsd(2) ) .and.
     >        ( rsdnm(3) .lt. tolrsd(3) ) .and.
     >        ( rsdnm(4) .lt. tolrsd(4) ) .and.
     >        ( rsdnm(5) .lt. tolrsd(5) ) ) then
c            if (ipr .eq. 1 ) then
c               write (*,1004) istep
c            end if
            return
         end if
 
      end do
 
      call timer_stop(t_total)
      maxtime= timer_read(t_total)
      return
      
 1001 format (1x/5x,'pseudo-time SSOR iteration no.=',i4/)
 1004 format (1x/1x,'convergence was achieved after ',i4,
     >   ' pseudo-time steps' )
 1006 format (1x/1x,'RMS-norm of SSOR-iteration correction ',
     > 'for first pde  = ',1pe12.5/,
     > 1x,'RMS-norm of SSOR-iteration correction ',
     > 'for second pde = ',1pe12.5/,
     > 1x,'RMS-norm of SSOR-iteration correction ',
     > 'for third pde  = ',1pe12.5/,
     > 1x,'RMS-norm of SSOR-iteration correction ',
     > 'for fourth pde = ',1pe12.5/,
     > 1x,'RMS-norm of SSOR-iteration correction ',
     > 'for fifth pde  = ',1pe12.5)
 1007 format (1x/1x,'RMS-norm of steady-state residual for ',
     > 'first pde  = ',1pe12.5/,
     > 1x,'RMS-norm of steady-state residual for ',
     > 'second pde = ',1pe12.5/,
     > 1x,'RMS-norm of steady-state residual for ',
     > 'third pde  = ',1pe12.5/,
     > 1x,'RMS-norm of steady-state residual for ',
     > 'fourth pde = ',1pe12.5/,
     > 1x,'RMS-norm of steady-state residual for ',
     > 'fifth pde  = ',1pe12.5)
 
      end
