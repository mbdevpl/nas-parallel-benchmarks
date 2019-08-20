
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine x_solve

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c this function performs the solution of the approximate factorization
c step in the x-direction for all five matrix components
c simultaneously. The Thomas algorithm is employed to solve the
c systems for the x-lines. Boundary conditions are non-periodic
c---------------------------------------------------------------------

       include 'header.h'

       integer k

c---------------------------------------------------------------------
c---------------------------------------------------------------------

       if (timeron) call timer_start(t_xsolve)
!$omp parallel do default(shared) private(k)
       do  k = 1, nz2
          call x_solve_k(k)
       end do
       if (timeron) call timer_stop(t_xsolve)

c---------------------------------------------------------------------
c      Do the block-diagonal inversion          
c---------------------------------------------------------------------
       if (timeron) call timer_start(t_ninvr)
       call ninvr
       if (timeron) call timer_stop(t_ninvr)

       return
       end
    






