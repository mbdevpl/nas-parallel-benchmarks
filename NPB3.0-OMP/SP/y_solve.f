
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine y_solve

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c this function performs the solution of the approximate factorization
c step in the y-direction for all five matrix components
c simultaneously. The Thomas algorithm is employed to solve the
c systems for the y-lines. Boundary conditions are non-periodic
c---------------------------------------------------------------------

       include 'header.h'

       integer k


c---------------------------------------------------------------------
c---------------------------------------------------------------------

       if (timeron) call timer_start(t_ysolve)
!$omp parallel do default(shared) private(k)
       do  k = 1, nz2
          call y_solve_k(k)
       end do
       if (timeron) call timer_stop(t_ysolve)


       if (timeron) call timer_start(t_pinvr)
       call pinvr
       if (timeron) call timer_stop(t_pinvr)

       return
       end







