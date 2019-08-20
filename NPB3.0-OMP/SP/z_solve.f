
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine z_solve

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c this function performs the solution of the approximate factorization
c step in the z-direction for all five matrix components
c simultaneously. The Thomas algorithm is employed to solve the
c systems for the z-lines. Boundary conditions are non-periodic
c---------------------------------------------------------------------

       include 'header.h'

       integer j


c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c Prepare for z-solve, array redistribution
c---------------------------------------------------------------------

       if (timeron) call timer_start(t_zsolve)
!$omp parallel do default(shared) private(j)
       do   j = 1, ny2
          call z_solve_j(j)
       end do
       if (timeron) call timer_stop(t_zsolve)

       if (timeron) call timer_start(t_tzetar)
       call tzetar
       if (timeron) call timer_stop(t_tzetar)

       return
       end







