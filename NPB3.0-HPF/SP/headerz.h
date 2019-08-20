c---------------------------------------------------------------------
c   Temp. arrays for redistribution
c---------------------------------------------------------------------
      double precision 
     >   wsz     (0:IMAX/2*2, 0:JMAX-1, 0:KMAX-1    ),
     >   rho_iz  (0:IMAX/2*2, 0:JMAX-1, 0:KMAX-1    ),
     >   speedz  (0:IMAX/2*2, 0:JMAX-1, 0:KMAX-1    ),
     >   rhsz    (5,  0:IMAX/2*2, 0:JMAX-1, 0:KMAX-1)

      common /work_fields/  wsz, rho_iz, speedz, rhsz

!hpf$    template gridy(0:JMAX-1)
!hpf$    distribute(block) :: gridy
!hpf$    align(*,*,:,*) with gridy :: rhsz
!hpf$    align(*,:,*) with gridy :: wsz, rho_iz, speedz
