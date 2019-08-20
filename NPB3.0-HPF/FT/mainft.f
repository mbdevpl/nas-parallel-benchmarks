!-------------------------------------------------------------------------!
!                                                                         !
!        N  A  S     P A R A L L E L     B E N C H M A R K S  3.0         !
!                                                                         !
!                         H P F     V E R S I O N                         !
!                                                                         !
!                                   F T                                   !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    This benchmark is an HPF version of the NPB FT code.                 !
!                                                                         !
!    Permission to use, copy, distribute and modify this software         !
!    for any purpose with or without fee is hereby granted.  We           !
!    request, however, that all derived work reference the NAS            !
!    Parallel Benchmarks 3.0. This software is provided "as is"           !
!    without express or implied warranty.                                 !
!                                                                         !
!    Information on NPB 3.0, including the technical report, the          !
!    original specifications, source code, results and information        !
!    on how to submit new results, is available at:                       !
!                                                                         !
!           http://www.nas.nasa.gov/Software/NPB/                         !
!                                                                         !
!    Send comments or suggestions to  npb@nas.nasa.gov                    !
!                                                                         !
!          NAS Parallel Benchmarks Group                                  !
!          NASA Ames Research Center                                      !
!          Mail Stop: T27A-1                                              !
!          Moffett Field, CA   94035-1000                                 !
!                                                                         !
!          E-mail:  npb@nas.nasa.gov                                      !
!          Fax:     (650) 604-3957                                        !
!                                                                         !
!-------------------------------------------------------------------------!

c---------------------------------------------------------------------
c
c Authors: D. Bailey
c          W. Saphir
c HPF version: 
c          M.Frumkin
c
c---------------------------------------------------------------------

c---------------------------------------------------------------------

c---------------------------------------------------------------------
c FT benchmark
c---------------------------------------------------------------------
      program mainft
         include 'npbparams.h'

	 integer nprocs, niter, fstatus
	 character*8  class
         double precision tmax, mflops
	 logical timers_enabled
	 
         open (unit=2,file='timer.flag',status='old',iostat=fstatus)
         if (fstatus .eq. 0) then
            timers_enabled = .true.
            close(2)
         else
            timers_enabled = .false.
         endif

	 niter=niter_default
         nprocs=number_of_processors()

      	 write(*, 1000)
      	 write(*, 1001) nx, ny, nz
      	 write(*, 1002) niter
      	 write(*, 1003) nprocs

 1000    format(//,' NAS Parallel Benchmarks (NPB3.0-HPF)',
     >          ' - FT Benchmark', /)
 1001    format(' Size                : ', i3, 'x', i3, 'x', i3)
 1002    format(' Iterations          :     ', i7)
 1003    format(' Number of processes :     ', i7/)

         class(1:1)='U'
	 call getclass(class)
!
         call appft (niter, tmax, verified, timers_enabled)
!
        if( tmax .ne. 0. ) then
           mflops = 1.0d-6*float(ntotal) *
     >             (14.8157+7.19641*log(float(ntotal))
     >          +  (5.23518+7.21113*log(float(ntotal)))*niter)
     >                 /tmax
        else
           mflops = 0.0
        endif
        call print_results('FT', class(1:1), nx, 
     >  ny, nz, niter, nprocs,
     >  tmax, mflops, '          floating point', 
     >  verified, npbversion,compiletime, cs1, cs2, cs3, cs4, cs5, 
     >  cs6, '(none)')
!
      end
      
      subroutine getclass( class)
        include 'npbparams.h'
	character*8  class
        if ((nx .eq. 64) .and. (ny .eq. 64) .and.                 
     &      (nz .eq. 64) .and. (niter_default .eq. 6)) then
          class(1:1)='S'
        else if ((nx .eq. 128) .and. (ny .eq. 128) .and.                 
     &           (nz .eq. 32) .and. (niter_default .eq. 6)) then
          class(1:1)='W'
        else if ((nx .eq. 256) .and. (ny .eq. 256) .and.                 
     &           (nz .eq. 128) .and. (niter_default .eq. 6)) then
          class(1:1)='A'
        else if ((nx .eq. 512) .and. (ny .eq. 256) .and.                 
     &           (nz .eq. 256) .and. (niter_default .eq. 20)) then
          class(1:1)='B'
        else if ((nx .eq. 512) .and. (ny .eq. 512) .and.                 
     &           (nz .eq. 512) .and. (niter_default .eq. 20)) then
          class(1:1)='C'
	endif
	
	return
      end
      
