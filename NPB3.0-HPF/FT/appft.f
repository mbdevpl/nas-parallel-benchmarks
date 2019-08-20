      extrinsic (HPF) subroutine appft(niter,tmax,verified,tflag)
         include  'global.h'

         integer, intent(in) :: niter
         double precision, intent(out) :: tmax
         logical, intent(out) :: verified
         logical, intent(in) :: tflag
         
         interface 
               extrinsic(HPF)
     >           subroutine compute_initial_conditions(x,d1,d2,d3)
                 integer, intent(in) :: d1,d2,d3
                 double complex, dimension(:,:,:), intent(inout) :: x
!hpf$    template grid(d3+1)
!hpf$    distribute(block) :: grid
!hpf$    align(*,*,:) with grid :: x
               end subroutine compute_initial_conditions
               
               extrinsic(HPF)
     >           subroutine fftXYZ(sign,x,x3,exp1,exp2,exp3,n1,n2,n3)
                 include 'fftXYZ.h'
               end subroutine fftXYZ

               extrinsic (HPF)
     >           subroutine CalculateChecksum(Sum,iterN,u,d1,d2,d3)
                 integer, intent (in) :: iterN
                 double complex, dimension(:,:,:), intent(in) :: u
!hpf$    template grid(nz+1)
!hpf$    distribute(block) :: grid
!hpf$    align(*,*,:) with grid :: u
                 double complex, intent (out) :: Sum
                 integer, intent (in) :: d1, d2, d3
               end subroutine CalculateChecksum
         end interface
!
! Local variables
!
         integer i, j, k, kt, n12, n22, n32
         real*8 ap, scale
         double complex checksum(niter_default)
         
         double complex xtr(ny+1,nx,nz)
!hpf$    template gridx(nz)
!hpf$    distribute(block) :: gridx
!hpf$    align (*,*,:) with gridx(:) :: xtr

         double complex xnt(ny+1,nz,nx), y(ny+1,nz,nx)
         real*8 twiddle(ny+1,nz,nx)
!hpf$    template gridy(ny+1,nz,nx)
!hpf$    distribute(*,*,block) :: gridy
!hpf$    align (:,:,:) with gridy(:,:,:) :: xnt, y, twiddle

         double complex exp1(nx), exp2(ny), exp3(nz)
!hpf$    align exp1(*) with gridy(*,*,*)
!hpf$    align exp2(*) with gridy(*,*,*)
!hpf$    align exp3(*) with gridy(*,*,*)
         
         do i=1,15
           call timer_clear(i)
	 end do         
!
! Initialize AP, and X1
!
         scale = 1.d0 / dble(nx * ny * nz)
         ap = - 4.d0 * alpha * pi ** 2
         n12 = nx/2
         n22 = ny/2
         n32 = nz/2

         call timer_start(2)         
         call compute_initial_conditions(xtr,ny,nx,nz)
         
         call CompExp( nx, exp1 )
         call CompExp( ny, exp2 )
         call CompExp( nz, exp3 )           
         call fftXYZ(1, xtr, y, exp1, exp2, exp3,ny,nx,nz)
         call timer_stop(2)      

         do i=1,15
           if (i.ne.2) call timer_clear(i)
	 end do         

         call timer_start(1)
         if (tflag) call timer_start(13)
!HPF$    independent
         do i = 1, nx
	   ii = i-1-((i-1)/n12)*nx
	   ii2 = ii*ii
           do k = 1, nz
	     kk = k-1-((k-1)/n32)*nz
	     ik2 = ii2 + kk*kk
             do j = 1, ny
	         jj = j-1-((j-1)/n22)*ny
                 twiddle(j,k,i) = 
     *                exp( ap*dble(jj*jj + ik2) )
               end do
            end do
         end do
         if (tflag) call timer_stop(13)      
         
         if (tflag) call timer_start(12)
         call compute_initial_conditions(xtr,ny,nx,nz)             
         if (tflag) call timer_stop(12)      
         if (tflag) call timer_start(15)      
         call fftXYZ(1, xtr, y, exp2, exp1, exp3,ny,nx,nz)
         if (tflag) call timer_stop(15)      

         do kt = 1, niter
	   if (tflag) call timer_start(11)      
           do i = 1, nx
             do k = 1, nz
               do j = 1, ny
                   xnt(j,k,i)=y(j,k,i)*twiddle(j,k,i)**kt
                 end do
              end do
           end do
           if (tflag) call timer_stop(11)      
           if (tflag) call timer_start(15)      
           call fftXYZ(-1,xnt,xtr,exp2,exp3,exp1,ny,nz,nx)
           if (tflag) call timer_stop(15)      
           if (tflag) call timer_start(10)      
           call CalculateChecksum(checksum(kt),kt,xtr,ny+1,nx,nz)           
           if (tflag) call timer_stop(10)      
         end do
!
! Verification test.
!
         if (tflag) call timer_start(14)      
         call verify(nx,ny,nz,niter,checksum,verified)
         if (tflag) call timer_stop(14)      
         call timer_stop(1)
                                  
         tmax = timer_read(1)
	 if (.not.tflag) return

         print *, 'FT time =                     ', tmax
         print *, 'WarmUp time =                 ', timer_read(2)
         print *, 'ffXYZ body time =             ', timer_read(3)
         print *, 'Swarztrauber body time =      ', timer_read(4)
         print *, 'Redistribution time =         ', timer_read(5)
         print *, 'Transposition time =          ', timer_read(6)
         print *, 'X time =                      ', timer_read(7)
         print *, 'Y time =                      ', timer_read(8)
         print *, 'Z time =                      ', timer_read(9)
         print *, 'CalculateChecksum =           ', timer_read(10)
         print *, 'evolve =                      ', timer_read(11)
         print *, 'compute_initial_conditions =  ', timer_read(12)
         print *, 'twiddle =                     ', timer_read(13)
         print *, 'verify =                      ', timer_read(14)
         print *, 'fftXYZ =                      ', timer_read(15)

         print *, 'Benchmark time = ', tmax, ' seconds'

         return
      end
