c---------------------------------------------------------------------   

      pure extrinsic (hpf_local) subroutine 
     >              Swarztrauber(is,m,len,n,x,xd1,exponent)
      implicit none
      include 'global.h'

c---------------------------------------------------------------------
c   Computes NY N-point complex-to-complex FFTs of X using an algorithm due
c   to Swarztrauber.  X is both the input and the output array, while Y is a 
c   scratch array.  It is assumed that N = 2^M.  Before calling 
c   Swarztrauber to 
c   perform FFTs
c---------------------------------------------------------------------
  
      integer, intent(in) :: is,m,len,n,xd1
      double complex, dimension(:,:), intent(inout) ::  x
      double complex, dimension(:), intent(in) :: exponent
      integer i,j,l,mx      
      double complex u1,x11,x21
      integer k, n1,li,lj,lk,ku,i11,i12,i21,i22
      integer BlockStart,BlockEnd
      
      call timer_start(4)
c---------------------------------------------------------------------
c   Perform one variant of the Stockham FFT.
c---------------------------------------------------------------------


cc     This is an extra loop for fitting a slab into cash
cc     Origin 2000 cash size 32K = 4096 of 8 byte words
      fftblock = 4*4096/n
      
      do BlockStart = 1, len, fftblock
        BlockEnd = BlockStart + fftblock - 1
        if ( BlockEnd .gt. len) then
	  BlockEnd = len
	endif
      do l = 1, m, 2
        n1 = n / 2
        lk = 2 ** (l - 1)
        li = 2 ** (m - l)
        lj = 2 * lk
        ku = li + 1

        do i = 0, li - 1
          i11 = i * lk + 1
          i12 = i11 + n1
          i21 = i * lj + 1
          i22 = i21 + lk
        
          if (is .ge. 1) then
            u1 = exponent(ku+i)
          else
            u1 = dconjg (exponent(ku+i))
          endif
          do k = 0, lk - 1
            do j = BlockStart, BlockEnd
              x11 = x(j,i11+k)
              x21 = x(j,i12+k)
              scr(j,i21+k) = x11 + x21
              scr(j,i22+k) = u1 * (x11 - x21)
            end do
          end do
        end do

        if (l .eq. m) then
          x(BlockStart:BlockEnd,1:n) = scr(BlockStart:BlockEnd,1:n)
        else
          n1 = n / 2
          lk = 2 ** l
          li = 2 ** (m - l - 1)
          lj = 2 * lk
          ku = li + 1

          do i = 0, li - 1
            i11 = i * lk + 1
            i12 = i11 + n1
            i21 = i * lj + 1
            i22 = i21 + lk
        
            if (is .ge. 1) then
              u1 = exponent(ku+i)
            else
              u1 = dconjg (exponent(ku+i))
            endif
            do k = 0, lk - 1
              do j = BlockStart, BlockEnd
                x11 = scr(j,i11+k)
                x21 = scr(j,i12+k)
                x(j,i21+k) = x11 + x21
                x(j,i22+k) = u1 * (x11 - x21)
              end do
            end do
          end do
        endif
      end do
      end do

      call timer_stop(4)

      return
      end

c---------------------------------------------------------------------

      extrinsic(HPF)
     > subroutine fftXYZ(sign,x,x3,exp1,exp2,exp3,n1,n2,n3)
       include 'fftXYZ.h'
       include 'global.h'
       
         interface
           pure extrinsic (hpf_local) subroutine 
     >            Swarztrauber(is, m,len,n, x,xd1,exponent)
               integer, intent(in) :: is,m,len,n,xd1
               double complex, dimension(:,:), intent(inout) ::  x
               double complex, dimension(:), intent(in) :: exponent
           end subroutine Swarztrauber
           
         end interface
       double complex x2(n2+1,n1,n3)
!hpf$    template grid2(n3)
!hpf$    distribute(block) :: grid2
!hpf$    align(*,*,:) with grid2 :: x2
       integer i, j, k, log
         call timer_start(3)        
	 log = ilog2( n2 )
         call timer_start(7)        
!HPF$ INDEPENDENT         
         do k = 1, n3
            call Swarztrauber(sign,log,n1,n2,x(:,:,k),n1,exp2) 
         end do
         call timer_stop(7)
	 log = ilog2( n1)
         call timer_start(8) 
!HPF$ INDEPENDENT
         do k = 1, n3
         call timer_start(6) 
           do j = 1, n2
           do i = 1, n1
             plane(j,i) = x(i,j,k)
           end do
           end do
         call timer_stop(6)
           call Swarztrauber(sign,log,n2,n1,plane,maxdim,exp1)     
         call timer_start(6) 
           do j = 1, n2
           do i = 1, n1
             x(i,j,k)=plane(j,i)
           end do
           end do
         call timer_stop(6)
         end do
         call timer_stop(8)

         call timer_start(5)
         do k = 1, n3
         do i = 1, n2
         do j = 1, n1
           x3(j,k,i) = x(j,i,k)
         end do
         end do
         end do
         call timer_stop(5)
         
       	 log = ilog2(n3)
         call timer_start(9)

!HPF$ independent
         do k = 1, n2
           call Swarztrauber(sign,log,n1,n3,x3(:,:,k),n1,exp3)
         end do                 
         call timer_stop(9)

         call timer_stop(3)

         return
      end
