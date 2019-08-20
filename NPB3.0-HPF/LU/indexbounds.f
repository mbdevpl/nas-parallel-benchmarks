      subroutine InitSectionsN
      implicit none
      include 'applu.incl'
      include 'indexbounds.h'

        xsize = isiz1
        ysize = isiz2
        zsize = isiz3

        section1 = xsize + 2
        section2 = xsize + ysize - 1
        section3 = xsize + ysize + zsize - 3

      return
      end

      subroutine InitLoopBoundsLocalN( k, LC)
      implicit none
      include 'indexbounds.h'
      include 'loopcontrol.h'
      integer, intent(in) :: k
      type(LoopControlStruct), intent(inout) :: LC

        LC%section1 = section1
        LC%section2 = section2
        LC%section3 = section3

        LC%xsize = xsize
        LC%ysize = ysize
        LC%zsize = zsize

        LC%updown = 0
        LC%updowne = 0
        LC%topseg = 0
        LC%botseg = 0

            if ( k .le. section1 ) then
              jfirst = 2
              jlast = k-4
            else if ( k .le. section2 ) then
              jfirst = 2
              jlast = ysize-1
            else
              jfirst = k-section2+1
              jlast = xsize-1
            endif

      return
      end


         pure extrinsic (hpf_local)
     >      subroutine GetIBoundsLocalN( j, k, LC)

      implicit none
      include 'indexbounds.h'
      include 'loopcontrol.h'

      integer, intent(in) :: j,k
      type(LoopControlStruct), intent(inout) :: LC

             if ( k .le. LC%section1) then
               LC%first = 2
               LC%last = jlast + 2 - j
             else if ( k .le. LC%section2) then
               if ( j .lt. k - LC%xsize ) then
                 LC%first = k - j - LC%xsize + 1
                 LC%last = LC%zsize - 1
               else
                 LC%first = 2
                 LC%last =  k - j - 2
               endif
             else
               LC%first = k - j - LC%xsize + 1
               LC%last = LC%xsize - 1
             endif

      return
      end

      extrinsic (hpf) subroutine SetLoopBounds
      implicit none
      include 'applu.incl'
      include 'indexbounds.h'
      include 'loopcontrol.h'
      integer sumdim, k, j
      type(LoopControlStruct) :: loopControl

      sumdim = isiz1+isiz2+isiz3-3
      do k = 6, sumdim
      call InitLoopBoundsLocalN( k, loopControl)
         jlow(k) = jfirst
         jhigh(k) = jlast
         do j= jfirst, jlast
           call GetIBoundsLocalN(j,k,loopControl)
           ilow(k,j) = loopControl%first
           ihigh(k,j) = loopControl%last
         end do
      end do

      return
      end




