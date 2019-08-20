      include 'grid.h'      
      integer, intent (in) :: m1k, m2k, m3k, m1j, m2j, m3j,k
      double precision, intent (in) :: r(:,:,:)
      double precision, intent (out) :: s(:,:,:)
c!hpf$    align(*,*,:) with grid :: r,s
!hpf$    align(*,:,:) with grid :: r,s
