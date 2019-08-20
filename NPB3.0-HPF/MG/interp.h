      include 'grid.h'
      integer, intent (in) :: mm1, mm2, mm3, n1, n2, n3,k
      double precision, intent (in) :: z(:,:,:)
      double precision, intent (inout) :: u(:,:,:)     
cc!hpf$    align(*,*,:) with grid :: z, u
!hpf$    align(*,:,:) with grid :: z, u
