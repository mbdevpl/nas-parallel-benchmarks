      include 'grid.h'
      integer, intent (in) :: n1, n2, n3
      double precision, intent (in) :: v(:,:,:,:)
      double precision, intent (inout) :: u(:,:,:,:),r(:,:,:,:)
      double precision, intent (in) :: a(0:3),c(0:3)
c!hpf$    align(*,*,:,*) with grid :: u,v,r
!hpf$    align(*,:,:,*) with grid :: u,v,r
