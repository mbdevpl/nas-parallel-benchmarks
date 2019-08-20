      include 'grid.h'
      integer, intent (in) :: n1,n2,n3,k
      double precision, intent (in) ::  u(:,:,:),v(:,:,:),a(0:3)
      double precision, intent (out) :: r(:,:,:)
c!hpf$ align(*,*,:) with grid :: u,v,r
!hpf$ align(*,:,:) with grid :: u,v,r
