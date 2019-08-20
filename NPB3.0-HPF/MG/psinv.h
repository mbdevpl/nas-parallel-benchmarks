      include 'grid.h'
      integer, intent (in) :: n1,n2,n3,k
      double precision, intent (in) :: r(:,:,:),c(0:3)
      double precision, intent (inout) ::  u(:,:,:)
cc!hpf$ align(*,*,:) with grid :: r,u         
!hpf$ align(*,:,:) with grid :: r,u         
