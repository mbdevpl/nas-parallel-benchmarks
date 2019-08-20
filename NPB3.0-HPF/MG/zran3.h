      include 'grid.h'      
      integer, intent(in) :: n1, n2, n3, k, nx, ny
      double precision, intent(inout) :: z(:,:,:)
c!hpf$    align(*,*,:) with grid :: z
