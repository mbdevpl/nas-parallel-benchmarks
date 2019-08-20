         integer, intent(in) :: sign,n1,n2,n3
!
         double complex, dimension(:,:,:), intent(inout) :: x
!hpf$    template grid1(size(x,3))
!hpf$    distribute(block) :: grid1
!hpf$    align(*,*,:) with *grid1 :: x

         double complex, dimension(:,:,:), intent(out) :: x3
!hpf$    template grid3(size(x3,3))
!hpf$    distribute(block) :: grid3
!hpf$    align(*,*,:) with *grid3 :: x3
         double complex, dimension(:), intent(in) :: exp1, exp2, exp3
