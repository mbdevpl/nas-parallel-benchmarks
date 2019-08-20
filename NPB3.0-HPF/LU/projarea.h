!hpf$    template ProjArea(isiz1/2*2+1,isiz2/2*2+1)
!hpf$    distribute(*,block) ::ProjArea
c!hpf$    distribute(block,block) ::ProjArea
         integer jlow(isiz1+isiz2+isiz3+3), jhigh(isiz1+isiz2+isiz3+3)
         integer ilow(isiz1+isiz2+isiz3+3,isiz2/2*2+1),
     >	         ihigh(isiz1+isiz2+isiz3+3,isiz2/2*2+1)
         common/indbounds/ jlow, jhigh, ilow, ihigh
!hpf$    distribute(*,*) :: ilow, ihigh
!hpf$    distribute(*) :: jlow, jhigh
