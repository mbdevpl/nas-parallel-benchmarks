c#ifdef classS
cc!hpf$    template grid(36)
c!hpf$    template grid(36,36)
c#endif
c#ifdef classW
cc!hpf$    template grid(68)
c!hpf$    template grid(68,68)
c#endif
c#ifdef classA
c!hpf$    template grid(262)
!hpf$    template grid(262,262)
c#endif
c#ifdef classB
c!hpf$    template grid(262)
c!hpf$    template grid(262,262)
c#endif
c#ifdef classC
c!hpf$    template grid(516)
c!hpf$    template grid(516,516)
c#endif
cc!hpf$    distribute(block) :: grid
!hpf$    distribute(*,block) :: grid

