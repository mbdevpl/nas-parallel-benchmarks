      integer, intent(in) ::  ldmx, ldmy, ldmz, l
      double precision, intent(in) :: omega
      double precision, intent(inout) :: 
     >       v(5, ldmx/2*2+1, ldmy/2*2+1, ldmz)

           double precision, intent(in) :: 
     >        d( 5, 5, ldmx/2*2+1, ldmy),
     >        udx( 5, 5, ldmx/2*2+1, ldmy),
     >        udy( 5, 5, ldmx/2*2+1, ldmy),
     >        udz( 5, 5, ldmx/2*2+1, ldmy )

      include 'npbparams.h'
      include 'projarea.h'       
 
!hpf$    align(*,:,:,*) with ProjArea :: v
!hpf$    align(*,*,:,:) with ProjArea :: udz,udy,udx,d
