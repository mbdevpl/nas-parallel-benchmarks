!
! FT verification routine.
!
      subroutine verify(n1,n2,n3,nt,cksum,verified)
        implicit none
        include 'npbparams.h'
       logical verified
!
! Arguments.
!
         integer, intent(in) :: n1, n2, n3, nt
         double complex, dimension(niter_default), intent(in) :: cksum
!
! Local variables.
!
         logical success
         logical temp(niter_default)
         integer kt,i
         double complex cexpd(20)
         real*8 epsilon
!
! Initialize tolerance level and success flag.
!
         epsilon = 1.0e-12
         success = .true.
!
         if ((n1 .eq. 64) .and. (n2 .eq. 64) .and.                 
     &            (n3 .eq. 64) .and. (nt .eq. 6)) then
!
! Class S reference values.
!
            cexpd(1) = dcmplx(554.6087004964D0, 484.5363331978D0)
            cexpd(2) = dcmplx(554.6385409189D0, 486.5304269511D0)
            cexpd(3) = dcmplx(554.6148406171D0, 488.3910722336D0)
            cexpd(4) = dcmplx(554.5423607415D0, 490.1273169046D0)
            cexpd(5) = dcmplx(554.4255039624D0, 491.7475857993D0)
            cexpd(6) = dcmplx(554.2683411902D0, 493.2597244941D0)
            
         else if ((n1 .eq. 128) .and. (n2 .eq. 128) .and.                 
     &            (n3 .eq. 32) .and. (nt .eq. 6)) then
!
! Class W reference values.
!
            cexpd(1) = dcmplx(567.3612178944D0, 529.3246849175D0)
            cexpd(2) = dcmplx(563.1436885271D0, 528.2149986629D0)
            cexpd(3) = dcmplx(559.4024089970D0, 527.0996558037D0)
            cexpd(4) = dcmplx(556.0698047020D0, 526.0027904925D0)
            cexpd(5) = dcmplx(553.0898991250D0, 524.9400845633D0)
            cexpd(6) = dcmplx(550.4159734538D0, 523.9212247086D0)
!
         else if ((n1 .eq. 256) .and. (n2 .eq. 256) .and.               
     &            (n3 .eq. 128) .and. (nt .eq. 6)) then
!
! Class A reference values.
!
            cexpd(1) = dcmplx(504.6735008193D0, 511.4047905510D0)
            cexpd(2) = dcmplx(505.9412319734D0, 509.8809666433D0)
            cexpd(3) = dcmplx(506.9376896287D0, 509.8144042213D0)
            cexpd(4) = dcmplx(507.7892868474D0, 510.1336130759D0)
            cexpd(5) = dcmplx(508.5233095391D0, 510.4914655194D0)
            cexpd(6) = dcmplx(509.1487099959D0, 510.7917842803D0)
!
         else if ((n1 .eq. 512) .and. (n2 .eq. 256) .and.               
     &            (n3 .eq. 256) .and. (nt .eq. 20)) then
!
! Class B reference values.
!
             cexpd(1) = dcmplx(5.177643571579D+02,5.077803458597D+02)
             cexpd(2) = dcmplx(5.154521291263D+02,5.088249431599D+02)
             cexpd(3) = dcmplx(5.146409228649D+02,5.096208912659D+02)
             cexpd(4) = dcmplx(5.142378756213D+02,5.101023387619D+02)
             cexpd(5) = dcmplx(5.139626667737D+02,5.103976610617D+02)
             cexpd(6) = dcmplx(5.137423460082D+02,5.105948019802D+02)
             cexpd(7) = dcmplx(5.135547056878D+02,5.107404165783D+02)
             cexpd(8) = dcmplx(5.133910925466D+02,5.108576573661D+02)
             cexpd(9) = dcmplx(5.132470705390D+02,5.109577278523D+02)
            cexpd(10) = dcmplx(5.131197729984D+02,5.110460304483D+02)
            cexpd(11) = dcmplx(5.130070319283D+02,5.111252433800D+02)
            cexpd(12) = dcmplx(5.129070537032D+02,5.111968077718D+02)
            cexpd(13) = dcmplx(5.128182883502D+02,5.112616233064D+02)
            cexpd(14) = dcmplx(5.127393733383D+02,5.113203605551D+02)
            cexpd(15) = dcmplx(5.126691062020D+02,5.113735928093D+02)
            cexpd(16) = dcmplx(5.126064276004D+02,5.114218460548D+02)
            cexpd(17) = dcmplx(5.125504076570D+02,5.114656139760D+02)
            cexpd(18) = dcmplx(5.125002331720D+02,5.115053595966D+02)
            cexpd(19) = dcmplx(5.124551951846D+02,5.115415130407D+02)
            cexpd(20) = dcmplx(5.124146770029D+02,5.115744692211D+02)
         else if ((n1 .eq. 512) .and. (n2 .eq. 512) .and.               
     &            (n3 .eq. 512) .and. (nt .eq. 20)) then
!
! Class C reference values.
!
             cexpd(1) = dcmplx(5.195078707457D+02, 5.149019699238D+02)
             cexpd(2) = dcmplx(5.155422171134D+02, 5.127578201997D+02)
             cexpd(3) = dcmplx(5.144678022222D+02, 5.122251847514D+02)
             cexpd(4) = dcmplx(5.140150594328D+02, 5.121090289018D+02)
             cexpd(5) = dcmplx(5.137550426810D+02, 5.121143685824D+02)
             cexpd(6) = dcmplx(5.135811056728D+02, 5.121496764568D+02)
             cexpd(7) = dcmplx(5.134569343165D+02, 5.121870921893D+02)
             cexpd(8) = dcmplx(5.133651975661D+02, 5.122193250322D+02)
             cexpd(9) = dcmplx(5.132955192805D+02, 5.122454735794D+02)
            cexpd(10) = dcmplx(5.132410471738D+02, 5.122663649603D+02)
            cexpd(11) = dcmplx(5.131971141679D+02, 5.122830879827D+02)
            cexpd(12) = dcmplx(5.131605205716D+02, 5.122965869718D+02)
            cexpd(13) = dcmplx(5.131290734194D+02, 5.123075927445D+02)
            cexpd(14) = dcmplx(5.131012720314D+02, 5.123166486553D+02)
            cexpd(15) = dcmplx(5.130760908195D+02, 5.123241541685D+02)
            cexpd(16) = dcmplx(5.130528295923D+02, 5.123304037599D+02)
            cexpd(17) = dcmplx(5.130310107773D+02, 5.123356167976D+02)
            cexpd(18) = dcmplx(5.130103090133D+02, 5.123399592211D+02)
            cexpd(19) = dcmplx(5.129905029333D+02, 5.123435588985D+02)
            cexpd(20) = dcmplx(5.129714421109D+02, 5.123465164008D+02)
!
         else
!
            write (*,    120) 'NOT DONE'
            success = .false.
!
         end if
!
! Verification test for results.
!
	 verified=.false.
         if (success) then

            success = all (abs ((cksum(1:nt) - cexpd(1:nt)) /
     &                     cexpd(1:nt)) .le. epsilon)

            if (success) then
	       verified=.true.
               write (*,    120) 'PASSED'
            else
               write (*,    120) 'FAILED'
            end if

  120       format (' Verification test for FT ', a)
         end if
!
         return
      end
