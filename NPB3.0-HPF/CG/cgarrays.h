      double precision           z(na+2),
     >                        	 x(na+2),
     >                        	 a(nz),
     >                        	 p(na+2),
     >                        	 pd(na+2),
     >                        	 q(na+2),
     >                        	 r(na+2)
      double precision adistr(rowmaxnz/2*2+1,na+2)
      integer colidxd(rowmaxnz/2*2+1,na+2)
!hpf$ distribute (block) :: z
!hpf$ align (:) with z :: q,r,x,pd
!hpf$ align p(*) with z(*)
!hpf$ align (*,:) with z :: adistr, colidxd
