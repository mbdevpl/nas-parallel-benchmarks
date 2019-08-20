
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine jacld( l)

c---------------------------------------------------------------------
c---------------------------------------------------------------------


c---------------------------------------------------------------------
c   compute the lower triangular part of the jacobian matrix
c---------------------------------------------------------------------

      implicit none

      include 'applu.incl'
      
c---------------------------------------------------------------------
c  input parameters
c---------------------------------------------------------------------
      integer l

c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------
      integer i, j, idx, k 
      double precision  r43
      double precision  c1345
      double precision  c34
      double precision  tmp1, tmp2, tmp3
      double precision  u0(5), ux(5), uy(5), uz(5)
      double precision  u0l(5,isiz1/2*2+1,isiz2),
     >                  uxl(5,isiz1/2*2+1,isiz2),
     >                  uyl(5,isiz1/2*2+1,isiz2),
     >                  uzl(5,isiz1/2*2+1,isiz2),
     >                  uylt(5,isiz1/2*2+1,isiz2)
      common/workarr/ u0l, uxl, uyl, uzl, uylt
!hpf$    align(*,:,:) with ProjArea :: u0l, uxl, uyl, uzl, uylt  

      r43 = ( 4.0d+00 / 3.0d+00 )
      c1345 = c1 * c3 * c4 * c5
      c34 = c3 * c4
 
!HPF$ independent, new(k)
         do j= jlow(l), jhigh(l) 
           do i=ilow(l,j), ihigh(l,j)
             k = l-i-j	 
             do idx = 1, 5
               u0l(idx,i,j) = u(idx,i,j,k)
               uxl(idx,i,j) = u(idx,i-1,j,k)
               uzl(idx,i,j) = u(idx,i,j,k-1)               
             end do           
           end do
	 end do

!HPF$ independent, new(k)
         do j= jlow(l), jhigh(l) 
           do i=ilow(l,j), ihigh(l,j)
             k = l-i-j	 
             do idx = 1, 5
               uyl(idx,i,j-1) = u(idx,i,j-1,k)
             end do           
           end do
	 end do
        uylt(:,jlow(l):jhigh(l),jlow(l):jhigh(l)) =
     >                uyl(:,jlow(l):jhigh(l),jlow(l)-1:jhigh(l)-1)
     
!HPF$ independent, new( u0, ux, uy, uz,tmp1,tmp2,tmp3,k)
         do j= jlow(l), jhigh(l) 
           do i=ilow(l,j), ihigh(l,j)
	 
               do idx = 1, 5
                 u0(idx) = u0l(idx,i,j)
                 ux(idx) = uxl(idx,i,j)
                 uy(idx) = uylt(idx,i,j)
                 uz(idx) = uzl(idx,i,j)               
               end do           
c---------------------------------------------------------------------
c   form the block daigonal
c---------------------------------------------------------------------
c               tmp1 = rho_i(i,j,k)
               tmp1 = 1.0d+00 / u0(1)
               tmp2 = tmp1 * tmp1
               tmp3 = tmp1 * tmp2

               d(1,1,i,j) =  1.0d+00
     >                       + dt * 2.0d+00 * (   tx1 * dx1
     >                                          + ty1 * dy1
     >                                          + tz1 * dz1 )
               d(1,2,i,j) =  0.0d+00
               d(1,3,i,j) =  0.0d+00
               d(1,4,i,j) =  0.0d+00
               d(1,5,i,j) =  0.0d+00

               d(2,1,i,j) = -dt * 2.0d+00
     >          * (  tx1 * r43 + ty1 + tz1  )
     >          * c34 * tmp2 * u0(2)
               d(2,2,i,j) =  1.0d+00
     >          + dt * 2.0d+00 * c34 * tmp1 
     >          * (  tx1 * r43 + ty1 + tz1 )
     >          + dt * 2.0d+00 * (   tx1 * dx2
     >                             + ty1 * dy2
     >                             + tz1 * dz2  )
               d(2,3,i,j) = 0.0d+00
               d(2,4,i,j) = 0.0d+00
               d(2,5,i,j) = 0.0d+00

               d(3,1,i,j) = -dt * 2.0d+00
     >           * (  tx1 + ty1 * r43 + tz1  )
     >           * c34 * tmp2 * u0(3)
               d(3,2,i,j) = 0.0d+00
               d(3,3,i,j) = 1.0d+00
     >         + dt * 2.0d+00 * c34 * tmp1
     >              * (  tx1 + ty1 * r43 + tz1 )
     >         + dt * 2.0d+00 * (  tx1 * dx3
     >                           + ty1 * dy3
     >                           + tz1 * dz3 )
               d(3,4,i,j) = 0.0d+00
               d(3,5,i,j) = 0.0d+00

               d(4,1,i,j) = -dt * 2.0d+00
     >           * (  tx1 + ty1 + tz1 * r43  )
     >           * c34 * tmp2 * u0(4)
               d(4,2,i,j) = 0.0d+00
               d(4,3,i,j) = 0.0d+00
               d(4,4,i,j) = 1.0d+00
     >         + dt * 2.0d+00 * c34 * tmp1
     >              * (  tx1 + ty1 + tz1 * r43 )
     >         + dt * 2.0d+00 * (  tx1 * dx4
     >                           + ty1 * dy4
     >                           + tz1 * dz4 )
               d(4,5,i,j) = 0.0d+00

               d(5,1,i,j) = -dt * 2.0d+00
     >  * ( ( ( tx1 * ( r43*c34 - c1345 )
     >     + ty1 * ( c34 - c1345 )
     >     + tz1 * ( c34 - c1345 ) ) * ( u0(2) ** 2 )
     >   + ( tx1 * ( c34 - c1345 )
     >     + ty1 * ( r43*c34 - c1345 )
     >     + tz1 * ( c34 - c1345 ) ) * ( u0(3) ** 2 )
     >   + ( tx1 * ( c34 - c1345 )
     >     + ty1 * ( c34 - c1345 )
     >     + tz1 * ( r43*c34 - c1345 ) ) * ( u0(4) ** 2 )
     >      ) * tmp3
     >   + ( tx1 + ty1 + tz1 ) * c1345 * tmp2 * u0(5) )

               d(5,2,i,j) = dt * 2.0d+00 * tmp2 * u0(2)
     > * ( tx1 * ( r43*c34 - c1345 )
     >   + ty1 * (     c34 - c1345 )
     >   + tz1 * (     c34 - c1345 ) )
               d(5,3,i,j) = dt * 2.0d+00 * tmp2 * u0(3)
     > * ( tx1 * ( c34 - c1345 )
     >   + ty1 * ( r43*c34 -c1345 )
     >   + tz1 * ( c34 - c1345 ) )
               d(5,4,i,j) = dt * 2.0d+00 * tmp2 * u0(4)
     > * ( tx1 * ( c34 - c1345 )
     >   + ty1 * ( c34 - c1345 )
     >   + tz1 * ( r43*c34 - c1345 ) )
               d(5,5,i,j) = 1.0d+00
     >   + dt * 2.0d+00 * ( tx1  + ty1 + tz1 ) * c1345 * tmp1
     >   + dt * 2.0d+00 * (  tx1 * dx5
     >                    +  ty1 * dy5
     >                    +  tz1 * dz5 )

c---------------------------------------------------------------------
c   form the first block sub-diagonal
c---------------------------------------------------------------------
c               tmp1 = rho_i(i,j,k-1)
               tmp1 = 1.0d+00 / uz(1)
               tmp2 = tmp1 * tmp1
               tmp3 = tmp1 * tmp2

               a(1,1,i,j) = - dt * tz1 * dz1
               a(1,2,i,j) =   0.0d+00
               a(1,3,i,j) =   0.0d+00
               a(1,4,i,j) = - dt * tz2
               a(1,5,i,j) =   0.0d+00

               a(2,1,i,j) = - dt * tz2
     >           * ( - ( uz(2)*uz(4) ) * tmp2 )
     >           - dt * tz1 * ( - c34 * tmp2 * uz(2) )
               a(2,2,i,j) = - dt * tz2 * ( uz(4) * tmp1 )
     >           - dt * tz1 * c34 * tmp1
     >           - dt * tz1 * dz2 
               a(2,3,i,j) = 0.0d+00
               a(2,4,i,j) = - dt * tz2 * ( uz(2) * tmp1 )
               a(2,5,i,j) = 0.0d+00

               a(3,1,i,j) = - dt * tz2
     >           * ( - ( uz(3)*uz(4) ) * tmp2 )
     >           - dt * tz1 * ( - c34 * tmp2 * uz(3) )
               a(3,2,i,j) = 0.0d+00
               a(3,3,i,j) = - dt * tz2 * ( uz(4) * tmp1 )
     >           - dt * tz1 * ( c34 * tmp1 )
     >           - dt * tz1 * dz3
               a(3,4,i,j) = - dt * tz2 * ( uz(3) * tmp1 )
               a(3,5,i,j) = 0.0d+00

               a(4,1,i,j) = - dt * tz2
     >        * ( - ( uz(4) * tmp1 ) ** 2
     >            + 0.50d+00 * c2
     >            * ( ( uz(2) * uz(2)
     >                + uz(3) * uz(3)
     >                + uz(4) * uz(4) ) * tmp2 ) )
     >        - dt * tz1 * ( - r43 * c34 * tmp2 * uz(4) )
               a(4,2,i,j) = - dt * tz2
     >             * ( - c2 * ( uz(2) * tmp1 ) )
               a(4,3,i,j) = - dt * tz2
     >             * ( - c2 * ( uz(3) * tmp1 ) )
               a(4,4,i,j) = - dt * tz2 * ( 2.0d+00 - c2 )
     >             * ( uz(4) * tmp1 )
     >             - dt * tz1 * ( r43 * c34 * tmp1 )
     >             - dt * tz1 * dz4
               a(4,5,i,j) = - dt * tz2 * c2

               a(5,1,i,j) = - dt * tz2
     >     * ( ( c2 * (  uz(2) * uz(2)
     >                 + uz(3) * uz(3)
     >                 + uz(4) * uz(4) ) * tmp2
     >       - c1 * ( uz(5) * tmp1 ) )
     >            * ( uz(4) * tmp1 ) )
     >       - dt * tz1
     >       * ( - ( c34 - c1345 ) * tmp3 * (uz(2)**2)
     >           - ( c34 - c1345 ) * tmp3 * (uz(3)**2)
     >           - ( r43*c34 - c1345 )* tmp3 * (uz(4)**2)
     >          - c1345 * tmp2 * uz(5) )
               a(5,2,i,j) = - dt * tz2
     >       * ( - c2 * ( uz(2)*uz(4) ) * tmp2 )
     >       - dt * tz1 * ( c34 - c1345 ) * tmp2 * uz(2)
               a(5,3,i,j) = - dt * tz2
     >       * ( - c2 * ( uz(3)*uz(4) ) * tmp2 )
     >       - dt * tz1 * ( c34 - c1345 ) * tmp2 * uz(3)
               a(5,4,i,j) = - dt * tz2
     >       * ( c1 * ( uz(5) * tmp1 )
     >       - 0.50d+00 * c2
     >       * ( (  uz(2)*uz(2)
     >            + uz(3)*uz(3)
     >            + 3.0d+00*uz(4)*uz(4) ) * tmp2 ) )
     >       - dt * tz1 * ( r43*c34 - c1345 ) * tmp2 * uz(4)
               a(5,5,i,j) = - dt * tz2
     >       * ( c1 * ( uz(4) * tmp1 ) )
     >       - dt * tz1 * c1345 * tmp1
     >       - dt * tz1 * dz5

c---------------------------------------------------------------------
c   form the second block sub-diagonal
c---------------------------------------------------------------------
c               tmp1 = rho_i(i,j-1,k)
               tmp1 = 1.0d+00 / uy(1)
               tmp2 = tmp1 * tmp1
               tmp3 = tmp1 * tmp2

               b(1,1,i,j) = - dt * ty1 * dy1
               b(1,2,i,j) =   0.0d+00
               b(1,3,i,j) = - dt * ty2
               b(1,4,i,j) =   0.0d+00
               b(1,5,i,j) =   0.0d+00

               b(2,1,i,j) = - dt * ty2
     >           * ( - ( uy(2)*uy(3) ) * tmp2 )
     >           - dt * ty1 * ( - c34 * tmp2 * uy(2) )
               b(2,2,i,j) = - dt * ty2 * ( uy(3) * tmp1 )
     >          - dt * ty1 * ( c34 * tmp1 )
     >          - dt * ty1 * dy2
               b(2,3,i,j) = - dt * ty2 * ( uy(2) * tmp1 )
               b(2,4,i,j) = 0.0d+00
               b(2,5,i,j) = 0.0d+00

               b(3,1,i,j) = - dt * ty2
     >           * ( - ( uy(3) * tmp1 ) ** 2
     >      + 0.50d+00 * c2 * ( (  uy(2) * uy(2)
     >                           + uy(3) * uy(3)
     >                           + uy(4) * uy(4) )
     >                          * tmp2 ) )
     >       - dt * ty1 * ( - r43 * c34 * tmp2 * uy(3) )
               b(3,2,i,j) = - dt * ty2
     >                   * ( - c2 * ( uy(2) * tmp1 ) )
               b(3,3,i,j) = - dt * ty2 * ( ( 2.0d+00 - c2 )
     >                   * ( uy(3) * tmp1 ) )
     >       - dt * ty1 * ( r43 * c34 * tmp1 )
     >       - dt * ty1 * dy3
               b(3,4,i,j) = - dt * ty2
     >                   * ( - c2 * ( uy(4) * tmp1 ) )
               b(3,5,i,j) = - dt * ty2 * c2

               b(4,1,i,j) = - dt * ty2
     >              * ( - ( uy(3)*uy(4) ) * tmp2 )
     >       - dt * ty1 * ( - c34 * tmp2 * uy(4) )
               b(4,2,i,j) = 0.0d+00
               b(4,3,i,j) = - dt * ty2 * ( uy(4) * tmp1 )
               b(4,4,i,j) = - dt * ty2 * ( uy(3) * tmp1 )
     >                        - dt * ty1 * ( c34 * tmp1 )
     >                        - dt * ty1 * dy4
               b(4,5,i,j) = 0.0d+00

               b(5,1,i,j) = - dt * ty2
     >          * ( ( c2 * (  uy(2) * uy(2)
     >                      + uy(3) * uy(3)
     >                      + uy(4) * uy(4) ) * tmp2
     >               - c1 * ( uy(5) * tmp1 ) )
     >          * ( uy(3) * tmp1 ) )
     >          - dt * ty1
     >          * ( - (     c34 - c1345 )*tmp3*(uy(2)**2)
     >              - ( r43*c34 - c1345 )*tmp3*(uy(3)**2)
     >              - (     c34 - c1345 )*tmp3*(uy(4)**2)
     >              - c1345*tmp2*uy(5) )
               b(5,2,i,j) = - dt * ty2
     >          * ( - c2 * ( uy(2)*uy(3) ) * tmp2 )
     >          - dt * ty1
     >          * ( c34 - c1345 ) * tmp2 * uy(2)
               b(5,3,i,j) = - dt * ty2
     >          * ( c1 * ( uy(5) * tmp1 )
     >          - 0.50d+00 * c2 
     >          * ( (  uy(2)*uy(2)
     >               + 3.0d+00 * uy(3)*uy(3)
     >               + uy(4)*uy(4) ) * tmp2 ) )
     >          - dt * ty1
     >          * ( r43*c34 - c1345 ) * tmp2 * uy(3)
               b(5,4,i,j) = - dt * ty2
     >          * ( - c2 * ( uy(3)*uy(4) ) * tmp2 )
     >          - dt * ty1 * ( c34 - c1345 ) * tmp2 * uy(4)
               b(5,5,i,j) = - dt * ty2
     >          * ( c1 * ( uy(3) * tmp1 ) )
     >          - dt * ty1 * c1345 * tmp1
     >          - dt * ty1 * dy5

c---------------------------------------------------------------------
c   form the third block sub-diagonal
c---------------------------------------------------------------------
c               tmp1 = rho_i(i-1,j,k)
               tmp1 = 1.0d+00 / ux(1)
               tmp2 = tmp1 * tmp1
               tmp3 = tmp1 * tmp2

               c(1,1,i,j) = - dt * tx1 * dx1
               c(1,2,i,j) = - dt * tx2
               c(1,3,i,j) =   0.0d+00
               c(1,4,i,j) =   0.0d+00
               c(1,5,i,j) =   0.0d+00

               c(2,1,i,j) = - dt * tx2
     >          * ( - ( ux(2) * tmp1 ) ** 2
     >     + c2 * 0.50d+00 * (  ux(2) * ux(2)
     >                        + ux(3) * ux(3)
     >                        + ux(4) * ux(4) ) * tmp2 )
     >          - dt * tx1 * ( - r43 * c34 * tmp2 * ux(2) )
               c(2,2,i,j) = - dt * tx2
     >          * ( ( 2.0d+00 - c2 ) * ( ux(2) * tmp1 ) )
     >          - dt * tx1 * ( r43 * c34 * tmp1 )
     >          - dt * tx1 * dx2
               c(2,3,i,j) = - dt * tx2
     >              * ( - c2 * ( ux(3) * tmp1 ) )
               c(2,4,i,j) = - dt * tx2
     >              * ( - c2 * ( ux(4) * tmp1 ) )
               c(2,5,i,j) = - dt * tx2 * c2 

               c(3,1,i,j) = - dt * tx2
     >              * ( - ( ux(2) * ux(3) ) * tmp2 )
     >         - dt * tx1 * ( - c34 * tmp2 * ux(3) )
               c(3,2,i,j) = - dt * tx2 * ( ux(3) * tmp1 )
               c(3,3,i,j) = - dt * tx2 * ( ux(2) * tmp1 )
     >          - dt * tx1 * ( c34 * tmp1 )
     >          - dt * tx1 * dx3
               c(3,4,i,j) = 0.0d+00
               c(3,5,i,j) = 0.0d+00

               c(4,1,i,j) = - dt * tx2
     >          * ( - ( ux(2)*ux(4) ) * tmp2 )
     >          - dt * tx1 * ( - c34 * tmp2 * ux(4) )
               c(4,2,i,j) = - dt * tx2 * ( ux(4) * tmp1 )
               c(4,3,i,j) = 0.0d+00
               c(4,4,i,j) = - dt * tx2 * ( ux(2) * tmp1 )
     >          - dt * tx1 * ( c34 * tmp1 )
     >          - dt * tx1 * dx4
               c(4,5,i,j) = 0.0d+00

               c(5,1,i,j) = - dt * tx2
     >          * ( ( c2 * (  ux(2) * ux(2)
     >                      + ux(3) * ux(3)
     >                      + ux(4) * ux(4) ) * tmp2
     >              - c1 * ( ux(5) * tmp1 ) )
     >          * ( ux(2) * tmp1 ) )
     >          - dt * tx1
     >          * ( - ( r43*c34 - c1345 ) * tmp3 * ( ux(2)**2 )
     >              - (     c34 - c1345 ) * tmp3 * ( ux(3)**2 )
     >              - (     c34 - c1345 ) * tmp3 * ( ux(4)**2 )
     >              - c1345 * tmp2 * ux(5) )
               c(5,2,i,j) = - dt * tx2
     >          * ( c1 * ( ux(5) * tmp1 )
     >             - 0.50d+00 * c2
     >             * ( (  3.0d+00*ux(2)*ux(2)
     >                  + ux(3)*ux(3)
     >                  + ux(4)*ux(4) ) * tmp2 ) )
     >           - dt * tx1
     >           * ( r43*c34 - c1345 ) * tmp2 * ux(2)
               c(5,3,i,j) = - dt * tx2
     >           * ( - c2 * ( ux(3)*ux(2) ) * tmp2 )
     >           - dt * tx1
     >           * (  c34 - c1345 ) * tmp2 * ux(3)
               c(5,4,i,j) = - dt * tx2
     >           * ( - c2 * ( ux(4)*ux(2) ) * tmp2 )
     >           - dt * tx1
     >           * (  c34 - c1345 ) * tmp2 * ux(4)
               c(5,5,i,j) = - dt * tx2
     >           * ( c1 * ( ux(2) * tmp1 ) )
     >           - dt * tx1 * c1345 * tmp1
     >           - dt * tx1 * dx5
        end do

      end do

      return
      end
