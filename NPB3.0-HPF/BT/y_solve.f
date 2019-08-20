c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine y_solve

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     Performs line solves in Y direction by first factoring
c     the block-tridiagonal matrix into an upper triangular matrix, 
c     and then performing back substitution to solve for the unknow
c     vectors of each line.  
c     
c     Make sure we treat elements zero to cell_size in the direction
c     of the sweep.
c---------------------------------------------------------------------

      include 'header.h'
      interface
        pure extrinsic (hpf_local) 
     >    subroutine matvec_sub(ablock,avec,bvec)
        double precision, dimension (:,:), intent(in):: ablock 
        double precision, dimension (:), intent(in) :: avec
        double precision, dimension (:), intent(inout) :: bvec
        end subroutine
	
        pure extrinsic (hpf_local)
     >    subroutine matmul_sub(ablock, bblock, cblock)
        implicit none
        double precision, dimension (:,:), intent(in):: ablock,bblock 
        double precision, dimension (:,:), intent(inout) :: cblock	
        end subroutine 

        pure extrinsic (hpf_local) subroutine binvcrhs( lhs,c,r )
          implicit none
          double precision, dimension (:,:), intent(inout):: lhs,c 
          double precision, dimension (:), intent(inout) :: r
        end subroutine      

        pure extrinsic (hpf_local) subroutine binvrhs(lhs,r)
        implicit none
        double precision, dimension (:,:), intent(inout):: lhs 
        double precision, dimension (:), intent(inout) :: r
        end subroutine 
	     
      end interface
      integer idx, jdx, kdx
      double precision pivot, coeff
      integer i, j, k, m, n, jsize

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      if (timeron) call timer_start(t_ysolve)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     This function computes the left hand side for the three y-factors   
c---------------------------------------------------------------------

      jsize = grid_points(2)-1

c---------------------------------------------------------------------
c     Compute the indices for storing the tri-diagonal matrix;
c     determine a (labeled f) and n jacobians for cell c
c---------------------------------------------------------------------
!HPF$ independent  new(tmp1, tmp2, tmp3, fjac, njac, lhs)
      do k = 1, grid_points(3)-2
         do i = 1, grid_points(1)-2
            do j = 0, jsize

               tmp1 = rho_i(i,j,k)
               tmp2 = tmp1 * tmp1
               tmp3 = tmp1 * tmp2

               fjac(1,1,j) = 0.0d+00
               fjac(1,2,j) = 0.0d+00
               fjac(1,3,j) = 1.0d+00
               fjac(1,4,j) = 0.0d+00
               fjac(1,5,j) = 0.0d+00

               fjac(2,1,j) = - ( u(2,i,j,k)*u(3,i,j,k) )
     >              * tmp2
               fjac(2,2,j) = u(3,i,j,k) * tmp1
               fjac(2,3,j) = u(2,i,j,k) * tmp1
               fjac(2,4,j) = 0.0d+00
               fjac(2,5,j) = 0.0d+00

               fjac(3,1,j) = - ( u(3,i,j,k)*u(3,i,j,k)*tmp2)
     >              + c2 * qs(i,j,k)
               fjac(3,2,j) = - c2 *  u(2,i,j,k) * tmp1
               fjac(3,3,j) = ( 2.0d+00 - c2 )
     >              *  u(3,i,j,k) * tmp1 
               fjac(3,4,j) = - c2 * u(4,i,j,k) * tmp1 
               fjac(3,5,j) = c2

               fjac(4,1,j) = - ( u(3,i,j,k)*u(4,i,j,k) )
     >              * tmp2
               fjac(4,2,j) = 0.0d+00
               fjac(4,3,j) = u(4,i,j,k) * tmp1
               fjac(4,4,j) = u(3,i,j,k) * tmp1
               fjac(4,5,j) = 0.0d+00

               fjac(5,1,j) = ( c2 * 2.0d0 * square(i,j,k)
     >              - c1 * u(5,i,j,k) )
     >              * u(3,i,j,k) * tmp2
               fjac(5,2,j) = - c2 * u(2,i,j,k)*u(3,i,j,k) 
     >              * tmp2
               fjac(5,3,j) = c1 * u(5,i,j,k) * tmp1 
     >              - c2 
     >              * ( qs(i,j,k)
     >              + u(3,i,j,k)*u(3,i,j,k) * tmp2 )
               fjac(5,4,j) = - c2 * ( u(3,i,j,k)*u(4,i,j,k) )
     >              * tmp2
               fjac(5,5,j) = c1 * u(3,i,j,k) * tmp1 

               njac(1,1,j) = 0.0d+00
               njac(1,2,j) = 0.0d+00
               njac(1,3,j) = 0.0d+00
               njac(1,4,j) = 0.0d+00
               njac(1,5,j) = 0.0d+00

               njac(2,1,j) = - c3c4 * tmp2 * u(2,i,j,k)
               njac(2,2,j) =   c3c4 * tmp1
               njac(2,3,j) =   0.0d+00
               njac(2,4,j) =   0.0d+00
               njac(2,5,j) =   0.0d+00

               njac(3,1,j) = - con43 * c3c4 * tmp2 * u(3,i,j,k)
               njac(3,2,j) =   0.0d+00
               njac(3,3,j) =   con43 * c3c4 * tmp1
               njac(3,4,j) =   0.0d+00
               njac(3,5,j) =   0.0d+00

               njac(4,1,j) = - c3c4 * tmp2 * u(4,i,j,k)
               njac(4,2,j) =   0.0d+00
               njac(4,3,j) =   0.0d+00
               njac(4,4,j) =   c3c4 * tmp1
               njac(4,5,j) =   0.0d+00

               njac(5,1,j) = - (  c3c4
     >              - c1345 ) * tmp3 * (u(2,i,j,k)**2)
     >              - ( con43 * c3c4
     >              - c1345 ) * tmp3 * (u(3,i,j,k)**2)
     >              - ( c3c4 - c1345 ) * tmp3 * (u(4,i,j,k)**2)
     >              - c1345 * tmp2 * u(5,i,j,k)

               njac(5,2,j) = (  c3c4 - c1345 ) * tmp2 * u(2,i,j,k)
               njac(5,3,j) = ( con43 * c3c4
     >              - c1345 ) * tmp2 * u(3,i,j,k)
               njac(5,4,j) = ( c3c4 - c1345 ) * tmp2 * u(4,i,j,k)
               njac(5,5,j) = ( c1345 ) * tmp1

            enddo

c---------------------------------------------------------------------
c     now joacobians set, so form left hand side in y direction
c---------------------------------------------------------------------
            lhs(:,:,:,0) = 0.0d0
            lhs(:,:,:,jsize) = 0.0d0
            lhs(1,1,2,0) = 1.0d0
            lhs(2,2,2,0) = 1.0d0
            lhs(3,3,2,0) = 1.0d0
            lhs(4,4,2,0) = 1.0d0
            lhs(5,5,2,0) = 1.0d0
            lhs(1,1,2,jsize) = 1.0d0
            lhs(2,2,2,jsize) = 1.0d0
            lhs(3,3,2,jsize) = 1.0d0
            lhs(4,4,2,jsize) = 1.0d0
            lhs(5,5,2,jsize) = 1.0d0

            do j = 1, jsize-1

               tmp1 = dt * ty1
               tmp2 = dt * ty2

               lhs(1,1,aa,j) = - tmp2 * fjac(1,1,j-1)
     >              - tmp1 * njac(1,1,j-1)
     >              - tmp1 * dy1 
               lhs(1,2,aa,j) = - tmp2 * fjac(1,2,j-1)
     >              - tmp1 * njac(1,2,j-1)
               lhs(1,3,aa,j) = - tmp2 * fjac(1,3,j-1)
     >              - tmp1 * njac(1,3,j-1)
               lhs(1,4,aa,j) = - tmp2 * fjac(1,4,j-1)
     >              - tmp1 * njac(1,4,j-1)
               lhs(1,5,aa,j) = - tmp2 * fjac(1,5,j-1)
     >              - tmp1 * njac(1,5,j-1)

               lhs(2,1,aa,j) = - tmp2 * fjac(2,1,j-1)
     >              - tmp1 * njac(2,1,j-1)
               lhs(2,2,aa,j) = - tmp2 * fjac(2,2,j-1)
     >              - tmp1 * njac(2,2,j-1)
     >              - tmp1 * dy2
               lhs(2,3,aa,j) = - tmp2 * fjac(2,3,j-1)
     >              - tmp1 * njac(2,3,j-1)
               lhs(2,4,aa,j) = - tmp2 * fjac(2,4,j-1)
     >              - tmp1 * njac(2,4,j-1)
               lhs(2,5,aa,j) = - tmp2 * fjac(2,5,j-1)
     >              - tmp1 * njac(2,5,j-1)

               lhs(3,1,aa,j) = - tmp2 * fjac(3,1,j-1)
     >              - tmp1 * njac(3,1,j-1)
               lhs(3,2,aa,j) = - tmp2 * fjac(3,2,j-1)
     >              - tmp1 * njac(3,2,j-1)
               lhs(3,3,aa,j) = - tmp2 * fjac(3,3,j-1)
     >              - tmp1 * njac(3,3,j-1)
     >              - tmp1 * dy3 
               lhs(3,4,aa,j) = - tmp2 * fjac(3,4,j-1)
     >              - tmp1 * njac(3,4,j-1)
               lhs(3,5,aa,j) = - tmp2 * fjac(3,5,j-1)
     >              - tmp1 * njac(3,5,j-1)

               lhs(4,1,aa,j) = - tmp2 * fjac(4,1,j-1)
     >              - tmp1 * njac(4,1,j-1)
               lhs(4,2,aa,j) = - tmp2 * fjac(4,2,j-1)
     >              - tmp1 * njac(4,2,j-1)
               lhs(4,3,aa,j) = - tmp2 * fjac(4,3,j-1)
     >              - tmp1 * njac(4,3,j-1)
               lhs(4,4,aa,j) = - tmp2 * fjac(4,4,j-1)
     >              - tmp1 * njac(4,4,j-1)
     >              - tmp1 * dy4
               lhs(4,5,aa,j) = - tmp2 * fjac(4,5,j-1)
     >              - tmp1 * njac(4,5,j-1)

               lhs(5,1,aa,j) = - tmp2 * fjac(5,1,j-1)
     >              - tmp1 * njac(5,1,j-1)
               lhs(5,2,aa,j) = - tmp2 * fjac(5,2,j-1)
     >              - tmp1 * njac(5,2,j-1)
               lhs(5,3,aa,j) = - tmp2 * fjac(5,3,j-1)
     >              - tmp1 * njac(5,3,j-1)
               lhs(5,4,aa,j) = - tmp2 * fjac(5,4,j-1)
     >              - tmp1 * njac(5,4,j-1)
               lhs(5,5,aa,j) = - tmp2 * fjac(5,5,j-1)
     >              - tmp1 * njac(5,5,j-1)
     >              - tmp1 * dy5

               lhs(1,1,bb,j) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(1,1,j)
     >              + tmp1 * 2.0d+00 * dy1
               lhs(1,2,bb,j) = tmp1 * 2.0d+00 * njac(1,2,j)
               lhs(1,3,bb,j) = tmp1 * 2.0d+00 * njac(1,3,j)
               lhs(1,4,bb,j) = tmp1 * 2.0d+00 * njac(1,4,j)
               lhs(1,5,bb,j) = tmp1 * 2.0d+00 * njac(1,5,j)

               lhs(2,1,bb,j) = tmp1 * 2.0d+00 * njac(2,1,j)
               lhs(2,2,bb,j) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(2,2,j)
     >              + tmp1 * 2.0d+00 * dy2
               lhs(2,3,bb,j) = tmp1 * 2.0d+00 * njac(2,3,j)
               lhs(2,4,bb,j) = tmp1 * 2.0d+00 * njac(2,4,j)
               lhs(2,5,bb,j) = tmp1 * 2.0d+00 * njac(2,5,j)

               lhs(3,1,bb,j) = tmp1 * 2.0d+00 * njac(3,1,j)
               lhs(3,2,bb,j) = tmp1 * 2.0d+00 * njac(3,2,j)
               lhs(3,3,bb,j) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(3,3,j)
     >              + tmp1 * 2.0d+00 * dy3
               lhs(3,4,bb,j) = tmp1 * 2.0d+00 * njac(3,4,j)
               lhs(3,5,bb,j) = tmp1 * 2.0d+00 * njac(3,5,j)

               lhs(4,1,bb,j) = tmp1 * 2.0d+00 * njac(4,1,j)
               lhs(4,2,bb,j) = tmp1 * 2.0d+00 * njac(4,2,j)
               lhs(4,3,bb,j) = tmp1 * 2.0d+00 * njac(4,3,j)
               lhs(4,4,bb,j) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(4,4,j)
     >              + tmp1 * 2.0d+00 * dy4
               lhs(4,5,bb,j) = tmp1 * 2.0d+00 * njac(4,5,j)

               lhs(5,1,bb,j) = tmp1 * 2.0d+00 * njac(5,1,j)
               lhs(5,2,bb,j) = tmp1 * 2.0d+00 * njac(5,2,j)
               lhs(5,3,bb,j) = tmp1 * 2.0d+00 * njac(5,3,j)
               lhs(5,4,bb,j) = tmp1 * 2.0d+00 * njac(5,4,j)
               lhs(5,5,bb,j) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(5,5,j) 
     >              + tmp1 * 2.0d+00 * dy5

               lhs(1,1,cc,j) =  tmp2 * fjac(1,1,j+1)
     >              - tmp1 * njac(1,1,j+1)
     >              - tmp1 * dy1
               lhs(1,2,cc,j) =  tmp2 * fjac(1,2,j+1)
     >              - tmp1 * njac(1,2,j+1)
               lhs(1,3,cc,j) =  tmp2 * fjac(1,3,j+1)
     >              - tmp1 * njac(1,3,j+1)
               lhs(1,4,cc,j) =  tmp2 * fjac(1,4,j+1)
     >              - tmp1 * njac(1,4,j+1)
               lhs(1,5,cc,j) =  tmp2 * fjac(1,5,j+1)
     >              - tmp1 * njac(1,5,j+1)

               lhs(2,1,cc,j) =  tmp2 * fjac(2,1,j+1)
     >              - tmp1 * njac(2,1,j+1)
               lhs(2,2,cc,j) =  tmp2 * fjac(2,2,j+1)
     >              - tmp1 * njac(2,2,j+1)
     >              - tmp1 * dy2
               lhs(2,3,cc,j) =  tmp2 * fjac(2,3,j+1)
     >              - tmp1 * njac(2,3,j+1)
               lhs(2,4,cc,j) =  tmp2 * fjac(2,4,j+1)
     >              - tmp1 * njac(2,4,j+1)
               lhs(2,5,cc,j) =  tmp2 * fjac(2,5,j+1)
     >              - tmp1 * njac(2,5,j+1)

               lhs(3,1,cc,j) =  tmp2 * fjac(3,1,j+1)
     >              - tmp1 * njac(3,1,j+1)
               lhs(3,2,cc,j) =  tmp2 * fjac(3,2,j+1)
     >              - tmp1 * njac(3,2,j+1)
               lhs(3,3,cc,j) =  tmp2 * fjac(3,3,j+1)
     >              - tmp1 * njac(3,3,j+1)
     >              - tmp1 * dy3
               lhs(3,4,cc,j) =  tmp2 * fjac(3,4,j+1)
     >              - tmp1 * njac(3,4,j+1)
               lhs(3,5,cc,j) =  tmp2 * fjac(3,5,j+1)
     >              - tmp1 * njac(3,5,j+1)

               lhs(4,1,cc,j) =  tmp2 * fjac(4,1,j+1)
     >              - tmp1 * njac(4,1,j+1)
               lhs(4,2,cc,j) =  tmp2 * fjac(4,2,j+1)
     >              - tmp1 * njac(4,2,j+1)
               lhs(4,3,cc,j) =  tmp2 * fjac(4,3,j+1)
     >              - tmp1 * njac(4,3,j+1)
               lhs(4,4,cc,j) =  tmp2 * fjac(4,4,j+1)
     >              - tmp1 * njac(4,4,j+1)
     >              - tmp1 * dy4
               lhs(4,5,cc,j) =  tmp2 * fjac(4,5,j+1)
     >              - tmp1 * njac(4,5,j+1)

               lhs(5,1,cc,j) =  tmp2 * fjac(5,1,j+1)
     >              - tmp1 * njac(5,1,j+1)
               lhs(5,2,cc,j) =  tmp2 * fjac(5,2,j+1)
     >              - tmp1 * njac(5,2,j+1)
               lhs(5,3,cc,j) =  tmp2 * fjac(5,3,j+1)
     >              - tmp1 * njac(5,3,j+1)
               lhs(5,4,cc,j) =  tmp2 * fjac(5,4,j+1)
     >              - tmp1 * njac(5,4,j+1)
               lhs(5,5,cc,j) =  tmp2 * fjac(5,5,j+1)
     >              - tmp1 * njac(5,5,j+1)
     >              - tmp1 * dy5

            enddo

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     performs guaussian elimination on this cell.
c     
c     assumes that unpacking routines for non-first cells 
c     preload C' and rhs' from previous cell.
c     
c     assumed send happens outside this routine, but that
c     c'(JMAX) and rhs'(JMAX) will be sent to next cell
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     multiply c(i,0,k) by b_inverse and copy back to c
c     multiply rhs(0) by b_inverse(0) and copy to rhs
c---------------------------------------------------------------------
c            call binvcrhs( lhs(:,:,bb,0),
c     >                        lhs(:,:,cc,0),
c     >                        rhs(:,i,0,k) )
           do kdx = 1, 5           
             pivot = 1.00d0/lhs(kdx,kdx,bb,0)
             rhs(kdx,i,0,k) = rhs(kdx,i,0,k)*pivot
             do jdx = kdx+1, 5
               lhs(kdx,jdx,bb,0) = lhs(kdx,jdx,bb,0)*pivot
             end do
             do jdx = 1, 5
               lhs(kdx,jdx,cc,0) = lhs(kdx,jdx,cc,0)*pivot
             end do
             
             
             do jdx = 1, kdx-1
               coeff = lhs(jdx,kdx,bb,0)
               rhs(jdx,i,0,k) = rhs(jdx,i,0,k) - rhs(kdx,i,0,k)*coeff
               do idx = kdx+1, 5
                 lhs(jdx,idx,bb,0) = lhs(jdx,idx,bb,0) - 
     >               lhs(kdx,idx,bb,0)*coeff
               end do 
               do idx = 1, 5
                 lhs(jdx,idx,cc,0) = lhs(jdx,idx,cc,0) - 
     >               lhs(kdx,idx,cc,0)*coeff
               end do 

             end do 

             do jdx = kdx+1,5
               coeff = lhs(jdx,kdx,bb,0)
               rhs(jdx,i,0,k) = rhs(jdx,i,0,k) - rhs(kdx,i,0,k)*coeff
               do idx = kdx+1, 5
                 lhs(jdx,idx,bb,0) = lhs(jdx,idx,bb,0) - 
     >               lhs(kdx,idx,bb,0)*coeff
               end do 
               do idx = 1, 5
                 lhs(jdx,idx,cc,0) = lhs(jdx,idx,cc,0) - 
     >               lhs(kdx,idx,cc,0)*coeff
               end do 
             end do                           
           end do 
c---------------------------------------------------------------------
c     begin inner most do loop
c     do all the elements of the cell unless last 
c---------------------------------------------------------------------
            do j=1,jsize-1

c---------------------------------------------------------------------
c     subtract A*lhs_vector(j-1) from lhs_vector(j)
c     
c     rhs(j) = rhs(j) - A*rhs(j-1)
c---------------------------------------------------------------------
c               call matvec_sub(lhs(:,:,aa,j),
c     >                         rhs(:,i,j-1,k),rhs(:,i,j,k))

c---------------------------------------------------------------------
c     B(j) = B(j) - C(j-1)*A(j)
c---------------------------------------------------------------------
c               call matmul_sub(lhs(:,:,aa,j),
c     >                         lhs(:,:,cc,j-1),
c     >                         lhs(:,:,bb,j))

           do kdx = 1, 5           
             do jdx = 1, 5
              rhs(kdx,i,j,k) = rhs(kdx,i,j,k) -                   
     >                     lhs(kdx,jdx,aa,j)*rhs(jdx,i,j-1,k)
              lhs(kdx,jdx,bb,j) = lhs(kdx,jdx,bb,j) -     
     >               lhs(kdx,1,aa,j)*lhs(1,jdx,cc,j-1) -
     >               lhs(kdx,2,aa,j)*lhs(2,jdx,cc,j-1) -
     >               lhs(kdx,3,aa,j)*lhs(3,jdx,cc,j-1) -
     >               lhs(kdx,4,aa,j)*lhs(4,jdx,cc,j-1) -
     >               lhs(kdx,5,aa,j)*lhs(5,jdx,cc,j-1)    
             end do                                                
           end do 
c---------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(i,1,k) by b_inverse(i,1,k) and copy to rhs
c---------------------------------------------------------------------
c               call binvcrhs( lhs(:,:,bb,j),
c     >                        lhs(:,:,cc,j),
c     >                        rhs(:,i,j,k) )
           do kdx = 1, 5           
             pivot = 1.00d0/lhs(kdx,kdx,bb,j)
             rhs(kdx,i,j,k) = rhs(kdx,i,j,k)*pivot
             do jdx = kdx+1, 5
               lhs(kdx,jdx,bb,j) = lhs(kdx,jdx,bb,j)*pivot
             end do
             do jdx = 1, 5
               lhs(kdx,jdx,cc,j) = lhs(kdx,jdx,cc,j)*pivot   
             end do                                                  
             
             
             do jdx = 1, kdx-1
               coeff = lhs(jdx,kdx,bb,j)
               rhs(jdx,i,j,k) = rhs(jdx,i,j,k) - rhs(kdx,i,j,k)*coeff
               do idx = kdx+1, 5
                 lhs(jdx,idx,bb,j) = lhs(jdx,idx,bb,j) - 
     >                               lhs(kdx,idx,bb,j)*coeff
               end do 
                 lhs(jdx,1,cc,j) = lhs(jdx,1,cc,j) - 
     >                             lhs(kdx,1,cc,j)*coeff
                 lhs(jdx,2,cc,j) = lhs(jdx,2,cc,j) - 
     >                             lhs(kdx,2,cc,j)*coeff
                 lhs(jdx,3,cc,j) = lhs(jdx,3,cc,j) - 
     >                             lhs(kdx,3,cc,j)*coeff
                 lhs(jdx,4,cc,j) = lhs(jdx,4,cc,j) - 
     >                             lhs(kdx,4,cc,j)*coeff
                 lhs(jdx,5,cc,j) = lhs(jdx,5,cc,j) - 
     >                             lhs(kdx,5,cc,j)*coeff
             end do 

             do jdx = kdx+1,5
               coeff = lhs(jdx,kdx,bb,j)
               rhs(jdx,i,j,k) = rhs(jdx,i,j,k) - rhs(kdx,i,j,k)*coeff
               do idx = kdx+1, 5
                 lhs(jdx,idx,bb,j) = lhs(jdx,idx,bb,j) -    
     >                               lhs(kdx,idx,bb,j)*coeff
               end do                                               
                 lhs(jdx,1,cc,j) = lhs(jdx,1,cc,j) - 
     >                             lhs(kdx,1,cc,j)*coeff
                 lhs(jdx,2,cc,j) = lhs(jdx,2,cc,j) - 
     >                             lhs(kdx,2,cc,j)*coeff
                 lhs(jdx,3,cc,j) = lhs(jdx,3,cc,j) - 
     >                             lhs(kdx,3,cc,j)*coeff
                 lhs(jdx,4,cc,j) = lhs(jdx,4,cc,j) - 
     >                             lhs(kdx,4,cc,j)*coeff
                 lhs(jdx,5,cc,j) = lhs(jdx,5,cc,j) - 
     >                             lhs(kdx,5,cc,j)*coeff
             end do                          
           end do 

            end do


c---------------------------------------------------------------------
c     rhs(jsize) = rhs(jsize) - A*rhs(jsize-1)
c---------------------------------------------------------------------
c            call matvec_sub(lhs(:,:,aa,jsize),
c     >                         rhs(:,i,jsize-1,k),rhs(:,i,jsize,k))
c
c---------------------------------------------------------------------
c     B(jsize) = B(jsize) - C(jsize-1)*A(jsize)
c     call matmul_sub(aa,i,jsize,k,c,
c     $              cc,i,jsize-1,k,c,bb,i,jsize,k)
c---------------------------------------------------------------------
c            call matmul_sub(lhs(:,:,aa,jsize),
c     >                         lhs(:,:,cc,jsize-1),
c     >                         lhs(:,:,bb,jsize))
           do kdx = 1, 5           
             do jdx = 1, 5
               rhs(kdx,i,jsize,k) = rhs(kdx,i,jsize,k) - 
     >              lhs(kdx,jdx,aa,jsize)*rhs(jdx,i,jsize-1,k)
     
               lhs(kdx,jdx,bb,jsize) = lhs(kdx,jdx,bb,jsize) - 
     >            lhs(kdx,1,aa,jsize)*lhs(1,jdx,cc,jsize-1) -
     >            lhs(kdx,2,aa,jsize)*lhs(2,jdx,cc,jsize-1) -
     >            lhs(kdx,3,aa,jsize)*lhs(3,jdx,cc,jsize-1) -
     >            lhs(kdx,4,aa,jsize)*lhs(4,jdx,cc,jsize-1) -
     >            lhs(kdx,5,aa,jsize)*lhs(5,jdx,cc,jsize-1)
             end do 
           end do 

c---------------------------------------------------------------------
c     multiply rhs(jsize) by b_inverse(jsize) and copy to rhs
c---------------------------------------------------------------------
c            call binvrhs( lhs(:,:,bb,jsize),
c     >                       rhs(:,i,jsize,k) )
          do kdx = 1, 5 
             pivot = 1.00d0/lhs(kdx,kdx,bb,jsize)
             rhs(kdx,i,jsize,k) = rhs(kdx,i,jsize,k)*pivot
             do jdx = kdx+1, 5
               lhs(kdx,jdx,bb,jsize) = lhs(kdx,jdx,bb,jsize)*pivot
             end do       
             
             do jdx = 1, kdx-1
               coeff = lhs(jdx,kdx,bb,jsize)
               rhs(jdx,i,jsize,k) = rhs(jdx,i,jsize,k) - 
     >                              rhs(kdx,i,jsize,k)*coeff
               do idx = kdx+1, 5
                 lhs(jdx,idx,bb,jsize) = lhs(jdx,idx,bb,jsize) - 
     >                                   lhs(kdx,idx,bb,jsize)*coeff
               end do 
             end do 

             do jdx = kdx+1,5
               coeff = lhs(jdx,kdx,bb,jsize)
               rhs(jdx,i,jsize,k) = rhs(jdx,i,jsize,k) - 
     >                              rhs(kdx,i,jsize,k)*coeff
               do idx = kdx+1, 5
                 lhs(jdx,idx,bb,jsize) = lhs(jdx,idx,bb,jsize) - 
     >                                   lhs(kdx,idx,bb,jsize)*coeff
               end do 
             end do                           
           end do  


c---------------------------------------------------------------------
c     back solve: if last cell, then generate U(jsize)=rhs(jsize)
c     else assume U(jsize) is loaded in un pack backsub_info
c     so just use it
c     after call u(jstart) will be sent to next cell
c---------------------------------------------------------------------
      
            do j=jsize-1,0,-1
               do m=1,BLOCK_SIZE
                  do n=1,BLOCK_SIZE
                     rhs(m,i,j,k) = rhs(m,i,j,k) 
     >                    - lhs(m,n,cc,j)*rhs(n,i,j+1,k)
                  enddo
               enddo
            enddo

         enddo
      enddo
      if (timeron) call timer_stop(t_ysolve)

      return
      end


