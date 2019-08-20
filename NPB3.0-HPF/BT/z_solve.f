c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine z_solve

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     Performs line solves in Z direction by first factoring
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

      double precision rhsz(5,0:IMAX/2*2, 0:JMAX-1,0:KMAX-1),
     >                   uz(5, 0:IMAX/2*2,0:JMAX-1,0:KMAX-1)

!HPF$    template gridz(0:JMAX-1)
!HPF$    distribute(block) :: gridz
!HPF$    align(*,*,:,*) with gridz :: rhsz,uz

      integer i, j, k, m, n, ksize

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      if (timeron) call timer_start(t_rdis1)
      uz(:,:,0:JMAX-1,:)=u(:,:,0:JMAX-1,:)
      rhsz(:,:,0:JMAX-1,:)=rhs(:,:,0:JMAX-1,:)
      if (timeron) call timer_stop(t_rdis1)

      if (timeron) call timer_start(t_zsolve)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     This function computes the left hand side for the three z-factors
c---------------------------------------------------------------------

      ksize = grid_points(3)-1
c---------------------------------------------------------------------
c     Compute the indices for storing the block-diagonal matrix;
c     determine c (labeled f) and s jacobians
c---------------------------------------------------------------------
!HPF$ independent  new(tmp1, tmp2, tmp3, njac, fjac, lhs)
      do j = 1, grid_points(2)-2
         do i = 1, grid_points(1)-2
      	    do k = 0, ksize
               tmp1 = 1.0d+00 / uz(1,i,j,k)
               tmp2 = tmp1 * tmp1
               tmp3 = tmp1 * tmp2

               fjac(1,1,k) = 0.0d+00
               fjac(1,2,k) = 0.0d+00
               fjac(1,3,k) = 0.0d+00
               fjac(1,4,k) = 1.0d+00
               fjac(1,5,k) = 0.0d+00
               fjac(2,1,k) = - ( uz(2,i,j,k)*uz(4,i,j,k) )
     >              * tmp2
               fjac(2,2,k) = uz(4,i,j,k) * tmp1
               fjac(2,3,k) = 0.0d+00
               fjac(2,4,k) = uz(2,i,j,k) * tmp1
               fjac(2,5,k) = 0.0d+00

               fjac(3,1,k) = - ( uz(3,i,j,k)*uz(4,i,j,k) )
     >              * tmp2
               fjac(3,2,k) = 0.0d+00
               fjac(3,3,k) = uz(4,i,j,k) * tmp1
               fjac(3,4,k) = uz(3,i,j,k) * tmp1
               fjac(3,5,k) = 0.0d+00

c               fjac(4,1,k) = - (uz(4,i,j,k)*uz(4,i,j,k) * tmp2 )
c     >              + c2 * qs(i,j,k)

               fjac(4,1,k) = - (uz(4,i,j,k)*uz(4,i,j,k) * tmp2 )
     >              + 0.50d+00 * c2 * ( (  uz(2,i,j,k) * uz(2,i,j,k)
     >              + uz(3,i,j,k) * uz(3,i,j,k)
     >              + uz(4,i,j,k) * uz(4,i,j,k) ) * tmp2 )
               fjac(4,2,k) = - c2 *  uz(2,i,j,k) * tmp1
               fjac(4,3,k) = - c2 *  uz(3,i,j,k) * tmp1
               fjac(4,4,k) = ( 2.0d+00 - c2 )
     >              *  uz(4,i,j,k) * tmp1
               fjac(4,5,k) = c2

c               fjac(5,1,k) = ( c2 * 2.0d0 * square(i,j,k)
c     >              - c1 * uz(5,i,j,k) )
c     >              * uz(4,i,j,k) * tmp2
               fjac(5,1,k) = ( c2 * (  uz(2,i,j,k) * uz(2,i,j,k)
     >              + uz(3,i,j,k) * uz(3,i,j,k)
     >              + uz(4,i,j,k) * uz(4,i,j,k) )
     >              * tmp2
     >              - c1 * ( uz(5,i,j,k) * tmp1 ) )
     >              * ( uz(4,i,j,k) * tmp1 )
               fjac(5,2,k) = - c2 * ( uz(2,i,j,k)*uz(4,i,j,k) )
     >              * tmp2
               fjac(5,3,k) = - c2 * ( uz(3,i,j,k)*uz(4,i,j,k) )
     >              * tmp2
c               fjac(5,4,k) = c1 * ( uz(5,i,j,k) * tmp1 )
c     >              - c2
c     >              * ( qs(i,j,k)
c     >              + uz(4,i,j,k)*uz(4,i,j,k) * tmp2 )
               fjac(5,4,k) = c1 * ( uz(5,i,j,k) * tmp1 )
     >              - 0.50d+00 * c2
     >              * ( (  uz(2,i,j,k)*uz(2,i,j,k)
     >              + uz(3,i,j,k)*uz(3,i,j,k)
     >              + 3.0d+00*uz(4,i,j,k)*uz(4,i,j,k) )
     >              * tmp2 )
               fjac(5,5,k) = c1 * uz(4,i,j,k) * tmp1

               njac(1,1,k) = 0.0d+00
               njac(1,2,k) = 0.0d+00
               njac(1,3,k) = 0.0d+00
               njac(1,4,k) = 0.0d+00
               njac(1,5,k) = 0.0d+00

               njac(2,1,k) = - c3c4 * tmp2 * uz(2,i,j,k)
               njac(2,2,k) =   c3c4 * tmp1
               njac(2,3,k) =   0.0d+00
               njac(2,4,k) =   0.0d+00
               njac(2,5,k) =   0.0d+00

               njac(3,1,k) = - c3c4 * tmp2 * uz(3,i,j,k)
               njac(3,2,k) =   0.0d+00
               njac(3,3,k) =   c3c4 * tmp1
               njac(3,4,k) =   0.0d+00
               njac(3,5,k) =   0.0d+00

               njac(4,1,k) = - con43 * c3c4 * tmp2 * uz(4,i,j,k)
               njac(4,2,k) =   0.0d+00
               njac(4,3,k) =   0.0d+00
               njac(4,4,k) =   con43 * c3 * c4 * tmp1
               njac(4,5,k) =   0.0d+00

               njac(5,1,k) = - (  c3c4
     >              - c1345 ) * tmp3 * (uz(2,i,j,k)**2)
     >              - ( c3c4 - c1345 ) * tmp3 * (uz(3,i,j,k)**2)
     >              - ( con43 * c3c4
     >              - c1345 ) * tmp3 * (uz(4,i,j,k)**2)
     >              - c1345 * tmp2 * uz(5,i,j,k)

               njac(5,2,k) = (  c3c4 - c1345 ) * tmp2 * uz(2,i,j,k)
               njac(5,3,k) = (  c3c4 - c1345 ) * tmp2 * uz(3,i,j,k)
               njac(5,4,k) = ( con43 * c3c4
     >              - c1345 ) * tmp2 * uz(4,i,j,k)
               njac(5,5,k) = ( c1345 )* tmp1

            enddo

c---------------------------------------------------------------------
c     now jacobians set, so form left hand side in z direction
c---------------------------------------------------------------------
            lhs(:,:,:,0) = 0.0d0
            lhs(:,:,:,ksize) = 0.0d0
            lhs(1,1,2,0) = 1.0d0
            lhs(2,2,2,0) = 1.0d0
            lhs(3,3,2,0) = 1.0d0
            lhs(4,4,2,0) = 1.0d0
            lhs(5,5,2,0) = 1.0d0
            lhs(1,1,2,ksize) = 1.0d0
            lhs(2,2,2,ksize) = 1.0d0
            lhs(3,3,2,ksize) = 1.0d0
            lhs(4,4,2,ksize) = 1.0d0
            lhs(5,5,2,ksize) = 1.0d0

      	    do k = 1, ksize-1

               tmp1 = dt * tz1
               tmp2 = dt * tz2

               lhs(1,1,aa,k) = - tmp2 * fjac(1,1,k-1)
     >              - tmp1 * njac(1,1,k-1)
     >              - tmp1 * dz1
               lhs(1,2,aa,k) = - tmp2 * fjac(1,2,k-1)
     >              - tmp1 * njac(1,2,k-1)
               lhs(1,3,aa,k) = - tmp2 * fjac(1,3,k-1)
     >              - tmp1 * njac(1,3,k-1)
               lhs(1,4,aa,k) = - tmp2 * fjac(1,4,k-1)
     >              - tmp1 * njac(1,4,k-1)
               lhs(1,5,aa,k) = - tmp2 * fjac(1,5,k-1)
     >              - tmp1 * njac(1,5,k-1)

               lhs(2,1,aa,k) = - tmp2 * fjac(2,1,k-1)
     >              - tmp1 * njac(2,1,k-1)
               lhs(2,2,aa,k) = - tmp2 * fjac(2,2,k-1)
     >              - tmp1 * njac(2,2,k-1)
     >              - tmp1 * dz2
               lhs(2,3,aa,k) = - tmp2 * fjac(2,3,k-1)
     >              - tmp1 * njac(2,3,k-1)
               lhs(2,4,aa,k) = - tmp2 * fjac(2,4,k-1)
     >              - tmp1 * njac(2,4,k-1)
               lhs(2,5,aa,k) = - tmp2 * fjac(2,5,k-1)
     >              - tmp1 * njac(2,5,k-1)

               lhs(3,1,aa,k) = - tmp2 * fjac(3,1,k-1)
     >              - tmp1 * njac(3,1,k-1)
               lhs(3,2,aa,k) = - tmp2 * fjac(3,2,k-1)
     >              - tmp1 * njac(3,2,k-1)
               lhs(3,3,aa,k) = - tmp2 * fjac(3,3,k-1)
     >              - tmp1 * njac(3,3,k-1)
     >              - tmp1 * dz3
               lhs(3,4,aa,k) = - tmp2 * fjac(3,4,k-1)
     >              - tmp1 * njac(3,4,k-1)
               lhs(3,5,aa,k) = - tmp2 * fjac(3,5,k-1)
     >              - tmp1 * njac(3,5,k-1)

               lhs(4,1,aa,k) = - tmp2 * fjac(4,1,k-1)
     >              - tmp1 * njac(4,1,k-1)
               lhs(4,2,aa,k) = - tmp2 * fjac(4,2,k-1)
     >              - tmp1 * njac(4,2,k-1)
               lhs(4,3,aa,k) = - tmp2 * fjac(4,3,k-1)
     >              - tmp1 * njac(4,3,k-1)
               lhs(4,4,aa,k) = - tmp2 * fjac(4,4,k-1)
     >              - tmp1 * njac(4,4,k-1)
     >              - tmp1 * dz4
               lhs(4,5,aa,k) = - tmp2 * fjac(4,5,k-1)
     >              - tmp1 * njac(4,5,k-1)

               lhs(5,1,aa,k) = - tmp2 * fjac(5,1,k-1)
     >              - tmp1 * njac(5,1,k-1)
               lhs(5,2,aa,k) = - tmp2 * fjac(5,2,k-1)
     >              - tmp1 * njac(5,2,k-1)
               lhs(5,3,aa,k) = - tmp2 * fjac(5,3,k-1)
     >              - tmp1 * njac(5,3,k-1)
               lhs(5,4,aa,k) = - tmp2 * fjac(5,4,k-1)
     >              - tmp1 * njac(5,4,k-1)
               lhs(5,5,aa,k) = - tmp2 * fjac(5,5,k-1)
     >              - tmp1 * njac(5,5,k-1)
     >              - tmp1 * dz5

               lhs(1,1,bb,k) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(1,1,k)
     >              + tmp1 * 2.0d+00 * dz1
               lhs(1,2,bb,k) = tmp1 * 2.0d+00 * njac(1,2,k)
               lhs(1,3,bb,k) = tmp1 * 2.0d+00 * njac(1,3,k)
               lhs(1,4,bb,k) = tmp1 * 2.0d+00 * njac(1,4,k)
               lhs(1,5,bb,k) = tmp1 * 2.0d+00 * njac(1,5,k)

               lhs(2,1,bb,k) = tmp1 * 2.0d+00 * njac(2,1,k)
               lhs(2,2,bb,k) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(2,2,k)
     >              + tmp1 * 2.0d+00 * dz2
               lhs(2,3,bb,k) = tmp1 * 2.0d+00 * njac(2,3,k)
               lhs(2,4,bb,k) = tmp1 * 2.0d+00 * njac(2,4,k)
               lhs(2,5,bb,k) = tmp1 * 2.0d+00 * njac(2,5,k)

               lhs(3,1,bb,k) = tmp1 * 2.0d+00 * njac(3,1,k)
               lhs(3,2,bb,k) = tmp1 * 2.0d+00 * njac(3,2,k)
               lhs(3,3,bb,k) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(3,3,k)
     >              + tmp1 * 2.0d+00 * dz3
               lhs(3,4,bb,k) = tmp1 * 2.0d+00 * njac(3,4,k)
               lhs(3,5,bb,k) = tmp1 * 2.0d+00 * njac(3,5,k)

               lhs(4,1,bb,k) = tmp1 * 2.0d+00 * njac(4,1,k)
               lhs(4,2,bb,k) = tmp1 * 2.0d+00 * njac(4,2,k)
               lhs(4,3,bb,k) = tmp1 * 2.0d+00 * njac(4,3,k)
               lhs(4,4,bb,k) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(4,4,k)
     >              + tmp1 * 2.0d+00 * dz4
               lhs(4,5,bb,k) = tmp1 * 2.0d+00 * njac(4,5,k)

               lhs(5,1,bb,k) = tmp1 * 2.0d+00 * njac(5,1,k)
               lhs(5,2,bb,k) = tmp1 * 2.0d+00 * njac(5,2,k)
               lhs(5,3,bb,k) = tmp1 * 2.0d+00 * njac(5,3,k)
               lhs(5,4,bb,k) = tmp1 * 2.0d+00 * njac(5,4,k)
               lhs(5,5,bb,k) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(5,5,k)
     >              + tmp1 * 2.0d+00 * dz5

               lhs(1,1,cc,k) =  tmp2 * fjac(1,1,k+1)
     >              - tmp1 * njac(1,1,k+1)
     >              - tmp1 * dz1
               lhs(1,2,cc,k) =  tmp2 * fjac(1,2,k+1)
     >              - tmp1 * njac(1,2,k+1)
               lhs(1,3,cc,k) =  tmp2 * fjac(1,3,k+1)
     >              - tmp1 * njac(1,3,k+1)
               lhs(1,4,cc,k) =  tmp2 * fjac(1,4,k+1)
     >              - tmp1 * njac(1,4,k+1)
               lhs(1,5,cc,k) =  tmp2 * fjac(1,5,k+1)
     >              - tmp1 * njac(1,5,k+1)

               lhs(2,1,cc,k) =  tmp2 * fjac(2,1,k+1)
     >              - tmp1 * njac(2,1,k+1)
               lhs(2,2,cc,k) =  tmp2 * fjac(2,2,k+1)
     >              - tmp1 * njac(2,2,k+1)
     >              - tmp1 * dz2
               lhs(2,3,cc,k) =  tmp2 * fjac(2,3,k+1)
     >              - tmp1 * njac(2,3,k+1)
               lhs(2,4,cc,k) =  tmp2 * fjac(2,4,k+1)
     >              - tmp1 * njac(2,4,k+1)
               lhs(2,5,cc,k) =  tmp2 * fjac(2,5,k+1)
     >              - tmp1 * njac(2,5,k+1)

               lhs(3,1,cc,k) =  tmp2 * fjac(3,1,k+1)
     >              - tmp1 * njac(3,1,k+1)
               lhs(3,2,cc,k) =  tmp2 * fjac(3,2,k+1)
     >              - tmp1 * njac(3,2,k+1)
               lhs(3,3,cc,k) =  tmp2 * fjac(3,3,k+1)
     >              - tmp1 * njac(3,3,k+1)
     >              - tmp1 * dz3
               lhs(3,4,cc,k) =  tmp2 * fjac(3,4,k+1)
     >              - tmp1 * njac(3,4,k+1)
               lhs(3,5,cc,k) =  tmp2 * fjac(3,5,k+1)
     >              - tmp1 * njac(3,5,k+1)

               lhs(4,1,cc,k) =  tmp2 * fjac(4,1,k+1)
     >              - tmp1 * njac(4,1,k+1)
               lhs(4,2,cc,k) =  tmp2 * fjac(4,2,k+1)
     >              - tmp1 * njac(4,2,k+1)
               lhs(4,3,cc,k) =  tmp2 * fjac(4,3,k+1)
     >              - tmp1 * njac(4,3,k+1)
               lhs(4,4,cc,k) =  tmp2 * fjac(4,4,k+1)
     >              - tmp1 * njac(4,4,k+1)
     >              - tmp1 * dz4
               lhs(4,5,cc,k) =  tmp2 * fjac(4,5,k+1)
     >              - tmp1 * njac(4,5,k+1)

               lhs(5,1,cc,k) =  tmp2 * fjac(5,1,k+1)
     >              - tmp1 * njac(5,1,k+1)
               lhs(5,2,cc,k) =  tmp2 * fjac(5,2,k+1)
     >              - tmp1 * njac(5,2,k+1)
               lhs(5,3,cc,k) =  tmp2 * fjac(5,3,k+1)
     >              - tmp1 * njac(5,3,k+1)
               lhs(5,4,cc,k) =  tmp2 * fjac(5,4,k+1)
     >              - tmp1 * njac(5,4,k+1)
               lhs(5,5,cc,k) =  tmp2 * fjac(5,5,k+1)
     >              - tmp1 * njac(5,5,k+1)
     >              - tmp1 * dz5

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
c     c'(KMAX) and rhs'(KMAX) will be sent to next cell.
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     outer most do loops - sweeping in i direction
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     multiply c(i,j,0) by b_inverse and copy back to c
c     multiply rhs(0) by b_inverse(0) and copy to rhs
c---------------------------------------------------------------------
c            call binvcrhs( lhs(:,:,bb,0),
c     >                        lhs(:,:,cc,0),
c     >                        rhsz(:,i,j,0) )
           do kdx = 1, 5
             pivot = 1.00d0/lhs(kdx,kdx,bb,0)
             rhsz(kdx,i,j,0) = rhsz(kdx,i,j,0)*pivot
             do jdx = kdx+1, 5
               lhs(kdx,jdx,bb,0) = lhs(kdx,jdx,bb,0)*pivot
             end do
             do jdx = 1, 5
               lhs(kdx,jdx,cc,0) = lhs(kdx,jdx,cc,0)*pivot
             end do

             do jdx = 1, kdx-1
               coeff = lhs(jdx,kdx,bb,0)
               rhsz(jdx,i,j,0) = rhsz(jdx,i,j,0) - rhsz(kdx,i,j,0)*coeff
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
               rhsz(jdx,i,j,0) = rhsz(jdx,i,j,0) - rhsz(kdx,i,j,0)*coeff
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
      	    do k=1,ksize-1

c---------------------------------------------------------------------
c     subtract A*lhs_vector(k-1) from lhs_vector(k)
c
c     rhs(k) = rhs(k) - A*rhs(k-1)
c---------------------------------------------------------------------
c               call matvec_sub(lhs(:,:,aa,k),
c     >                         rhsz(:,i,j,k-1),rhsz(:,i,j,k))

c---------------------------------------------------------------------
c     B(k) = B(k) - C(k-1)*A(k)
c     call matmul_sub(aa,i,j,k,c,cc,i,j,k-1,c,bb,i,j,k)
c---------------------------------------------------------------------
c               call matmul_sub(lhs(:,:,aa,k),
c     >                         lhs(:,:,cc,k-1),
c     >                         lhs(:,:,bb,k))

           do kdx = 1, 5
             do jdx = 1, 5
               rhsz(kdx,i,j,k) = rhsz(kdx,i,j,k) -
     >                     lhs(kdx,jdx,aa,k)*rhsz(jdx,i,j,k-1)
               lhs(kdx,jdx,bb,k) = lhs(kdx,jdx,bb,k) -
     >               lhs(kdx,1,aa,k)*lhs(1,jdx,cc,k-1) -
     >               lhs(kdx,2,aa,k)*lhs(2,jdx,cc,k-1) -
     >               lhs(kdx,3,aa,k)*lhs(3,jdx,cc,k-1) -
     >               lhs(kdx,4,aa,k)*lhs(4,jdx,cc,k-1) -
     >               lhs(kdx,5,aa,k)*lhs(5,jdx,cc,k-1)
             end do
           end do
c---------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(i,j,1) by b_inverse(i,j,1) and copy to rhs
c---------------------------------------------------------------------
c               call binvcrhs( lhs(:,:,bb,k),
c     >                        lhs(:,:,cc,k),
c     >                        rhsz(:,i,j,k) )

           do kdx = 1, 5
             pivot = 1.00d0/lhs(kdx,kdx,bb,k)
             rhsz(kdx,i,j,k) = rhsz(kdx,i,j,k)*pivot
             do jdx = kdx+1, 5
               lhs(kdx,jdx,bb,k) = lhs(kdx,jdx,bb,k)*pivot
             end do
             do jdx = 1, 5
               lhs(kdx,jdx,cc,k) = lhs(kdx,jdx,cc,k)*pivot    !!!
             end do                                                     !!!


             do jdx = 1, kdx-1
             coeff = lhs(jdx,kdx,bb,k)
             rhsz(jdx,i,j,k) = rhsz(jdx,i,j,k) - rhsz(kdx,i,j,k)*coeff
               do idx = kdx+1, 5
                 lhs(jdx,idx,bb,k) = lhs(jdx,idx,bb,k) -
     >               lhs(kdx,idx,bb,k)*coeff
               end do
c               do idx = 1, 5
c                 lhs(jdx,idx,cc,k) = lhs(jdx,idx,cc,k) -
c     >               lhs(kdx,idx,cc,k)*coeff
c               end do
                 lhs(jdx,1,cc,k) = lhs(jdx,1,cc,k) -
     >               lhs(kdx,1,cc,k)*coeff
                 lhs(jdx,2,cc,k) = lhs(jdx,2,cc,k) -
     >               lhs(kdx,2,cc,k)*coeff
                 lhs(jdx,3,cc,k) = lhs(jdx,3,cc,k) -
     >               lhs(kdx,3,cc,k)*coeff
                 lhs(jdx,4,cc,k) = lhs(jdx,4,cc,k) -
     >               lhs(kdx,4,cc,k)*coeff
                 lhs(jdx,5,cc,k) = lhs(jdx,5,cc,k) -
     >               lhs(kdx,5,cc,k)*coeff
             end do

             do jdx = kdx+1,5
             coeff = lhs(jdx,kdx,bb,k)
             rhsz(jdx,i,j,k) = rhsz(jdx,i,j,k) - rhsz(kdx,i,j,k)*coeff
               do idx = kdx+1, 5
                 lhs(jdx,idx,bb,k) = lhs(jdx,idx,bb,k) -     !!!
     >               lhs(kdx,idx,bb,k)*coeff
               end do                                                  !!!
c               do idx = 1, 5
c                 lhs(jdx,idx,cc,k) = lhs(jdx,idx,cc,k) -
c     >               lhs(kdx,idx,cc,k)*coeff
c               end do
                 lhs(jdx,1,cc,k) = lhs(jdx,1,cc,k) -
     >               lhs(kdx,1,cc,k)*coeff
                 lhs(jdx,2,cc,k) = lhs(jdx,2,cc,k) -
     >               lhs(kdx,2,cc,k)*coeff
                 lhs(jdx,3,cc,k) = lhs(jdx,3,cc,k) -
     >               lhs(kdx,3,cc,k)*coeff
                 lhs(jdx,4,cc,k) = lhs(jdx,4,cc,k) -
     >               lhs(kdx,4,cc,k)*coeff
                 lhs(jdx,5,cc,k) = lhs(jdx,5,cc,k) -
     >               lhs(kdx,5,cc,k)*coeff
             end do
           end do

            enddo

c---------------------------------------------------------------------
c     Now finish up special cases for last cell
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
c---------------------------------------------------------------------
c            call matvec_sub(lhs(:,:,aa,ksize),
c     >                         rhsz(:,i,j,ksize-1),rhsz(:,i,j,ksize))

c---------------------------------------------------------------------
c     B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
c     call matmul_sub(aa,i,j,ksize,c,
c     $              cc,i,j,ksize-1,c,bb,i,j,ksize)
c---------------------------------------------------------------------
c            call matmul_sub(lhs(:,:,aa,ksize),
c     >                         lhs(:,:,cc,ksize-1),
c     >                         lhs(:,:,bb,ksize))

           do kdx = 1, 5
             do jdx = 1, 5
               rhsz(kdx,i,j,ksize) = rhsz(kdx,i,j,ksize) -
     >              lhs(kdx,jdx,aa,ksize)*rhsz(jdx,i,j,ksize-1)

               lhs(kdx,jdx,bb,ksize) = lhs(kdx,jdx,bb,ksize) -
     >            lhs(kdx,1,aa,ksize)*lhs(1,jdx,cc,ksize-1) -
     >            lhs(kdx,2,aa,ksize)*lhs(2,jdx,cc,ksize-1) -
     >            lhs(kdx,3,aa,ksize)*lhs(3,jdx,cc,ksize-1) -
     >            lhs(kdx,4,aa,ksize)*lhs(4,jdx,cc,ksize-1) -
     >            lhs(kdx,5,aa,ksize)*lhs(5,jdx,cc,ksize-1)
             end do
           end do
c---------------------------------------------------------------------
c     multiply rhs(ksize) by b_inverse(ksize) and copy to rhs
c---------------------------------------------------------------------
c            call binvrhs( lhs(:,:,bb,ksize),
c     >                       rhsz(:,i,j,ksize) )

          do kdx = 1, 5
             pivot = 1.00d0/lhs(kdx,kdx,bb,ksize)
             rhsz(kdx,i,j,ksize) = rhsz(kdx,i,j,ksize)*pivot
             do jdx = kdx+1, 5
               lhs(kdx,jdx,bb,ksize) = lhs(kdx,jdx,bb,ksize)*pivot
             end do

             do jdx = 1, kdx-1
               coeff = lhs(jdx,kdx,bb,ksize)
               rhsz(jdx,i,j,ksize) = rhsz(jdx,i,j,ksize) -
     >                               rhsz(kdx,i,j,ksize)*coeff
               do idx = kdx+1, 5
                 lhs(jdx,idx,bb,ksize) = lhs(jdx,idx,bb,ksize) -
     >                                   lhs(kdx,idx,bb,ksize)*coeff
               end do
             end do

             do jdx = kdx+1,5
               coeff = lhs(jdx,kdx,bb,ksize)
               rhsz(jdx,i,j,ksize) = rhsz(jdx,i,j,ksize) -
     >                               rhsz(kdx,i,j,ksize)*coeff
               do idx = kdx+1, 5
                 lhs(jdx,idx,bb,ksize) = lhs(jdx,idx,bb,ksize) -
     >                                   lhs(kdx,idx,bb,ksize)*coeff
               end do
             end do
           end do

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     back solve: if last cell, then generate U(ksize)=rhs(ksize)
c     else assume U(ksize) is loaded in un pack backsub_info
c     so just use it
c     after call u(kstart) will be sent to next cell
c---------------------------------------------------------------------

      	    do k=ksize-1,0,-1
               do m=1,BLOCK_SIZE
                  do n=1,BLOCK_SIZE
                     rhsz(m,i,j,k) = rhsz(m,i,j,k)
     >                    - lhs(m,n,cc,k)*rhsz(n,i,j,k+1)
                  enddo
               enddo
            enddo

         enddo
      enddo
      if (timeron) call timer_stop(t_zsolve)
      if (timeron) call timer_start(t_rdis2)
      rhs(:,:,0:JMAX-1,:)=rhsz(:,:,0:JMAX-1,:)
      if (timeron) call timer_stop(t_rdis2)

      return
      end
