!
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine that fixes the su3 links. This subroutine needs to
!     called after a certain amount thermalisation has been done, the
!     purpose being to keep the links within the SU(3) algebra. This
!     subroutines forces the condition of unity U*Udag=I. In addition
!     it is called by ReadLinks and Write links to recover the full
!     SU(3) matrices.
!
!     Author: Richard Woloshyn
!             Parallelized by Frederi! D.R. Bonnet: September 1998
!             Supervisors: Derek Leinweber and Tony Williams
!
      MODULE GS_FIXSU3

      CONTAINS

      subroutine fixsu3(ur,ui)

      USE GS_LATTICESIZE

      implicit none

!     global variables

      integer,parameter                                       :: nc=3
      integer,parameter                                       :: mu=4

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

!     local variables

      integer,dimension(5)                                    :: yvector
!HPF$ DISTRIBUTE yvector(*)
      double precision,dimension(nx,ny,nz,nt)                 :: normr,normi
!HPF$ DISTRIBUTE normr(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE normi(*,*,BLOCK,BLOCK)
      integer                                                 :: jc,imu

!
!     First create an array to be looped over below
!     y(1) = 1 , y(2) = 2 , y(3) = 3 , y(4) = 1 , y(5) = 2

      do jc=1,nc
         yvector(jc) = jc
      end do
      do jc=4,5
         yvector(jc) = jc - 3
      end do

!     We'll do a loop here to save memory demands

      do imu=1,mu
!
!     first normalise first row
!
         normr = sqrt( ur(:,:,:,:,imu,1,1)**2 + ur(:,:,:,:,imu,1,2)**2 +
     &                 ui(:,:,:,:,imu,1,1)**2 + ui(:,:,:,:,imu,1,2)**2 +
     &                 ur(:,:,:,:,imu,1,3)**2 + ui(:,:,:,:,imu,1,3)**2 )

         do jc=1,nc
            ur(:,:,:,:,imu,1,jc) = ur(:,:,:,:,imu,1,jc) / normr(:,:,:,:)
            ui(:,:,:,:,imu,1,jc) = ui(:,:,:,:,imu,1,jc) / normr(:,:,:,:)
         end do
!
!     now compute row2 - (row2 dot row1)*row1
!
         normr = ur(:,:,:,:,imu,2,1) * ur(:,:,:,:,imu,1,1) +
     &           ui(:,:,:,:,imu,2,1) * ui(:,:,:,:,imu,1,1) +
     &           ur(:,:,:,:,imu,2,2) * ur(:,:,:,:,imu,1,2) +
     &           ui(:,:,:,:,imu,2,2) * ui(:,:,:,:,imu,1,2) +
     &           ur(:,:,:,:,imu,2,3) * ur(:,:,:,:,imu,1,3) +
     &           ui(:,:,:,:,imu,2,3) * ui(:,:,:,:,imu,1,3)

         normi = ui(:,:,:,:,imu,2,1) * ur(:,:,:,:,imu,1,1) -
     &           ur(:,:,:,:,imu,2,1) * ui(:,:,:,:,imu,1,1) +
     &           ui(:,:,:,:,imu,2,2) * ur(:,:,:,:,imu,1,2) -
     &           ur(:,:,:,:,imu,2,2) * ui(:,:,:,:,imu,1,2) +
     &           ui(:,:,:,:,imu,2,3) * ur(:,:,:,:,imu,1,3) -
     &           ur(:,:,:,:,imu,2,3) * ui(:,:,:,:,imu,1,3)

         do jc=1,nc
            ur(:,:,:,:,imu,2,jc) = ur(:,:,:,:,imu,2,jc) -
     &           (  normr * ur(:,:,:,:,imu,1,jc)
     &            - normi * ui(:,:,:,:,imu,1,jc) )
            ui(:,:,:,:,imu,2,jc) = ui(:,:,:,:,imu,2,jc) -
     &           (  normr * ui(:,:,:,:,imu,1,jc)
     &            + normi * ur(:,:,:,:,imu,1,jc) )
         end do
!
!     Now normalise the second row
!
         normr = sqrt( ur(:,:,:,:,imu,2,1)**2 + ui(:,:,:,:,imu,2,1)**2 +
     &                 ur(:,:,:,:,imu,2,2)**2 + ui(:,:,:,:,imu,2,2)**2 +
     &                 ur(:,:,:,:,imu,2,3)**2 + ui(:,:,:,:,imu,2,3)**2 )

         do jc=1,nc
            ur(:,:,:,:,imu,2,jc) = ur(:,:,:,:,imu,2,jc) / normr(:,:,:,:)
            ui(:,:,:,:,imu,2,jc) = ui(:,:,:,:,imu,2,jc) / normr(:,:,:,:)
         end do
!
!     now generate row3 = (row1 cross row2)^*
!

         do jc=1,nc
            ur(:,:,:,:,imu,3,jc) =
     &           ur(:,:,:,:,imu,1,yvector(jc+1)) * ur(:,:,:,:,imu,2,yvector(jc+2)) -
     &           ui(:,:,:,:,imu,1,yvector(jc+1)) * ui(:,:,:,:,imu,2,yvector(jc+2)) -
     &           ur(:,:,:,:,imu,1,yvector(jc+2)) * ur(:,:,:,:,imu,2,yvector(jc+1)) +
     &           ui(:,:,:,:,imu,1,yvector(jc+2)) * ui(:,:,:,:,imu,2,yvector(jc+1))
            ui(:,:,:,:,imu,3,jc) =
     &         - ur(:,:,:,:,imu,1,yvector(jc+1)) * ui(:,:,:,:,imu,2,yvector(jc+2)) -
     &           ui(:,:,:,:,imu,1,yvector(jc+1)) * ur(:,:,:,:,imu,2,yvector(jc+2)) +
     &           ur(:,:,:,:,imu,1,yvector(jc+2)) * ui(:,:,:,:,imu,2,yvector(jc+1)) +
     &           ui(:,:,:,:,imu,1,yvector(jc+2)) * ur(:,:,:,:,imu,2,yvector(jc+1))
         end do

      end do

      return

      end subroutine fixsu3

      END MODULE GS_FIXSU3
