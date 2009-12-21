c Subroutine that simulates ape smearing taking the product of links in the opposite 
c  direction to 'staples'
c
c Made into a module by F. Bissey Dec. 2003
c
      MODULE L_ape_smear

      CONTAINS

      subroutine ape_smear(ur,ui,alpha)

      USE GS_LATTICESIZE
      USE GS_NEWFIXSU3
      USE GS_SQUARES

c     global variables

      integer,parameter                                         :: nc=3 !sigma,color
      integer,parameter                                         :: nd=4 !direction


      double precision,dimension(nx,ny,nz,nt,nd,nc,nc)          :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

      double precision                                          :: alpha          ! smearing fraction

c  local variables

      integer                                                   :: xhat

      double precision,dimension(nx,ny,nz,nt,nc,nc)             ::stapler,staplei !smeared links
!HPF$ DISTRIBUTE stapler(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE staplei(*,*,BLOCK,BLOCK,*,*)

      double precision,dimension(nx,ny,nz,nt,nd,nc,nc)          :: ur_prm,ui_prm
!HPF$ DISTRIBUTE ur_prm(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui_prm(*,*,BLOCK,BLOCK,*,*,*)


      integer                                                   :: ic,jc          !counters
      double precision                                          :: g              ! g = alpha/(2*(nd-1))

c     starting the execution commands
      g = alpha/(2*(nd-1))

      do xhat=1,nd

         call squares(ur,ui,stapler,staplei,xhat,.true.)
 
         do ic=1,nc
           do jc=1,nc

               ur_prm(:,:,:,:,xhat,ic,jc) = (1-alpha)*ur(:,:,:,:,xhat,ic,jc) + g*stapler(:,:,:,:,jc,ic)
               ui_prm(:,:,:,:,xhat,ic,jc) = (1-alpha)*ui(:,:,:,:,xhat,ic,jc) - g*staplei(:,:,:,:,jc,ic)

            end do
         end do
      end do 

c     ur = ur_prm
c     ui = ui_prm

      call Newfixsu3(ur,ui,ur_prm,ui_prm,8)

      end subroutine ape_smear

      END MODULE L_ape_smear
