c     subroutine to calculate values for T shaped wilson loops, WT*, 
c     with "that" being the time direction
c     Adapted from wilsonloop.f and the original shapeTloop.f by F. Bissey Dec. 2003
c     This particular version initially limited itself to 2 kinds of T shapes.
c
c     It is also supplemented with the facilities to change the plane of computation:
c     xdir and ydir define the plane. 
c
      MODULE L_DQLOOPS

      include'DiQuarkfiles/loopsize.h'

      Contains

      subroutine t_loops(xdir,ydir,ur,ui,WTp,WTm,WT4p,WT4m)

      USE L_baryonParam
      USE L_epsilonIndex
      USE L_product

      implicit none

      integer                                                 :: xdir,ydir
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui        !links
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
c     T-shapes are labelled p(plus) and m(minus):
c     p: ( ixq; iyq;-iyq)
c     m: (-ixq; iyq;-iyq)
      double precision,dimension(nx,ny,nz,nt,0:ldq,nL(that))  :: WTp,WTm ! T-shape wilson loops
!HPF$ DISTRIBUTE WTp(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE WTm(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,0:ldq,nL(that))  :: WT4p,WT4m ! T-shape wilson loops
!HPF$ DISTRIBUTE WT4p(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE WT4m(*,*,BLOCK,BLOCK,*,*)
c
c  Local variables
c
      integer                                                 :: it           !size of wilson loops
      integer                                                 :: ic,jc,kc     !colour counters
      integer                                                 :: ixq
c     links products in xdir and ydir directions
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: upxr,upxi,upyr,upyi
!HPF$ DISTRIBUTE upxr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upxi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upyr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upyi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,1:ldq,nc,nc)     :: upxtr,upxti
!HPF$ DISTRIBUTE upxtr(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE upxti(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: upy2r,upy2i
!HPF$ DISTRIBUTE upy2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upy2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: tlr,tli    !time link products
!HPF$ DISTRIBUTE tlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: usr,usi      !shifted links
!HPF$ DISTRIBUTE usr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE usi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: ustr,usti    !shifted links
!HPF$ DISTRIBUTE ustr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE usti(*,*,BLOCK,BLOCK,*,*)
c
c  the following temporary array hold products of links during the calculation of loops
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: tmpr,tmpi
!HPF$ DISTRIBUTE tmpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmpi(*,*,BLOCK,BLOCK,*,*)
c
c  Stapes for the final baryon construction
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: sxpr,sxpi
!HPF$ DISTRIBUTE sxpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sxpi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: sxnr,sxni
!HPF$ DISTRIBUTE sxnr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sxni(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: sypr,sypi
!HPF$ DISTRIBUTE sypr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sypi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: synr,syni
!HPF$ DISTRIBUTE synr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE syni(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: syp4r,syp4i
!HPF$ DISTRIBUTE syp4r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE syp4i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: syn4r,syn4i
!HPF$ DISTRIBUTE syn4r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE syn4i(*,*,BLOCK,BLOCK,*,*)
c
c     Res holds the value for the wilson loop being calculated. This must be summed to gain
c     the average of the loop over the whole lattice
c
      double precision,dimension(nx,ny,nz,nt)                 :: Res
!HPF$ DISTRIBUTE Res(*,*,BLOCK,BLOCK)
c
c  Execution Begins
c
      upyr = ur(:,:,:,:,ydir,:,:)
      upyi = ui(:,:,:,:,ydir,:,:)

      call product(ur,ui,ydir,upy2r,upy2i,1)
      call product(ur,ui,ydir,upy2r,upy2i,2)
c
c  Prepare the xdir link
c
      do ixq = 1, ldq

         call product(ur,ui,xdir,upxr,upxi,ixq)
         upxtr(:,:,:,:,ixq,:,:) = upxr(:,:,:,:,:,:)
         upxti(:,:,:,:,ixq,:,:) = upxi(:,:,:,:,:,:)

      end do
c
c  create staple 1 in the positive y direction
c
      do it = 1,nL(that)

         call product(ur,ui,that,tlr,tli,it)

         tmpr = 0.0d0
         tmpi = 0.0d0

         ustr = cshift(tlr,dim=ydir,shift= 1)
         usti = cshift(tli,dim=ydir,shift= 1)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                + upyr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) - upyi(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                + upyr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) + upyi(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)

               end do
            end do
         end do
c
c     go backwards from the positive y direction, forming staple
c
         sypr = 0.0d0
         sypi = 0.0d0

         usr  = cshift(upyr,dim=that,shift= it)
         usi  = cshift(upyi,dim=that,shift= it)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  sypr(:,:,:,:,ic,jc) = sypr(:,:,:,:,ic,jc)
     &                + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc)

                  sypi(:,:,:,:,ic,jc) = sypi(:,:,:,:,ic,jc)
     &                - tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc)

               end do
            end do
         end do
c
c  create staple 1 in the negative y direction
c
         tmpr = 0.0d0
         tmpi = 0.0d0

         usr  = cshift(upyr,dim=ydir,shift=-1)
         usi  = cshift(upyi,dim=ydir,shift=-1)

         ustr = cshift(tlr,dim=ydir,shift=-1)
         usti = cshift(tli,dim=ydir,shift=-1)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                + usr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc) + usi(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                + usr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc) - usi(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)

               end do
            end do
         end do
c
c     go forwards from the negative y direction, forming staple
c
         synr = 0.0d0
         syni = 0.0d0

         usr  = cshift(usr, dim=that,shift= it)
         usi  = cshift(usi, dim=that,shift= it)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  synr(:,:,:,:,ic,jc) = synr(:,:,:,:,ic,jc)
     &                + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc)

                  syni(:,:,:,:,ic,jc) = syni(:,:,:,:,ic,jc)
     &                + tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc)

               end do
            end do
         end do
c
c-----------------------------------------------------------
c
c  create staple 2 in the positive y direction
c
         tmpr = 0.0d0
         tmpi = 0.0d0

         ustr = cshift(tlr,dim=ydir,shift= 2)
         usti = cshift(tli,dim=ydir,shift= 2)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                + upy2r(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) - upy2i(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                + upy2r(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) + upy2i(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)

               end do
            end do
         end do
c
c     go backwards from the positive y direction, forming staple
c
         syp4r = 0.0d0
         syp4i = 0.0d0

         usr  = cshift(upy2r,dim=that,shift= it)
         usi  = cshift(upy2i,dim=that,shift= it)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  syp4r(:,:,:,:,ic,jc) = syp4r(:,:,:,:,ic,jc)
     &                + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc)

                  syp4i(:,:,:,:,ic,jc) = syp4i(:,:,:,:,ic,jc)
     &                - tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc)

               end do
            end do
         end do
c
c  create staple 2 in the negative y direction
c
         tmpr = 0.0d0
         tmpi = 0.0d0

         usr  = cshift(upy2r,dim=ydir,shift=-2)
         usi  = cshift(upy2i,dim=ydir,shift=-2)

         ustr = cshift(tlr,dim=ydir,shift=-2)
         usti = cshift(tli,dim=ydir,shift=-2)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                + usr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc) + usi(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                + usr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc) - usi(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)

               end do
            end do
         end do
c
c     go forwards from the negative y direction, forming staple
c
         syn4r = 0.0d0
         syn4i = 0.0d0

         usr  = cshift(usr, dim=that,shift= it)
         usi  = cshift(usi, dim=that,shift= it)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  syn4r(:,:,:,:,ic,jc) = syn4r(:,:,:,:,ic,jc)
     &                + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc)

                  syn4i(:,:,:,:,ic,jc) = syn4i(:,:,:,:,ic,jc)
     &                + tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc)

               end do
            end do
         end do
c
c--------------------------------------------------------------
c
         do ixq = 1, ldq

            upxr(:,:,:,:,:,:) = upxtr(:,:,:,:,ixq,:,:)
            upxi(:,:,:,:,:,:) = upxti(:,:,:,:,ixq,:,:)
c
c  create staple in the positive x direction
c
            tmpr = 0.0d0
            tmpi = 0.0d0

            ustr = cshift(tlr,dim=xdir,shift= ixq)
            usti = cshift(tli,dim=xdir,shift= ixq)

            do ic = 1,nc
               do jc = 1,nc
                  do kc = 1,nc

                     tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                   + upxr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) - upxi(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)

                     tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                   + upxr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) + upxi(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)

                  end do
               end do
            end do
c
c     go backwards from the positive x direction, forming staple
c
            sxpr = 0.0d0
            sxpi = 0.0d0

            usr  = cshift(upxr,dim=that,shift= it)
            usi  = cshift(upxi,dim=that,shift= it)

            do ic = 1,nc
               do jc = 1,nc
                  do kc = 1,nc

                     sxpr(:,:,:,:,ic,jc) = sxpr(:,:,:,:,ic,jc)
     &                   + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc)

                     sxpi(:,:,:,:,ic,jc) = sxpi(:,:,:,:,ic,jc)
     &                   - tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc)

                  end do
               end do
            end do
c
c  create staple in the negative x direction
c
            tmpr = 0.0d0
            tmpi = 0.0d0

            usr  = cshift(upxr,dim=xdir,shift=-ixq)
            usi  = cshift(upxi,dim=xdir,shift=-ixq)

            ustr = cshift(tlr,dim=xdir,shift=-ixq)
            usti = cshift(tli,dim=xdir,shift=-ixq)

            do ic = 1,nc
               do jc = 1,nc
                  do kc = 1,nc

                     tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                   + usr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)  + usi(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)

                     tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                   + usr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)  - usi(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)

                  end do
               end do
            end do
c
c     go forwards from the negative x direction, forming staple
c
            sxnr = 0.0d0
            sxni = 0.0d0

            usr  = cshift(usr,dim=that,shift= it)
            usi  = cshift(usi,dim=that,shift= it)

            do ic = 1,nc
               do jc = 1,nc
                  do kc = 1,nc

                     sxnr(:,:,:,:,ic,jc) = sxnr(:,:,:,:,ic,jc)
     &                   + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc)

                     sxni(:,:,:,:,ic,jc) = sxni(:,:,:,:,ic,jc)
     &                   + tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc)

                  end do
               end do
            end do
c
c  Construct the baryon for quarks at +iyq, -iyq and +ixq : WTp
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sypr(:,:,:,:,la(ic),lap(ic)) * synr(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             - sypr(:,:,:,:,la(ic),lap(ic)) * syni(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * synr(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * syni(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WTp(:,:,:,:,ixq,it) = Res(:,:,:,:) / 6.0d0
c
c  Construct the baryon for quarks at +iyq, -iyq and -ixq : WTm
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sypr(:,:,:,:,la(ic),lap(ic)) * synr(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             - sypr(:,:,:,:,la(ic),lap(ic)) * syni(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * synr(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * syni(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WTm(:,:,:,:,ixq,it) = Res(:,:,:,:)  / 6.0d0
c
c  Construct the baryon for quarks at +iyq, -iyq and +ixq : WT4p
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               syp4r(:,:,:,:,la(ic),lap(ic)) * syn4r(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             - syp4r(:,:,:,:,la(ic),lap(ic)) * syn4i(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - syp4i(:,:,:,:,la(ic),lap(ic)) * syn4r(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - syp4i(:,:,:,:,la(ic),lap(ic)) * syn4i(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WT4p(:,:,:,:,ixq,it) = Res(:,:,:,:) / 6.0d0
c
c  Construct the baryon for quarks at +iyq, -iyq and -ixq : WTm
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               syp4r(:,:,:,:,la(ic),lap(ic)) * syn4r(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             - syp4r(:,:,:,:,la(ic),lap(ic)) * syn4i(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - syp4i(:,:,:,:,la(ic),lap(ic)) * syn4r(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - syp4i(:,:,:,:,la(ic),lap(ic)) * syn4i(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WT4m(:,:,:,:,ixq,it) = Res(:,:,:,:)  / 6.0d0

        end do      ! end ixq loop

      end do   ! end it loop

      return

      end subroutine t_loops

c========================================================================================================

      END MODULE L_DQLOOPS
