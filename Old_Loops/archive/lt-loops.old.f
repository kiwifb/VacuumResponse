c     subroutine to calculate values for T and L shaped wilson loops, WT* and WL*, 
c     with "that" being the time direction
c     Adapted from wilsonloop.f and the original shapeTloop.f by F. Bissey Dec. 2003
c
      MODULE LTLOOPS

      Contains

      subroutine lt_loops(upxr,upxi,upyr,upyi,ur,ui,ixq,iyq,WLpp,WLpm,WLmp,WLmm,WTpp,WTpm,WTmp,WTmm)

        USE baryonParam
        USE epsilonIndex
        USE product

        implicit none
        
        double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui        !links
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
        double precision,dimension(nx,ny,nz,nt,nc,nc)           :: upxr,upxi,upyr,upyi
!HPF$ DISTRIBUTE upxr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upxi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upyr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upyi(*,*,BLOCK,BLOCK,*,*)
        integer                                                 :: ixq,iyq
        double precision,dimension(nx,ny,nz,nt,nL(that))        :: WLpp,WLpm,WLmp,WLmm ! L-shape wilson loops
!HPF$ DISTRIBUTE WLpp(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WLpm(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WLmp(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WLmm(*,*,BLOCK,BLOCK,*)
        double precision,dimension(nx,ny,nz,nt,nL(that))        :: WTpp,WTpm,WTmp,WTmm ! T-shape wilson loops
!HPF$ DISTRIBUTE WTpp(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WTpm(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WTmp(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WTmm(*,*,BLOCK,BLOCK,*)
c
c  Local variables
c
        integer                                                 :: it           !size of wilson loops
        integer                                                 :: ic,jc,kc     !colour counters
c
        double precision,dimension(nx,ny,nz,nt,nc,nc)           :: uptr,upti    !time link products
!HPF$ DISTRIBUTE uptr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upti(*,*,BLOCK,BLOCK,*,*)
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
c
c     Res holds the value for the wilson loop being calculated. This must be summed to gain
c     the average of the loop over the whole lattice
c
        double precision,dimension(nx,ny,nz,nt)                 :: Res
!HPF$ DISTRIBUTE Res(*,*,BLOCK,BLOCK)
c
c  Execution Begins
c
        do it = 1,nL(that)

           call product(ur,ui,that,uptr,upti,it)
c
c-----------------------------------------------------------------------------------------
c
c  create staple in the positive x direction
c
           tmpr = 0.0d0
           tmpi = 0.0d0

           ustr = cshift(uptr(:,:,:,:,:,:),dim=xhat,shift=ixq) 
           usti = cshift(upti(:,:,:,:,:,:),dim=xhat,shift=ixq) 
           
           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc

                    tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                   + upxr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) 
     &                   - upxi(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)

                    tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                   + upxr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) 
     &                   + upxi(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)

                 end do
              end do
           end do
c
c     go backwards from the positive x direction, forming staple
c
           sxpr = 0.0d0
           sxpi = 0.0d0

           usr  = cshift(upxr(:,:,:,:,:,:),dim=that,shift=it) 
           usi  = cshift(upxi(:,:,:,:,:,:),dim=that,shift=it) 
           
           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc

                    sxpr(:,:,:,:,ic,jc) = sxpr(:,:,:,:,ic,jc) 
     &                   + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc)
     &                   + tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc)

                    sxpi(:,:,:,:,ic,jc) = sxpi(:,:,:,:,ic,jc) 
     &                   - tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc)
     &                   + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc)

                 end do
              end do
           end do
c
c---------------------------------------------------------------------------
c
c  create staple in the negative x direction
c
           tmpr = 0.0d0
           tmpi = 0.0d0

           usr  = cshift(upxr(:,:,:,:,:,:),dim=xhat,shift=-ixq) 
           usi  = cshift(upxi(:,:,:,:,:,:),dim=xhat,shift=-ixq) 

           ustr = cshift(uptr(:,:,:,:,:,:),dim=xhat,shift=-ixq) 
           usti = cshift(upti(:,:,:,:,:,:),dim=xhat,shift=-ixq) 

           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc

                    tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                   + usr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc) 
     &                   + usi(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)

                    tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                   + usr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc) 
     &                   - usi(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)
        
                 end do
              end do
           end do
c
c     go forwards from the negative x direction, forming staple
c
           sxnr = 0.0d0
           sxni = 0.0d0

           usr  = cshift(cshift(upxr(:,:,:,:,:,:), dim=that, shift=it), dim=xhat, shift=-ixq)
           usi  = cshift(cshift(upxi(:,:,:,:,:,:), dim=that, shift=it), dim=xhat, shift=-ixq)

           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc
                 
                    sxnr(:,:,:,:,ic,jc) = sxnr(:,:,:,:,ic,jc) 
     &                   + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc)
     &                   - tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc)
                       
                    sxni(:,:,:,:,ic,jc) = sxni(:,:,:,:,ic,jc) 
     &                   + tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc)
     &                   + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc)
                   
                 end do
              end do
           end do
c
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
c
c  create staple in the positive y direction
c
           tmpr = 0.0d0
           tmpi = 0.0d0
           
           ustr = cshift(uptr(:,:,:,:,:,:),dim=yhat,shift=iyq) 
           usti = cshift(upti(:,:,:,:,:,:),dim=yhat,shift=iyq) 
           
           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc
                    
                    tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                   + upyr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) 
     &                   - upyi(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)
                       
                    tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                   + upyr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) 
     &                   + upyi(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)
                   
                 end do
              end do
           end do
c
c     go backwards from the positive y direction, forming staple
c
           sypr = 0.0d0
           sypi = 0.0d0
           
           usr  = cshift(upyr(:,:,:,:,:,:),dim=that,shift=it) 
           usi  = cshift(upyi(:,:,:,:,:,:),dim=that,shift=it) 
              
           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc
                    
                    sypr(:,:,:,:,ic,jc) = sypr(:,:,:,:,ic,jc) 
     &                   + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc)
     &                   + tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc)
                       
                    sypi(:,:,:,:,ic,jc) = sypi(:,:,:,:,ic,jc) 
     &                   - tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc)
     &                   + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc)
                  
                 end do
              end do
           end do
c
c---------------------------------------------------------------------------
c
c  create staple in the negative y direction
c
           tmpr = 0.0d0
           tmpi = 0.0d0
           
           usr  = cshift(upyr(:,:,:,:,:,:),dim=yhat,shift=-iyq) 
           usi  = cshift(upyi(:,:,:,:,:,:),dim=yhat,shift=-iyq) 
           
           ustr = cshift(uptr(:,:,:,:,:,:),dim=yhat,shift=-iyq) 
           usti = cshift(upti(:,:,:,:,:,:),dim=yhat,shift=-iyq) 
           
           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc

                    tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                   + usr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc) 
     &                   + usi(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)

                    tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                   + usr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc) 
     &                   - usi(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)

                 end do
              end do
           end do
c
c     go forwards from the negative y direction, forming staple
c
           synr = 0.0d0
           syni = 0.0d0
           
           usr  = cshift(cshift(upyr(:,:,:,:,:,:), dim=that, shift=it), dim=yhat, shift=-iyq)
           usi  = cshift(cshift(upyi(:,:,:,:,:,:), dim=that, shift=it), dim=yhat, shift=-iyq)
          
           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc
                    
                    synr(:,:,:,:,ic,jc) = synr(:,:,:,:,ic,jc) 
     &                   + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc)
     &                   - tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc)
                   
                    syni(:,:,:,:,ic,jc) = syni(:,:,:,:,ic,jc) 
     &                   + tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc)
     &                   + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc)
                
                 end do
              end do
           end do
c
c==========================================================================
c==========================================================================
c==========================================================================
c
c  Construct the baryon for quarks at +iyq, -iyq and +ixq : WTpp
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
           WTpp(:,:,:,:,it) = Res(:,:,:,:) / 6.0d0
c
c  Construct the baryon for quarks at +iyq, -iyq and -ixq : WTmm
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
           WTmm(:,:,:,:,it) = Res(:,:,:,:)  / 6.0d0
c
c  Construct the baryon for quarks at +ixq, -ixq and +iyq : WTmp
c
           Res = 0.0d0

           do ic = 1,36
           
              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sxpr(:,:,:,:,la(ic),lap(ic)) * sxnr(:,:,:,:,lb(ic),lbp(ic)) * sypr(:,:,:,:,lc(ic),lcp(ic))
     &             - sxpr(:,:,:,:,la(ic),lap(ic)) * sxni(:,:,:,:,lb(ic),lbp(ic)) * sypi(:,:,:,:,lc(ic),lcp(ic))
     &             - sxpi(:,:,:,:,la(ic),lap(ic)) * sxnr(:,:,:,:,lb(ic),lbp(ic)) * sypi(:,:,:,:,lc(ic),lcp(ic))
     &             - sxpi(:,:,:,:,la(ic),lap(ic)) * sxni(:,:,:,:,lb(ic),lbp(ic)) * sypr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WTmp(:,:,:,:,it) = Res(:,:,:,:) / 6.0d0
c
c  Construct the baryon for quarks at +ixq, -ixq and -iyq : WTpm
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sxpr(:,:,:,:,la(ic),lap(ic)) * sxnr(:,:,:,:,lb(ic),lbp(ic)) * synr(:,:,:,:,lc(ic),lcp(ic))
     &             - sxpr(:,:,:,:,la(ic),lap(ic)) * sxni(:,:,:,:,lb(ic),lbp(ic)) * syni(:,:,:,:,lc(ic),lcp(ic))
     &             - sxpi(:,:,:,:,la(ic),lap(ic)) * sxnr(:,:,:,:,lb(ic),lbp(ic)) * syni(:,:,:,:,lc(ic),lcp(ic))
     &             - sxpi(:,:,:,:,la(ic),lap(ic)) * sxni(:,:,:,:,lb(ic),lbp(ic)) * synr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WTpm(:,:,:,:,it) = Res(:,:,:,:)  / 6.0d0
c
c  Construct the baryon for quarks at +iyq, 0 and +ixq : WLpp
c
           Res = 0.0d0

           do ic = 1,36
           
              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sypr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             - sypr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WLpp(:,:,:,:,it) = Res(:,:,:,:) / 6.0d0
c
c  Construct the baryon for quarks at +iyq, 0 and -ixq : WLmp
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sypr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             - sypr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WLmp(:,:,:,:,it) = Res(:,:,:,:)  / 6.0d0
c
c  Construct the baryon for quarks at -iyq, 0 and ixq : WLpm
c
           Res = 0.0d0

           do ic = 1,36
           
              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               synr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             - synr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - syni(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - syni(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WLpm(:,:,:,:,it) = Res(:,:,:,:)  / 6.0d0
c
c  Construct the baryon for quarks at -iyq, 0 and -ixq : WLmm
c
           Res = 0.0d0

           do ic = 1,36
           
              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               synr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             - synr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - syni(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - syni(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WLmm(:,:,:,:,it) = Res(:,:,:,:)  / 6.0d0
c

        end do                  !end L and T loop

        return

      end subroutine lt_loops

      END MODULE LTLOOPS
