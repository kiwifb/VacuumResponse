!     subroutine to calculate values for Dual-diquark shaped wilson loops, WPL* and WOP*,
!     with "that" being the time direction
!     Adapted from wilsonloop.f and the original shapeTloop.f by F. Bissey Dec. 2003
!     This particular version initially limited itself to 2 kinds of T shapes.
!
!     It is also supplemented with the facilities to change the plane of computation:
!     xdir and ydir define the plane.
!
      MODULE L_LOOPS

      Contains

      subroutine Planar_loops(xdir,ydir,ur,ui,WPL)

      USE L_baryonParam
      USE L_epsilonIndex
      USE L_product

      implicit none

      include'loopsize.f'

      integer                                                 :: xdir,ydir
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui        !links
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt,1:ysep,1:dqsep,nL(that))  :: WPL ! DDQ-shape wilson loops
!HPF$ DISTRIBUTE WPL(*,*,BLOCK,BLOCK,*,*,*)
!
!  Local variables
!
      integer                                                 :: it           !size of wilson loops
      integer                                                 :: ic,jc,kc     !colour counters
      integer                                                 :: idq, iys
!     links products in xdir (diquark sep.) and ydir (q-q sep.) directions
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: upr,upi,xlr,xli
!HPF$ DISTRIBUTE upr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE xlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE xli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc,1:nL(that)) :: tlr,tli !time link products table
!HPF$ DISTRIBUTE tlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: llr,lli      !staples lower links
!HPF$ DISTRIBUTE llr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE lli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: slr,sli      !time shifted lower links
!HPF$ DISTRIBUTE slr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: sxlr,sxli    !time shifted x-links
!HPF$ DISTRIBUTE slr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: stlr,stli  !shifted time links
!HPF$ DISTRIBUTE stlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE stli(*,*,BLOCK,BLOCK,*,*)
!
!  the following temporary array hold products of links during the calculation of loops
!
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: tmpr,tmpi
!HPF$ DISTRIBUTE tmpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmpi(*,*,BLOCK,BLOCK,*,*)
!
!  Staples for the final baryon construction
!
!         1     3
!         |dqsep| ysep
!         |-----|   +
!         |     | ysep
!         2     4
!
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: str,sti ! intermediate staple result
!HPF$ DISTRIBUTE str(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sti(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: st1r,st1i,st2r,st2i
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: st3r,st3i,st4r,st4i
!HPF$ DISTRIBUTE st1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE st1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE st2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE st2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE st3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE st3i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE st4r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE st4i(*,*,BLOCK,BLOCK,*,*)
!
!     Res holds the value for the wilson loop being calculated. This must be summed to gain
!     the average of the loop over the whole lattice
!
      double precision,dimension(nx,ny,nz,nt)                 :: Res
!HPF$ DISTRIBUTE Res(*,*,BLOCK,BLOCK)
!
!  Execution Begins
!
!  preparing the time link table
!
      do it = 1, nL(that)

         call product(ur,ui,that,tmpr,tmpi,it)
         tlr(:,:,:,:,:,:,it) = tmpr(:,:,:,:,:,:)
         tli(:,:,:,:,:,:,it) = tmpi(:,:,:,:,:,:)

      end do ! it
!
!  looping over the diquark separation
!
      do idq = 1, dqsep

        call product(ur,ui,xdir,xlr,xli,idq)
!
!     looping over the quark separation parameter
!
        do iys = 1, ysep

          call product(ur,ui,ydir,upr,upi,iys)
!
!  time loop
!
          do it = 1, nL(that)
!
!       -------------staple 1-------------
!       prepare the links by shifting up to ll, tl to st and time shifting ll to sl
!
            llr = upr
            lli = upi

            stlr(:,:,:,:,:,:) = cshift(tlr(:,:,:,:,:,:,it), dim=ydir, shift= iys)
            stli(:,:,:,:,:,:) = cshift(tli(:,:,:,:,:,:,it), dim=ydir, shift= iys)

      include 'ustaple.f'
            st1r = str
            st1i = sti
!
!       -------------staple 3-------------
!       prepare the links by shifting up to ll, tl to st and time shifting ll to sl
!
            llr = cshift(upr, dim=xdir, shift= idq)
            lli = cshift(upi, dim=xdir, shift= idq)

            stlr = cshift(stlr, dim=xdir, shift= idq)
            stli = cshift(stli, dim=xdir, shift= idq)

      include 'ustaple.f'
            st3r = str
            st3i = sti
!
!       -------------staple 2-------------
!       prepare the links by shifting up to ll, tl to st and time shifting ll to sl
!
            llr = cshift(upr, dim=ydir, shift=-iys)
            lli = cshift(upi, dim=ydir, shift=-iys)

            stlr(:,:,:,:,:,:) = cshift(tlr(:,:,:,:,:,:,it), dim=ydir, shift=-iys)
            stli(:,:,:,:,:,:) = cshift(tli(:,:,:,:,:,:,it), dim=ydir, shift=-iys)

      include 'dstaple.f'
            st2r = str
            st2i = sti
!
!       -------------staple 4-------------
!       prepare the links by shifting up to ll, tl to st and time shifting ll to sl
!
            llr = cshift(llr, dim=xdir, shift= idq)
            lli = cshift(lli, dim=xdir, shift= idq)

            stlr = cshift(stlr, dim=xdir, shift= idq)
            stli = cshift(stli, dim=xdir, shift= idq)

      include 'dstaple.f'
            st4r = str
            st4i = sti
!
!     time shifting x-link
!
            sxlr = cshift(xlr, dim=that , shift= it)
            sxli = cshift(xli, dim=that , shift= it)

            str = 0.0d0
            sti = 0.0d0

            do ic = 1, nc
              do jc = 1, nc
                do kc = 1, 36

                  str(:,:,:,:,ic,jc) = str(:,:,:,:,ic,jc) + fper(kc) * (                                                        &
                    st1r(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*xlr(:,:,:,:,lc(kc),jc)*sxlr(:,:,:,:,lcp(kc),ic) + &
                    st1r(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*xli(:,:,:,:,lc(kc),jc)*sxli(:,:,:,:,lcp(kc),ic) + &
                    st1r(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*xlr(:,:,:,:,lc(kc),jc)*sxli(:,:,:,:,lcp(kc),ic) + &
                    st1i(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*xlr(:,:,:,:,lc(kc),jc)*sxli(:,:,:,:,lcp(kc),ic) - &
                    st1r(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*xli(:,:,:,:,lc(kc),jc)*sxlr(:,:,:,:,lcp(kc),ic) - &
                    st1i(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*xli(:,:,:,:,lc(kc),jc)*sxlr(:,:,:,:,lcp(kc),ic) - &
                    st1i(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*xlr(:,:,:,:,lc(kc),jc)*sxlr(:,:,:,:,lcp(kc),ic) - &
                    st1i(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*xli(:,:,:,:,lc(kc),jc)*sxli(:,:,:,:,lcp(kc),ic)   &
                 )

                 sti(:,:,:,:,ic,jc) = sti(:,:,:,:,ic,jc) + fper(kc) * (                                                         &
                    st1i(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*xlr(:,:,:,:,lc(kc),jc)*sxlr(:,:,:,:,lcp(kc),ic) + &
                    st1r(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*xlr(:,:,:,:,lc(kc),jc)*sxlr(:,:,:,:,lcp(kc),ic) + &
                    st1r(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*xli(:,:,:,:,lc(kc),jc)*sxlr(:,:,:,:,lcp(kc),ic) - &
                    st1r(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*xlr(:,:,:,:,lc(kc),jc)*sxli(:,:,:,:,lcp(kc),ic) - &
                    st1i(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*xli(:,:,:,:,lc(kc),jc)*sxlr(:,:,:,:,lcp(kc),ic) + &
                    st1i(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*xlr(:,:,:,:,lc(kc),jc)*sxli(:,:,:,:,lcp(kc),ic) + &
                    st1i(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*xli(:,:,:,:,lc(kc),jc)*sxli(:,:,:,:,lcp(kc),ic) + &
                    st1r(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*xli(:,:,:,:,lc(kc),jc)*sxli(:,:,:,:,lcp(kc),ic)   &
                 )

                end do
              end do ! jc
            end do ! ic

            do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (                                                   &
                str(:,:,:,:,la(ic),lap(ic)) * st3r(:,:,:,:,lbp(ic),lb(ic)) * st4r(:,:,:,:,lcp(ic),lc(ic)) - &
                str(:,:,:,:,la(ic),lap(ic)) * st3i(:,:,:,:,lbp(ic),lb(ic)) * st4i(:,:,:,:,lcp(ic),lc(ic)) + &
                sti(:,:,:,:,la(ic),lap(ic)) * st3r(:,:,:,:,lbp(ic),lb(ic)) * st4i(:,:,:,:,lcp(ic),lc(ic)) + &
                sti(:,:,:,:,la(ic),lap(ic)) * st3i(:,:,:,:,lbp(ic),lb(ic)) * st4r(:,:,:,:,lcp(ic),lc(ic))   &
              )

            end do ! ic

            WPL(:,:,:,:,idq,iys,it) = Res(:,:,:,:)

          end do ! it
        end do ! iys
      end do ! idq

      return

      end subroutine Planar_loops

!========================================================================================================

      END MODULE L_LOOPS
