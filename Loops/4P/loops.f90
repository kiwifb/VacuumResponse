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

      subroutine Tetra_loops(ur,ui,WOP)

      USE L_baryonParam
      USE L_epsilonIndex
      USE L_product
      USE L_BOXES

      implicit none

      include'loopsize.f'

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui        !links
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))  :: WOP ! DDQ-shape wilson loops
!HPF$ DISTRIBUTE WPL(*,*,BLOCK,BLOCK,*,*,*)
!
!  Local variables
!
      integer                                                 :: it           !size of wilson loops
      integer                                                 :: ic,jc,kc     !colour counters
      integer                                                 :: ilp, mlp
!     links products in xdir (diquark sep.) and ydir (q-q sep.) directions
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: u1r,u1i,u2r,u2i ! xy plane: 1\ 2/ unit links
!HPF$ DISTRIBUTE u1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE u1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE u2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE u2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: u3r,u3i,u4r,u4i ! xz plane: \4 /3 unit links
!HPF$ DISTRIBUTE u3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE u3i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE u4r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE u4i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: b1r,b1i,b2r,b2i ! xy plane: 1\ 2/ link products
!HPF$ DISTRIBUTE b1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE b1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE b2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE b2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: b3r,b3i,b4r,b4i ! xz plane: \4 /3 link products
!HPF$ DISTRIBUTE b3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE b3i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE b4r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE b4i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc,1:nL(that)) :: tlr,tli !time link products table
!HPF$ DISTRIBUTE tlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: llr,lli      !staples lower links
!HPF$ DISTRIBUTE llr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE lli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: slr,sli      !time shifted lower links
!HPF$ DISTRIBUTE slr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: stlr,stli  !shifted time links
!HPF$ DISTRIBUTE stlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE stli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: jur,jui  ! unit matrix for the junction
!HPF$ DISTRIBUTE jur(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE jui(*,*,BLOCK,BLOCK,*,*)
!
!  the following temporary array hold products of links during the calculation of loops
!
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: tmpr,tmpi
!HPF$ DISTRIBUTE tmpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmpi(*,*,BLOCK,BLOCK,*,*)
!
!  Staples for the final baryon construction
!
!         1
!         |         (1,4) quarks
!      2--|--3
!         |         (2,3) antiquarks
!         4
!
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: str,sti ! intermediate staple result
!HPF$ DISTRIBUTE str(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sti(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: st1r,st1i,st2r,st2i  ! "xy plane" staples
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: st3r,st3i,st4r,st4i  ! "xz plane" staples
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
! preparing the "junction" between the two epsilon contractions
      call product(ur,ui,that,jur,jui,0)
!
! preparing the unit links in the various directions
!
      call boxm1p1(xhat,yhat,ur,ui,u1r,u1i)
      call boxm1m1(xhat,yhat,ur,ui,u2r,u2i)
      call boxp1p1(xhat,zhat,ur,ui,u3r,u3i)
      call boxp1m1(xhat,zhat,ur,ui,u4r,u4i)
!
!  looping over the separation from the center of the tretraedron
!
      do ilp = 1, nloop

        mlp = ilp - 1

        if ( ilp==1) then

          b1r = u1r
          b1i = u1i
          b2r = u2r
          b2i = u2i
          b3r = u3r
          b3i = u3i
          b4r = u4r
          b4i = u4i

        else

          call multboxes(xhat,yhat,-mlp, mlp,b1r,b1i,u1r,u1i)
          call multboxes(xhat,yhat,-mlp,-mlp,b2r,b2i,u2r,u2i)
          call multboxes(xhat,zhat, mlp, mlp,b3r,b3i,u3r,u3i)
          call multboxes(xhat,zhat, mlp,-mlp,b4r,b4i,u4r,u4i)

        end if
!
!  time loop
!
        do it = 1, nL(that)
!
!       -------------staple 1-------------
!       prepare the links by shifting up to ll, tl to st and time shifting ll to sl
!
          llr = b1r
          lli = b1i

          stlr(:,:,:,:,:,:) = cshift(cshift(tlr(:,:,:,:,:,:,it), dim=xhat, shift=-ilp), dim=yhat, shift= ilp)
          stli(:,:,:,:,:,:) = cshift(cshift(tli(:,:,:,:,:,:,it), dim=xhat, shift=-ilp), dim=yhat, shift= ilp)

      include 'ustaple.f90'
          st1r = str
          st1i = sti
!
!     -------------staple 2-------------
!     prepare the links by shifting up to ll, tl to st and time shifting ll to sl
!
          llr = b2r
          lli = b2i

          stlr = cshift(stlr, dim=yhat, shift=-2*ilp)
          stli = cshift(stli, dim=yhat, shift=-2*ilp)

      include 'ustaple.f90'
          st2r = str
          st2i = sti
!
!     -------------staple 3-------------
!     prepare the links by shifting up to ll, tl to st and time shifting ll to sl
!
          llr = b3r
          lli = b3i

          stlr(:,:,:,:,:,:) = cshift(cshift(tlr(:,:,:,:,:,:,it), dim=xhat, shift= ilp), dim=zhat, shift= ilp)
          stli(:,:,:,:,:,:) = cshift(cshift(tli(:,:,:,:,:,:,it), dim=xhat, shift= ilp), dim=zhat, shift= ilp)

      include 'ustaple.f90'
          st3r = str
          st3i = sti
!
!     -------------staple 4-------------
!     prepare the links by shifting up to ll, tl to st and time shifting ll to sl
!
          llr  = b4r
          lli  = b4i

          stlr = cshift(stlr, dim=zhat, shift=-2*ilp)
          stli = cshift(stli, dim=zhat, shift=-2*ilp)

      include 'ustaple.f90'
          st4r = str
          st4i = sti

          str = 0.0d0
          sti = 0.0d0
!
!         we can simplify since jui=0
!
          do ic = 1, nc
            do jc = 1, nc
              do kc = 1, 36

                str(:,:,:,:,ic,jc) = str(:,:,:,:,ic,jc) + fper(kc) * (                                                       &
                  st1r(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*jur(:,:,:,:,lc(kc),jc)*jur(:,:,:,:,lcp(kc),ic) + &
                  st1r(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*jui(:,:,:,:,lc(kc),jc)*jui(:,:,:,:,lcp(kc),ic) + &
                  st1r(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*jur(:,:,:,:,lc(kc),jc)*jui(:,:,:,:,lcp(kc),ic) + &
                  st1i(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*jur(:,:,:,:,lc(kc),jc)*jui(:,:,:,:,lcp(kc),ic) - &
                  st1r(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*jui(:,:,:,:,lc(kc),jc)*jur(:,:,:,:,lcp(kc),ic) - &
                  st1i(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*jui(:,:,:,:,lc(kc),jc)*jur(:,:,:,:,lcp(kc),ic) - &
                  st1i(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*jur(:,:,:,:,lc(kc),jc)*jur(:,:,:,:,lcp(kc),ic)  &
                 - st1i(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*jui(:,:,:,:,lc(kc),jc)*jui(:,:,:,:,lcp(kc),ic)   &
               )

               sti(:,:,:,:,ic,jc) = sti(:,:,:,:,ic,jc) + fper(kc) * (                                                        &
                  st1i(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*jur(:,:,:,:,lc(kc),jc)*jur(:,:,:,:,lcp(kc),ic) + &
                  st1r(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*jur(:,:,:,:,lc(kc),jc)*jur(:,:,:,:,lcp(kc),ic)  &
                +  st1r(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*jui(:,:,:,:,lc(kc),jc)*jur(:,:,:,:,lcp(kc),ic) - &
                  st1r(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*jur(:,:,:,:,lc(kc),jc)*jui(:,:,:,:,lcp(kc),ic) - &
                  st1i(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*jui(:,:,:,:,lc(kc),jc)*jur(:,:,:,:,lcp(kc),ic) + &
                  st1i(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*jur(:,:,:,:,lc(kc),jc)*jui(:,:,:,:,lcp(kc),ic) + &
                  st1i(:,:,:,:,la(kc),lap(kc))*st2r(:,:,:,:,lb(kc),lbp(kc))*jui(:,:,:,:,lc(kc),jc)*jui(:,:,:,:,lcp(kc),ic) + &
                  st1r(:,:,:,:,la(kc),lap(kc))*st2i(:,:,:,:,lb(kc),lbp(kc))*jui(:,:,:,:,lc(kc),jc)*jui(:,:,:,:,lcp(kc),ic)   &
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

          WOP(:,:,:,:,ilp,it) = Res(:,:,:,:)

        end do ! it
      end do ! ilp

      return

      end subroutine Tetra_loops

!========================================================================================================

      END MODULE L_LOOPS
