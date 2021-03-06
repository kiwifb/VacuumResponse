c
c     subroutine to calculate valules for Wilson loops, W, for 'VY' shape
!        adapted from 'Delta' shape
c     1.0 20060110 FB
!     2.0 FB071129 - ported to use boxes.f and hopefully correct a bug causing gauge
!        invariance problem in shape 7.
c
      MODULE L_LOOPS

      CONTAINS

c========================================================================================================

      subroutine delta_loops(xdir,ydir,ur,ui,W)
c
c     modules
c
      USE L_baryonParam
      USE L_epsilonIndex
      USE L_product
c
c     no implicit typing
c
      implicit none
c
c     global variables
c
#include"loopsize.f"
      integer                                                           :: xdir,ydir
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)                  :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
c
c     the five cases (1-nYloop) are stored in the second to last place in 'W'
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))          :: W
!HPF$ DISTRIBUTE W(*,*,BLOCK,BLOCK,*,*)
      integer                                                    :: time1, time2, count_rate
      integer                                                    :: timea, timeb
c
      include'loop_declarations.f'
c
c     ----------------
c     execution begins
c     ----------------
c
c     ---------------------------------
c     start of Wilson loop calculations
c     ---------------------------------
c
c     compute the time link loop
c
      do it = 1, nL(that)

         call product(ur,ui,that,tmp1r,tmp1i,it)
         tlr(:,:,:,:,:,:,it) = tmp1r(:,:,:,:,:,:)
         tli(:,:,:,:,:,:,it) = tmp1i(:,:,:,:,:,:)

      end do
c
c     CASE 1 - quarks at ( 1,0), (-1,-1) and (-1, 1)
c
      call system_clock(time1,count_rate)
      include'Cone.f'
      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'loops 1: ',(real(time2-time1)/real(count_rate)),' s'
c
c     CASE 2 - quarks at ( 2,0), (-1,-2) and (-1, 2)
c
      include'Ctwo.f'
      call system_clock(time1)
      write(*,'(a,f15.8,a)') 'loops 2: ',(real(time1-time2)/real(count_rate)),' s'
c
c     CASE 3 - quarks at ( 3,0), (-1,-2) and (-1, 2)
c
      include'Cthree.f'
      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'loops 3: ',(real(time2-time1)/real(count_rate)),' s'
c
c     CASE 4 - quarks at ( 3,0), (-2,-3) and (-2, 3)
c
      include'Cfour.f'
      call system_clock(time1)
      write(*,'(a,f15.8,a)') 'loops 4: ',(real(time1-time2)/real(count_rate)),' s'
c
c     CASE 5 - quarks at ( 4,0), (-3,-4) and (-3, 4)
c
      include'Cfive.f'
      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'loops 5: ',(real(time2-time1)/real(count_rate)),' s'
c
c     CASE 6 - quarks at ( 5,0), (-4,-5) and (-4, 5)
c
      include'Csix.f'
      call system_clock(time1)
      write(*,'(a,f15.8,a)') 'loops 6: ',(real(time1-time2)/real(count_rate)),' s'
c
c     CASE 7 - quarks at ( 7,0), (-4,-6) and (-4, 6)
c
      include'Cseven.f'
      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'loops 7: ',(real(time2-time1)/real(count_rate)),' s'

      return

      end subroutine delta_loops

      END MODULE L_LOOPS
