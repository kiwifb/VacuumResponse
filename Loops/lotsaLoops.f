cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Sundance's staples routine, used to calculate the contribution 
c     from the five loops.
c     Made into a module by F. Bissey Dec. 2003
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      MODULE L_LotsaLoops

      CONTAINS

      subroutine LotsaLoops(ur,ui,stapler,staplei,xhat,quad,loops,uzero,yhat,fcall)

      USE GS_LATTICESIZE
      USE GS_UU

      implicit none

c     local variables

      integer,parameter                                       :: nc=3
      integer,parameter                                       :: mu=4
      
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: stapler,staplei
!HPF$ DISTRIBUTE stapler(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE staplei(*,*,BLOCK,BLOCK,*,*)

      integer                                                 :: xhat
      integer                                                 :: loops
      double precision                                        :: uzero
      logical                                                 :: fcall

      integer                                                 :: quad
c      integer                                                 :: callcounter

c     The variable quad defines which quadrants the loop will be taken around, relative to
c     the starting link, for the purposes of calculating the topological charge, and also
c     as a consequence, for the action, replacing the local == .true. conditions which are
c     still used in the subroutine squares.
c     quad = 1 (top right quadrant) equivalent to local = false (OverImpWilsonAction and TopQandReconAction) 
c     quad = 2 (top left quadrant)                              (TopQandReconAction)
c     quad = 3 (bottom left quadrant)                           (TopQandReconAction)
c     quad = 4 (bottom right quadrant)                          (TopQandReconAction)
c     quad = 5 (all loops) equivalent to local = true           (OverImpCoolSweep)

c
c  Local Variables
c
c  Temporaries: for shifted and multilink objects 
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: x1r,x1i
!HPF$ DISTRIBUTE x1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE x1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: x2r,x2i
!HPF$ DISTRIBUTE x2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE x2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: x3r,x3i
!HPF$ DISTRIBUTE x3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE x3i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: y1r,y1i
!HPF$ DISTRIBUTE y1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE y1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: y2r,y2i
!HPF$ DISTRIBUTE y2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE y2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: y3r,y3i
!HPF$ DISTRIBUTE y3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE y3i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: t1r,t1i
!HPF$ DISTRIBUTE t1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: t2r,t2i
!HPF$ DISTRIBUTE t2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: linkr,linki
!HPF$ DISTRIBUTE linkr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE linki(*,*,BLOCK,BLOCK,*,*)

      integer                                                 :: yhat

      double precision                                        :: C1, C2, C3, C4, C5

      logical, parameter                                      :: report = .false.

      integer                                                 :: ic, kc

      double precision, dimension(nx,ny,nz,nt)                :: testactionsq
!HPF$ DISTRIBUTE testactionsq(*,*,BLOCK,BLOCK)

c
c====================================================================================
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c        Here we define the values of the constants which 
c        create an improved action.  C5 must be determined first as
c        other constants depend on it
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (.not.fcall) then      ! We are calculating the action
         if (loops == 0) then
c
c testing 1, 2 ,3 ...
c
            C5 = 1.0d0 
            C1 = 1.0d0
            C2 = 1.0d0
            C3 = 1.0d0
            C4 = 1.0d0        
c
         else if (loops == 1) then
c
c square
c
            C5 = 0.0d0
            C1 = 1.0d0
            C2 = 0.0d0
            C3 = 0.0d0
            C4 = 0.0d0
c
         else if (loops == 2) then
c
c square + rect
c
            C5 =  0.0d0
            C1 =  5.0d0/3.0d0
            C2 =  0.0d0
            C3 = -1.0d0/(12.0d0 * uzero**2)
            C4 =  0.0d0
c
         else if (loops == 3) then
c
c square + window + picture window
c
            C5 =  0.1d0
            C1 = 13.5d0/9.0d0
            C2 = -5.4d0/(144.0d0 * uzero**4)
            C3 =  0.0d0
            C4 =  0.0d0
c
         else if (loops == 4) then
c
c square + rect + window + ladder
c
            C5 =  0.0d0
            C1 = 19.0d0/9.0d0
            C2 =  1.0d0/(144.0d0 * uzero**4)
            C3 =-64.0d0/(360.0d0 * uzero**2)
            C4 =  1.0d0/( 90.0d0 * uzero**4)
c
         else if (loops == 5) then
c
c square + rect + window + ladder + picture window
c
            C5 =   0.05d0 
            C1 = (19.0d0 - C5* 55.0d0)/9.0d0
            C2 = ( 1.0d0 - C5* 64.0d0)/(144.0d0 * uzero**4)
            C3 =-(64.0d0 - C5*640.0d0)/(360.0d0 * uzero**2)
            C4 = ( 0.2d0 - C5*  2.0d0)/( 18.0d0 * uzero**4)
c
         else if (loops == 6) then
c
c ladder + window ( = elope action)
c
            C5 =  0.0d0
            C1 =  0.0d0
            C2 =  5.0d0/(16.0d0 * uzero**4)
            C3 =  0.0d0
            C4 = -4.0d0/( 9.0d0 * uzero**4)
c
c         else
c
c            pause
c
         end if
c
         C5 = C5/(81.0d0 * uzero**8)

      end if

      if (fcall) then                     ! we are calculating the topological charge
         if (loops == 0) then
c
c testing 1, 2 ,3 ...
c
            C5 = 1.0d0 
            C1 = 1.0d0
            C2 = 1.0d0
            C3 = 1.0d0
            C4 = 1.0d0        
c
         else if (loops == 1) then
c
c square
c
            C5 = 0.0d0
            C1 = 1.0d0
            C2 = 0.0d0
            C3 = 0.0d0
            C4 = 0.0d0
c
         else if (loops == 2) then
c
c square + rect
c
            C5 =  0.0d0
            C1 =  5.0d0/3.0d0
            C2 =  0.0d0
            C3 = -1.0d0/(6.0d0 * uzero**2)
            C4 =  0.0d0
c
         else if (loops == 3) then
c
c square + window + picture window
c
            C5 =  1.0d0/90.0d0 
            C1 = 13.5d0/9.0d0
            C2 = -5.4d0/(36.0d0 * uzero**4)
            C3 =  0.0d0
            C4 =  0.0d0
c
         else if (loops == 4) then
c
c square + rect + window + ladder ( = 4LTQ #1)
c
            C5 =    0.0d0
            C1 =   19.0d0/9.0d0
            C2 =    1.0d0/(36.0d0 * uzero**4)
            C3 = - 16.0d0/(45.0d0 * uzero**2)
            C4 =    1.0d0/( 30.0d0 * uzero**4)
c
         else if (loops == 5) then
c
c square + rect + window + ladder + picture window
c
            C5 = 1.0d0/180.0d0 
            C1 =       (19.0d0 - C5* 495.0d0)/9.0d0
            C2 =       ( 1.0d0 - C5* 576.0d0)/(36.0d0 * uzero**4)
            C3 = 16*(C5*90.0d0 - 1.0d0)/(45.0d0 * uzero**2)
            C4 =       ( 1.0d0 - C5* 90.0d0)/( 30.0d0 * uzero**4)
c
         else if (loops == 6) then
c
c rect + ladder + window + picture window ( = 4LTQ #2)
c
            C5 =  19.0d0/495.0d0
            C1 =       (19.0d0 - (19.0d0/495.0d0)* 495.0d0)/9.0d0
            C2 =       ( 1.0d0 - (19.0d0/495.0d0)* 576.0d0)/(36.0d0 * uzero**4)
            C3 = 16*((19.0d0/495.0d0)*90.0d0 - 1.0d0)/(45.0d0 * uzero**2)
            C4 =       ( 1.0d0 - (19.0d0/495.0d0)* 90.0d0)/( 30.0d0 * uzero**4)
c
         else if (loops == 7) then
 
c square + rect + window + picture window ( = 4LTQ #3)
c
            C5 =  1.0d0/576.0d0
            C1 =       (19.0d0 - (1.0d0/576.0d0)* 495.0d0)/9.0d0
            C2 =       ( 1.0d0 - (1.0d0/576.0d0)* 576.0d0)/(36.0d0 * uzero**4)
            C3 = 16*((1.0d0/576.0d0)*90.0d0 - 1.0d0)/(45.0d0 * uzero**2)
            C4 =       ( 1.0d0 - (1.0d0/576.0d0)* 90.0d0)/( 30.0d0 * uzero**4)

c
c         else
c
c            pause
c
         end if
c
         C5 = C5/(uzero**8)

      end if

c      callcounter = 0
c
c====================================================================================
c
c      stapler = 0.0d0
c      staplei = 0.0d0
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      if (loops /= 6) then
c     
c     Do the positive half of the (1x1) staple
c                                                        _
c                                                       | |
c      
         t1r = 0.0d0
         t1i = 0.0d0
         t2r = 0.0d0
         t2i = 0.0d0
   
         if (quad == 1 .or. quad == 2 .or. quad == 5) then

            y1r = cshift(ur(:,:,:,:,yhat,:,:),dim=xhat,shift=1)
            y1i = cshift(ui(:,:,:,:,yhat,:,:),dim=xhat,shift=1)
            x1r = cshift(ur(:,:,:,:,xhat,:,:),dim=yhat,shift=1)
            x1i = cshift(ui(:,:,:,:,xhat,:,:),dim=yhat,shift=1)

            linkr = ur(:,:,:,:,yhat,:,:)
            linki = ui(:,:,:,:,yhat,:,:)

            call UUdag(t1r,t1i,y1r,y1i,x1r,x1i)
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c       In the following multiplication we introduce a factor
c       C1 which will serve to create an improved action  
c
            t1r = C1*t1r
            t1i = C1*t1i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               
            call UUdag(stapler,staplei,t1r,t1i,linkr,linki)

            if (report) then
            
               testactionsq = 0.0d0
               
               do ic=1,nc
                     do kc=1,nc
               
                        testactionsq(:,:,:,:) = testactionsq(:,:,:,:)  
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic) 
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )
                    
                     end do
                  end do

                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for positive 1 x 1 loop'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C1 =', C1
c                  write(7,*)

               end if

            end if              !closes the quad = 1,2,5 condition   
c
c-------------------------------------------------------------------------------------
c
c     Do negative half of the (1x1) staple
c     
c                                                        |_|
 
            if (quad == 3 .or. quad == 4 .or. quad ==5 ) then
            
               y1r = cshift(cshift(ur(:,:,:,:,yhat,:,:),dim=xhat,shift = 1),
     &                      dim=yhat, shift = -1)
               y1i = cshift(cshift(ui(:,:,:,:,yhat,:,:),dim=xhat,shift = 1),
     &                      dim=yhat, shift = -1)
               x1r = cshift(ur(:,:,:,:,xhat,:,:),dim=yhat,shift = -1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:),dim=yhat,shift = -1)

               t1r = 0.0d0
               t1i = 0.0d0

               call UdagUdag(t1r,t1i,y1r,y1i,x1r,x1i)
               
c       Reassign a value to y1, to give the third link
                     
               y1r = cshift(ur(:,:,:,:,yhat,:,:), dim=yhat, shift = -1) 
               y1i = cshift(ui(:,:,:,:,yhat,:,:), dim=yhat, shift = -1)
                     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c       In the following multiplication we introduce a factor
c       C1 which will serve to create an improved action  
c
               t1r = C1*t1r
               t1i = C1*t1i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UU(stapler,staplei,t1r,t1i,y1r,y1i)
      
               testactionsq = 0.0d0
               
               do ic=1,nc
                  do kc=1,nc
               
                     testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)    
     &                    + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic) 
     &                    -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )
           
                  end do
               end do

c                callcounter = callcounter + 1
c                write(*,*)'That makes', callcounter, 'calls to UU'

               if (report) then

                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for negative 1 x 1 loop, i, yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C1 =', C1
c                  write(7,*)

               end if
c
            end if                      ! closes the quad == 3,4,5 condition
         end if                         ! closes the square
c
c====================================================================================
c
c        Constructing the 2 link products to be used through the subroutine
c
         if (loops /= 1) then
c
            x1r = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=1)
            x1i = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=1) 
   
            y1r = cshift(ur(:,:,:,:,yhat,:,:),dim=yhat,shift=1) 
            y1i = cshift(ui(:,:,:,:,yhat,:,:),dim=yhat,shift=1)
   
            x2r = 0.0d0
            x2i = 0.0d0
            y2r = 0.0d0
            y2i = 0.0d0
   
            linkr(:,:,:,:,:,:) = ur(:,:,:,:,xhat,:,:)
            linki(:,:,:,:,:,:) = ui(:,:,:,:,xhat,:,:)

            call UU(x2r,x2i,linkr,linki,x1r,x1i)

            linkr(:,:,:,:,:,:) = ur(:,:,:,:,yhat,:,:)
            linki(:,:,:,:,:,:) = ui(:,:,:,:,yhat,:,:)

            call UU(y2r,y2i,linkr,linki,y1r,y1i)

         end if
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c        Do the positive half of the (2x2) staple, update link at left
c
c                                                        __
c                                                       |  |
c                                                       | _|
c
         if (loops /= 1 .and. loops /= 2) then

            t1r = 0.0d0
            t1i = 0.0d0
            t2r = 0.0d0
            t2i = 0.0d0

c  Now we reuse x1r and x1i and use y1 as a temporary variable

            if (quad == 1 .or. quad == 5) then

               y1r = cshift(y2r(:,:,:,:,:,:),dim=xhat,shift=2)
               y1i = cshift(y2i(:,:,:,:,:,:),dim=xhat,shift=2) 
         
               call UU(t1r,t1i,x1r,x1i,y1r,y1i)

c  Now we redefine x1 to use it as a temporary variable

               x1r = cshift(x2r(:,:,:,:,:,:),dim=yhat,shift=2)
               x1i = cshift(x2i(:,:,:,:,:,:),dim=yhat,shift=2)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        In the following multiplication we introduce a factor
c        C2 which will serve to create an improved action
c
               t1r = C2*t1r
               t1i = C2*t1i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t2r,t2i,t1r,t1i,x1r,x1i)

               call UUdag(stapler,staplei,t2r,t2i,y2r,y2i)

               if (report) then
             
                  testactionsq = 0.0d0
     
                  do ic=1,nc
                     do kc=1,nc
       
                        testactionsq(:,:,:,:) = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )
         
                     end do
                  end do

                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)
c                  write(7,*)'the testactionsq for positive 2 x 2L loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C2 =', C2 
c                  write(7,*)

               end if
            end if                             !closes the quad == 1,5 condition 
c
c--------------------------------------------------------------------------------------------
c
c       Do the negative half of the (2x2) staple, update link at left
c
c                                                         _
c                                                       |  |
c                                                       |__|
c                          
            if (quad == 4 .or. quad == 5) then

               x1r = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=1)

               y1r = cshift(cshift(y2r(:,:,:,:,:,:),dim=xhat,shift=2),
     &                      dim= yhat, shift=-2)  
               y1i = cshift(cshift(y2i(:,:,:,:,:,:),dim=xhat,shift=2),
     &                      dim= yhat, shift=-2)
           
               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0
             
               call UUdag(t1r,t1i,x1r,x1i,y1r,y1i)

               x1r = cshift(x2r(:,:,:,:,:,:),dim=yhat,shift=-2)
               x1i = cshift(x2i(:,:,:,:,:,:),dim=yhat,shift=-2)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       In the following multiplication we introduce a factor
c       C2 which will serve to create an improved action
c
               t1r = C2*t1r
               t1i = C2*t1i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t2r,t2i,t1r,t1i,x1r,x1i)                  
 
               y1r = cshift(y2r(:,:,:,:,:,:),dim=yhat,shift=-2)
               y1i = cshift(y2i(:,:,:,:,:,:),dim=yhat,shift=-2)

               call UU(stapler,staplei,t2r,t2i,y1r,y1i)
        
               if (report) then
    
                  testactionsq = 0.0d0
    
                  do ic=1,nc
                     do kc=1,nc
      
                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )
    
                     end do
                  end do   
     
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)
c                  write(7,*)'the testactionsq for negative 2 x 2L loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C2 =', C2
c                  write(7,*)

               end if
            end if			! closes the quad == 4,5 condition   
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
c
c       Do the positive half of the (2x2) staple, update link at right
c                                                        __
c                                                       |  |
c                                                       |_ |
            if (quad == 2 .or. quad == 5 ) then 

               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0
     
               y1r = cshift(y2r(:,:,:,:,:,:),dim=xhat,shift=1)
               y1i = cshift(y2i(:,:,:,:,:,:),dim=xhat,shift=1)
           
               x1r = cshift(cshift(x2r(:,:,:,:,:,:),dim=yhat,shift=2),
     &                      dim = xhat, shift = -1)
               x1i = cshift(cshift(x2i(:,:,:,:,:,:),dim=yhat,shift=2),
     &                      dim = xhat, shift = -1)

               call UUdag(t1r,t1i,y1r,y1i,x1r,x1i)
 
               y1r = cshift(y2r(:,:,:,:,:,:),dim=xhat,shift=-1)
               y1i = cshift(y2i(:,:,:,:,:,:),dim=xhat,shift=-1)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       In the following multiplication we introduce a factor
c       C2 which will serve to create an improved action
c
               t1r = C2*t1r
               t1i = C2*t1i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t2r,t2i,t1r,t1i,y1r,y1i)

               x1r = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)

               call UU(stapler,staplei,t2r,t2i,x1r,x1i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:) = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic) 
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )
   
                     end do   
                  end do    

                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)
c                  write(7,*)'the testactionsq for positive 2 x 2R loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C2 =', C2
c                  write(7,*)

               end if
            end if                      !closes the quad == 2,5 condition         
c
c----------------------------------------------------------------------------------------------
c
c       Do the negative half of the (2x2) staple, update link at right
c                                                        _ 
c                                                       |  |
c                                                       |__|
            if (quad == 3 .or. quad == 5 ) then

               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0
                        
               y1r = cshift(cshift(y2r(:,:,:,:,:,:),dim=xhat,shift=1),
     &                      dim=yhat, shift=-2)
               y1i = cshift(cshift(y2i(:,:,:,:,:,:),dim=xhat,shift=1),
     &                      dim=yhat, shift=-2)               
         
               x1r = cshift(cshift(x2r(:,:,:,:,:,:),dim=yhat,shift=-2),
     &                      dim = xhat, shift = -1)
               x1i = cshift(cshift(x2i(:,:,:,:,:,:),dim=yhat,shift=-2),
     &                      dim = xhat, shift = -1)

               call UdagUdag(t1r,t1i,y1r,y1i,x1r,x1i)
         
               y1r = cshift(cshift(y2r(:,:,:,:,:,:),dim=xhat,shift=-1),
     &                      dim=yhat, shift=-2)
               y1i = cshift(cshift(y2i(:,:,:,:,:,:),dim=xhat,shift=-1),
     &                      dim=yhat, shift=-2)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       In the following multiplication we introduce a factor
c       C2 which will serve to create an improved action
c
               t1r = C2*t1r
               t1i = C2*t1i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UU(t2r,t2i,t1r,t1i,y1r,y1i)
            
               x1r = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)
                        
               call UU(stapler,staplei,t2r,t2i,x1r,x1i)

               if (report) then
  
                  testactionsq = 0.0d0
   
                  do ic=1,nc
                     do kc=1,nc
       
                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )
    
                     end do
                  end do
            
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)
c                  write(7,*)'the testactionsq for negative 2 x 2R loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C2 =', C2
c                  write(7,*)
 
               end if

            end if                ! closes the quad == 3,5 condition
         end if                   ! closes the if (loops /= 1 .and. loops /= 2)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
c
c     Do positive half of the (1x2) staple
c
c                                                        _
c                                                       | |
c                                                       | |
c
         if (loops == 0 .or. loops == 2 .or. loops == 4 .or. loops == 5) then

            if (quad == 1 .or. quad == 2 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0
   
               x1r = cshift(ur(:,:,:,:,xhat,:,:),dim=yhat,shift=2)
               x1i = cshift(ui(:,:,:,:,xhat,:,:),dim=yhat,shift=2)
    
               y1r = cshift(y2r(:,:,:,:,:,:),dim=xhat,shift=1)
               y1i = cshift(y2i(:,:,:,:,:,:),dim=xhat,shift=1)

               call UUdag(t1r,t1i,y1r,y1i,x1r,x1i)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       In the following multiplication we introduce a factor
c       C3 which will serve to create an improved action
c
               t1r = C3*t1r
               t1i = C3*t1i        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(stapler,staplei,t1r,t1i,y2r,y2i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:) = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )
     
                     end do
                  end do
      
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)
c                  write(7,*)'the testactionsq for positive 1 x 2 loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C3 =', C3
c                  write(7,*)

               end if
            end if                            ! closes the quad == 1,2,5 condition   
c
c-------------------------------------------------------------------------------------------
c
c     Do negative half of the (1x2) staple
c
c                                                       | |
c                                                       |_|

            if (quad == 3 .or. quad == 4 .or. quad == 5) then
    
               t1r = 0.0d0
               t1i = 0.0d0
 
               x1r = cshift(ur(:,:,:,:,xhat,:,:),dim=yhat,shift=-2)
               x1i = cshift(ui(:,:,:,:,xhat,:,:),dim=yhat,shift=-2)
                        
               y1r = cshift(cshift(y2r(:,:,:,:,:,:),dim=xhat,shift=1),
     &                      dim=yhat, shift=-2)
               y1i = cshift(cshift(y2i(:,:,:,:,:,:),dim=xhat,shift=1),
     &                      dim=yhat, shift=-2)

               call UdagUdag(t1r,t1i,y1r,y1i,x1r,x1i)

               y1r = cshift(y2r(:,:,:,:,:,:),dim=yhat, shift=-2)
               y1i = cshift(y2i(:,:,:,:,:,:),dim=yhat, shift=-2)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       In the following multiplication we introduce a factor
c       C3 which will serve to create an improved action
c
               t1r = C3*t1r
               t1i = C3*t1i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UU(stapler,staplei,t1r,t1i,y1r,y1i)

               if (report) then
 
                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:) = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )
 
                     end do
                  end do
               
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)
c                  write(7,*)'the testactionsq for negative 1 x 2 loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C3 =', C3
c                  write(7,*)
                  
               end if
            end if                ! closes the quad == 3,4,5 condition

c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Do positive half of the (2x1) staple, update link at left
c                                                        __
c                                                       | _|
            if (quad == 1 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0

               y1r = cshift(ur(:,:,:,:,yhat,:,:),dim=xhat,shift=2)
               y1i = cshift(ui(:,:,:,:,yhat,:,:),dim=xhat,shift=2)

               x1r = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=1)

               linkr(:,:,:,:,:,:) = ur(:,:,:,:,yhat,:,:)
               linki(:,:,:,:,:,:) = ui(:,:,:,:,yhat,:,:)

               call UU(t1r,t1i,x1r,x1i,y1r,y1i)       

               x1r = cshift(x2r(:,:,:,:,:,:),dim=yhat,shift =1)
               x1i = cshift(x2i(:,:,:,:,:,:),dim=yhat,shift =1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       In the following multiplication we introduce a factor
c       C3 which will serve to create an improved action
c
               t1r = C3*t1r
               t1i = C3*t1i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t2r,t2i,t1r,t1i,x1r,x1i)

               call UUdag(stapler,staplei,t2r,t2i,linkr,linki)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do   
     
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)
c                  write(7,*)'the testactionsq for positive 2 x 1L loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C3 =', C3
c                  write(7,*)

               end if
            end if                      ! closes the quad == 1,5 condition

c -----------------------------------------------------------------
c     Do negative half of the (2x1) staple, update link at left
c                                                         _
c                                                       |__|
            if (quad == 4 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0  
               t2r = 0.0d0
               t2i = 0.0d0
                                        
               y1r = cshift(cshift(ur(:,:,:,:,yhat,:,:),dim=xhat,shift=2),
     &                        dim=yhat, shift=-1)
               y1i = cshift(cshift(ui(:,:,:,:,yhat,:,:),dim=xhat,shift=2),
     &                        dim=yhat, shift=-1)
                        
               x1r = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=1)

               call UUdag(t1r,t1i,x1r,x1i,y1r,y1i)          

               x1r = cshift(x2r(:,:,:,:,:,:),dim=yhat,shift =-1)
               x1i = cshift(x2i(:,:,:,:,:,:),dim=yhat,shift =-1)  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
c       In the following multiplication we introduce a factor
c       C3 which will serve to create an improved action
c
               t1r = C3*t1r
               t1i = C3*t1i
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t2r,t2i,t1r,t1i,x1r,x1i)

               y1r = cshift(ur(:,:,:,:,yhat,:,:),dim=yhat,shift=-1)
               y1i = cshift(ui(:,:,:,:,yhat,:,:),dim=yhat,shift=-1)

               call UU(stapler,staplei,t2r,t2i,y1r,y1i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
        
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)
c                  write(7,*)'the testactionsq for negative 2 x 1L loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C3 =', C3
c                  write(7,*)

               end if
            end if                    ! closes the quad == 4,5 condition   
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Do positive half of the (2x1) staple, update link at right
c                                                        __
c                                                       |_ |

            if (quad == 2 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0
      
               y1r = cshift(ur(:,:,:,:,yhat,:,:),dim=xhat,shift=1)
               y1i = cshift(ui(:,:,:,:,yhat,:,:),dim=xhat,shift=1)
                
               x1r = cshift(cshift(x2r(:,:,:,:,:,:), dim=xhat, shift=-1),
     &                            dim = yhat, shift=1)
               x1i = cshift(cshift(x2i(:,:,:,:,:,:), dim=xhat, shift=-1),
     &                            dim = yhat, shift=1)

               call UUdag(t1r,t1i,y1r,y1i,x1r,x1i)

               y1r = cshift(ur(:,:,:,:,yhat,:,:),dim=xhat,shift=-1)
               y1i = cshift(ui(:,:,:,:,yhat,:,:),dim=xhat,shift=-1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
c       In the following multiplication we introduce a factor
c       C3 which will serve to create an improved action
c
               t1r = C3*t1r
               t1i = C3*t1i
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t2r,t2i,t1r,t1i,y1r,y1i)

               x1r = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)

               call UU(stapler,staplei,t2r,t2i,x1r,x1i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
        
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)
c                  write(7,*)'the testactionsq for positive 2 x 1R loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C3 =', C3
c                  write(7,*)

               end if
            end if                         ! closes the quad == 4,5 condition
c   
c ---------------------------------------------------------------------
c     Do negative half of the (2x1) staple, update link at right
c                                                        _ 
c                                                       |__|
            
            if (quad == 3 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0
                        
               y1r = cshift(cshift(ur(:,:,:,:,yhat,:,:),dim=xhat,shift=1),
     &                        dim=yhat,shift=-1)
               y1i = cshift(cshift(ui(:,:,:,:,yhat,:,:),dim=xhat,shift=1),
     &                        dim=yhat,shift=-1)

               x1r = cshift(cshift(x2r(:,:,:,:,:,:), dim=xhat, shift=-1),
     &                           dim = yhat, shift=-1)
               x1i = cshift(cshift(x2i(:,:,:,:,:,:), dim=xhat, shift=-1),
     &                           dim = yhat, shift=-1)
            
               call UdagUdag(t1r,t1i,y1r,y1i,x1r,x1i)

               y1r = cshift(cshift(ur(:,:,:,:,yhat,:,:),dim=xhat,shift=-1),
     &                   dim=yhat,shift=-1)
               y1i = cshift(cshift(ui(:,:,:,:,yhat,:,:),dim=xhat,shift=-1),
     &                   dim=yhat,shift=-1)            
        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
c       In the following multiplication we introduce a factor
c       C3 which will serve to create an improved action
c
               t1r = C3*t1r
               t1i = C3*t1i
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UU(t2r,t2i,t1r,t1i,y1r,y1i)
         
               x1r = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)
         
               call UU(stapler,staplei,t2r,t2i,x1r,x1i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
        
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)
c                  write(7,*)'the testactionsq for negative 2 x 1R loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C3 =', C3
c                  write(7,*)

               end if        
            end if                 ! closes the quad == 3,5 condition
         end if                    ! closes the 1x2 and 2x1 if block
c
c===============================================================================
c
c        Constructing the 3 link products in the x and y directions
c

         if (loops == 0 .or. loops == 3 .or. loops == 4 .or. loops == 5 .or. loops == 6) then

c           Preparing to construct the 3 link y-direction product
c           We start by shifting the single link to the third position

            y1r = cshift(ur(:,:,:,:,yhat,:,:), dim=yhat, shift= 2)
            y1i = cshift(ui(:,:,:,:,yhat,:,:), dim=yhat, shift= 2)

c           Constructing the 3 link product from the 2 link product and the cshifted 
c           single link

            y3r = 0.0d0
            y3i = 0.0d0

            call UU(y3r,y3i,y2r,y2i,y1r,y1i)
 
c           Preparing to calculate the 3 link x-product
c           We start by shifting the single link to the third position

            x1r = cshift(ur(:,:,:,:,xhat,:,:), dim=xhat, shift= 2)
            x1i = cshift(ui(:,:,:,:,xhat,:,:), dim=xhat, shift= 2)

c           Calculating the 3 link product from the 2 link product and the cshifted
c           single link

            x3r = 0.0d0
            x3i = 0.0d0

            call UU(x3r,x3i,x2r,x2i,x1r,x1i)

         end if

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Do positive half of the (1x3) staple
c                                                        _
c                                                       | |
c                                                       | |
c                                                       | |

         if (loops == 0 .or. loops == 4 .or. loops == 5 .or. loops == 6) then

            if (quad == 1 .or. quad == 2 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0  
               t2r = 0.0d0
               t2i = 0.0d0 

               x1r = cshift(ur(:,:,:,:,xhat,:,:), dim=yhat, shift= 3)
               x1i = cshift(ui(:,:,:,:,xhat,:,:), dim=yhat, shift= 3)

c        Now we use y1 as a temporary variable, as before

               y1r = cshift(y3r(:,:,:,:,:,:), dim=xhat, shift= 1)
               y1i = cshift(y3i(:,:,:,:,:,:), dim=xhat, shift= 1)

               call UUdag(t1r,t1i,y1r,y1i,x1r,x1i)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
c       In the following multiplication we introduce a factor
c       C4 which will serve to create an improved action
c       
               t1r = C4*t1r
               t1i = C4*t1i        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(stapler,staplei,t1r,t1i,y3r,y3i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:) 
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do   
                  end do   
         
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)
c                  write(7,*)'the testactionsq for positive 1 x 3 loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C4 =', C4
c                  write(7,*)

               end if
            end if                     ! closes the quad == 1,2,5 condition    
  
c ----------------------------------------------------------------------
c     Do negative half of the (1x3) staple
c                                                       | |
c                                                       | |
c                                                       |_|

            if (quad == 3 .or. quad == 4 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0
                
               x1r = cshift(ur(:,:,:,:,xhat,:,:), dim=yhat, shift= -3)
               x1i = cshift(ui(:,:,:,:,xhat,:,:), dim=yhat, shift= -3)

               y1r = cshift(cshift(y3r(:,:,:,:,:,:), dim=xhat, shift= 1),
     &                        dim=yhat, shift=-3)
               y1i = cshift(cshift(y3i(:,:,:,:,:,:), dim=xhat, shift= 1),
     &                        dim=yhat, shift=-3)
               
               call UdagUdag(t1r,t1i,y1r,y1i,x1r,x1i)

               y1r = cshift(y3r(:,:,:,:,:,:), dim=yhat, shift=-3)
               y1i = cshift(y3i(:,:,:,:,:,:), dim=yhat, shift=-3)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
c       In the following multiplication we introduce a factor
c       C4 which will serve to create an improved action
c       
               t1r = C4*t1r
               t1i = C4*t1i        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UU(stapler,staplei,t1r,t1i,y1r,y1i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
         
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for negative 1 x 3 loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C4 =', C4
c                  write(7,*)

               end if
            end if                        ! closes the quad == 3,4,5 condition

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Do positive half of the (3x1) staple, update link at left

c                                                        ___
c                                                       | __|
c        Now we use x1 as a temporary variable as above

            if (quad == 1 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0

               y1r = cshift(ur(:,:,:,:,yhat,:,:), dim=xhat, shift= 3)
               y1i = cshift(ui(:,:,:,:,yhat,:,:), dim=xhat, shift= 3)

               x1r = cshift(x2r(:,:,:,:,:,:), dim = xhat, shift=1)
               x1i = cshift(x2i(:,:,:,:,:,:), dim = xhat, shift=1)
  
               linkr(:,:,:,:,:,:) = ur(:,:,:,:,yhat,:,:)
               linki(:,:,:,:,:,:) = ui(:,:,:,:,yhat,:,:)

               call UU(t2r,t2i,x1r,x1i,y1r,y1i)

               x1r = cshift(x3r(:,:,:,:,:,:), dim=yhat, shift= 1)
               x1i = cshift(x3i(:,:,:,:,:,:), dim=yhat, shift= 1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
c       In the following multiplication we introduce a factor
c       C4 which will serve to create an improved action
c       
               t2r = C4*t2r
               t2i = C4*t2i        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t1r,t1i,t2r,t2i,x1r,x1i) 

               call UUdag(stapler,staplei,t1r,t1i,linkr,linki)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
            
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for positive 3 x 1L loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C4 =', C4         
c                  write(7,*)

               end if
            end if                        ! closes the quad == 1,5 condition

c -----------------------------------------------------------------------
c     Do negative half of the (3x1) staple, update link at left

c                                                         __
c                                                       |___|
            if (quad == 4 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0
        
               y1r = cshift(cshift(ur(:,:,:,:,yhat,:,:), dim=xhat, shift= 3),
     &                        dim=yhat, shift=-1)
               y1i = cshift(cshift(ui(:,:,:,:,yhat,:,:), dim=xhat, shift= 3),
     &                        dim=yhat, shift=-1)
                
               x1r = cshift(x2r(:,:,:,:,:,:), dim = xhat, shift=1)  
               x1i = cshift(x2i(:,:,:,:,:,:), dim = xhat, shift=1)  

               call UUdag(t2r,t2i,x1r,x1i,y1r,y1i)

               x1r = cshift(x3r(:,:,:,:,:,:), dim=yhat, shift= -1)
               x1i = cshift(x3i(:,:,:,:,:,:), dim=yhat, shift= -1)
         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
c       In the following multiplication we introduce a factor
c       C4 which will serve to create an improved action
c       
               t2r = C4*t2r
               t2i = C4*t2i        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t1r,t1i,t2r,t2i,x1r,x1i)
                
               y1r = cshift(ur(:,:,:,:,yhat,:,:), dim=yhat, shift=-1) 
               y1i = cshift(ui(:,:,:,:,yhat,:,:), dim=yhat, shift=-1) 

               call UU(stapler,staplei,t1r,t1i,y1r,y1i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
         
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for negative 3 x 1L loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C4 =', C4
c                  write(7,*)

               end if
            end if                           ! closes the quad == 4,5 condition
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Do positive half of the (3x1) staple, update link at middle

c                                                        ___
c                                                       |_ _|
            if (quad == 5) then
   
               t1r = 0.0d0  
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0

               x1r = cshift(ur(:,:,:,:,xhat,:,:), dim=xhat, shift= 1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:), dim=xhat, shift= 1)

               y1r = cshift(ur(:,:,:,:,yhat,:,:), dim=xhat, shift=2)
               y1i = cshift(ui(:,:,:,:,yhat,:,:), dim=xhat, shift=2)

               call UU(t1r,t1i,x1r,x1i,y1r,y1i)

               x1r = cshift(cshift(x3r(:,:,:,:,:,:), dim=yhat, shift=1), 
     &                        dim=xhat, shift= -1)
               x1i = cshift(cshift(x3i(:,:,:,:,:,:), dim=yhat, shift=1),
     &                        dim=xhat, shift= -1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
c       In the following multiplication we introduce a factor
c       C4 which will serve to create an improved action
c       
               t1r = C4*t1r
               t1i = C4*t1i        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t2r,t2i,t1r,t1i,x1r,x1i)         

               y1r = cshift(ur(:,:,:,:,yhat,:,:), dim=xhat, shift=-1)
               y1i = cshift(ui(:,:,:,:,yhat,:,:), dim=xhat, shift=-1)
 
               t1r = 0.0d0
               t1i = 0.0d0  

               call UUdag(t1r,t1i,t2r,t2i,y1r,y1i)

               x1r = cshift(ur(:,:,:,:,xhat,:,:), dim = xhat, shift = -1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:), dim = xhat, shift = -1)

               call UU(stapler,staplei,t1r,t1i,x1r,x1i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
         
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for positive 3 x 1M loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C4 =', C4
c                  write(7,*)

               end if
           
c ---------------------------------------------------------------------
c     Do negative half of the (3x1) staple, update link at middle
c                                                        _ _
c                                                       |___|
               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0
            
               x1r = cshift(ur(:,:,:,:,xhat,:,:), dim=xhat, shift= 1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:), dim=xhat, shift= 1)
               
               y1r = cshift(cshift(ur(:,:,:,:,yhat,:,:), dim=xhat, shift= 2),
     &                   dim=yhat, shift=-1)
               y1i = cshift(cshift(ui(:,:,:,:,yhat,:,:), dim=xhat, shift= 2),
     &                   dim=yhat, shift=-1)        

               call UUdag(t1r,t1i,x1r,x1i,y1r,y1i)
      
               x1r = cshift(cshift(x3r(:,:,:,:,:,:), dim=yhat, shift=-1),
     &              dim=xhat, shift= -1)
               x1i = cshift(cshift(x3i(:,:,:,:,:,:), dim=yhat, shift=-1),
     &              dim=xhat, shift= -1)
         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
c       In the following multiplication we introduce a factor
c       C4 which will serve to create an improved action
c       
               t1r = C4*t1r
               t1i = C4*t1i        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t2r,t2i,t1r,t1i,x1r,x1i)
                
               y1r = cshift(cshift(ur(:,:,:,:,yhat,:,:), dim=xhat, shift= -1),
     &                   dim=yhat, shift=-1)
               y1i = cshift(cshift(ui(:,:,:,:,yhat,:,:), dim=xhat, shift= -1),
     &                   dim=yhat, shift=-1)
      
               t1r = 0.0d0
               t1i = 0.0d0  

               call UU(t1r,t1i,t2r,t2i,y1r,y1i)              

               x1r = cshift(ur(:,:,:,:,xhat,:,:), dim = xhat, shift = -1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:), dim = xhat, shift = -1)

               call UU(stapler,staplei,t1r,t1i,x1r,x1i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
         
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for negative 3 x 1M loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C4 =', C4
c                  write(7,*)

               end if
            end if                      ! closes the quad == 5 condition
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Do positive half of the (3x1) staple, update link at right
c                                                        ___
c                                                       |__ |

            if (quad == 2 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0
                         
               x1r = cshift(cshift(x3r(:,:,:,:,:,:), dim=yhat, shift=1),
     &              dim=xhat, shift= -2)
               x1i = cshift(cshift(x3i(:,:,:,:,:,:), dim=yhat, shift=1),
     &              dim=xhat, shift= -2)

               y1r = cshift(ur(:,:,:,:,yhat,:,:), dim=xhat, shift=1)
               y1i = cshift(ui(:,:,:,:,yhat,:,:), dim=xhat, shift=1)

               call UUdag(t1r,t1i,y1r,y1i,x1r,x1i)

               y1r = cshift(ur(:,:,:,:,yhat,:,:), dim=xhat, shift=-2)
               y1i = cshift(ui(:,:,:,:,yhat,:,:), dim=xhat, shift=-2)
        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
c       In the following multiplication we introduce a factor
c       C4 which will serve to create an improved action
c       
               t1r = C4*t1r
               t1i = C4*t1i        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t2r,t2i,t1r,t1i,y1r,y1i)
   
               x1r = cshift(x2r(:,:,:,:,:,:), dim = xhat, shift=-2)
               x1i = cshift(x2i(:,:,:,:,:,:), dim = xhat, shift=-2)

               call UU(stapler,staplei,t2r,t2i,x1r,x1i)              

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
         
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for positive 3 x 1R loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C4 =', C4
c                  write(7,*)

               end if
            end if                         ! closes the quad == 2,5 condition
c
c --------------------------------------------------------------------
c     Do negative half of the (3x1) staple, update link at right
c                                                        __ 
c                                                       |___|
            if (quad == 3 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0
                         
               x1r = cshift(cshift(x3r(:,:,:,:,:,:), dim=yhat, shift=-1),
     &              dim=xhat, shift= -2)
               x1i = cshift(cshift(x3i(:,:,:,:,:,:), dim=yhat, shift=-1),
     &              dim=xhat, shift= -2)

               y1r = cshift(cshift(ur(:,:,:,:,yhat,:,:), dim=xhat, shift= 1),
     &                   dim=yhat, shift=-1)
               y1i = cshift(cshift(ui(:,:,:,:,yhat,:,:), dim=xhat, shift= 1),
     &                   dim=yhat, shift=-1)
               
               call UdagUdag(t1r,t1i,y1r,y1i,x1r,x1i)

               y1r = cshift(cshift(ur(:,:,:,:,yhat,:,:), dim=xhat, shift= -2),
     &                   dim=yhat, shift=-1)
               y1i = cshift(cshift(ui(:,:,:,:,yhat,:,:), dim=xhat, shift= -2),
     &                   dim=yhat, shift=-1)        

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
c       In the following multiplication we introduce a factor
c       C4 which will serve to create an improved action
c       
               t1r = C4*t1r
               t1i = C4*t1i        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UU(t2r,t2i,t1r,t1i,y1r,y1i)

               x1r = cshift(x2r(:,:,:,:,:,:), dim = xhat, shift=-2)
               x1i = cshift(x2i(:,:,:,:,:,:), dim = xhat, shift=-2)
                
               call UU(stapler,staplei,t2r,t2i,x1r,x1i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
         
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for negative 3 x 1R loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C4 =', C4         
c                  write(7,*)

               end if
            end if                        !closes the quad == 3,5 condition
         end if                        !closes the 1x3 block
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     Do positive half of the (3x3) staple, update link at left
c
c                                                        ___
c                                                       |   |
c                                                       |   |
c                                                       | __|
c
         if (loops == 0 .or. loops == 3 .or. loops == 5 .or. loops == 6) then

            if (quad == 1 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0

               y1r = cshift(y3r(:,:,:,:,:,:), dim=xhat, shift= 3)
               y1i = cshift(y3i(:,:,:,:,:,:), dim=xhat, shift= 3)

               x1r = cshift(x2r(:,:,:,:,:,:), dim = xhat, shift=1)
               x1i = cshift(x2i(:,:,:,:,:,:), dim = xhat, shift=1)

               call UU(t2r,t2i,x1r,x1i,y1r,y1i)

               x1r = cshift(x3r(:,:,:,:,:,:), dim=yhat, shift= 3)
               x1i = cshift(x3i(:,:,:,:,:,:), dim=yhat, shift= 3)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c       In the following multiplication we introduce a factor
c       C5 which will serve to create an improved action  
c
               t2r = C5*t2r
               t2i = C5*t2i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t1r,t1i,t2r,t2i,x1r,x1i)

               call UUdag(stapler,staplei,t1r,t1i,y3r,y3i)            

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
                  
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for positive 3 x 3L loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C5 =', C5         
c                  write(7,*)

               end if
            end if                     ! closes the quad == 1,5 condition
c
c ----------------------------------------------------------------------
c     Do negative half of the (3x3) staple, update link at left
c
c                                                         __
c                                                       |   |
c                                                       |   |
c                                                       |___|
            if (quad == 4 .or. quad == 5) then
   
               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0

               y1r = cshift(cshift(y3r(:,:,:,:,:,:), dim=xhat, shift= 3),
     &                dim=yhat, shift=-3)
               y1i = cshift(cshift(y3i(:,:,:,:,:,:), dim=xhat, shift= 3),
     &                dim=yhat, shift=-3)

               x1r = cshift(x2r(:,:,:,:,:,:), dim = xhat, shift=1)
               x1i = cshift(x2i(:,:,:,:,:,:), dim = xhat, shift=1)

               call UUdag(t2r,t2i,x1r,x1i,y1r,y1i)

               x1r = cshift(x3r(:,:,:,:,:,:), dim=yhat, shift= -3)
               x1i = cshift(x3i(:,:,:,:,:,:), dim=yhat, shift= -3)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c       In the following multiplication we introduce a factor
c       C5 which will serve to create an improved action  
c
               t2r = C5*t2r
               t2i = C5*t2i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t1r,t1i,t2r,t2i,x1r,x1i)

               y1r = cshift(y3r(:,:,:,:,:,:), dim=yhat, shift=-3) 
               y1i = cshift(y3i(:,:,:,:,:,:), dim=yhat, shift=-3) 

               call UU(stapler,staplei,t1r,t1i,y1r,y1i)       
 
               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
         
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for negative 3 x 3L loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C5 =', C5         
c                  write(7,*)

               end if
            end if                        ! closes the quad == 4,5 condition
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Do positive half of the (3x3) staple, update link at middle
c                                                        ___
c                                                       |   |
c                                                       |   |
c                                                       |_ _|
            if (quad == 5) then

               t1r = 0.0d0  
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0
         
               x1r = cshift(ur(:,:,:,:,xhat,:,:), dim=xhat, shift= 1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:), dim=xhat, shift= 1)

               y1r = cshift(y3r(:,:,:,:,:,:), dim=xhat, shift=2)
               y1i = cshift(y3i(:,:,:,:,:,:), dim=xhat, shift=2)

               call UU(t1r,t1i,x1r,x1i,y1r,y1i)

               x1r = cshift(cshift(x3r(:,:,:,:,:,:), dim=yhat, shift=3), 
     &                 dim=xhat, shift= -1)
               x1i = cshift(cshift(x3i(:,:,:,:,:,:), dim=yhat, shift=3),
     &                 dim=xhat, shift= -1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c       In the following multiplication we introduce a factor
c       C5 which will serve to create an improved action  
c
               t1r = C5*t1r
               t1i = C5*t1i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t2r,t2i,t1r,t1i,x1r,x1i)

               y1r = cshift(y3r(:,:,:,:,:,:), dim=xhat, shift=-1)
               y1i = cshift(y3i(:,:,:,:,:,:), dim=xhat, shift=-1)

               t1r = 0.0d0
               t1i = 0.0d0

               call UUdag(t1r,t1i,t2r,t2i,y1r,y1i)

               x1r = cshift(ur(:,:,:,:,xhat,:,:), dim = xhat, shift = -1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:), dim = xhat, shift = -1)

               call UU(stapler,staplei,t1r,t1i,x1r,x1i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
         
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for positive 3 x 3M loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C5 =', C5         
c                  write(7,*)

               end if

c ---------------------------------------------------------------------
c     Do negative half of the (3x3) staple, update link at middle
c                                                        _ _
c                                                       |   |
c                                                       |   |
c                                                       |___|
               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0

               x1r = cshift(ur(:,:,:,:,xhat,:,:), dim=xhat, shift= 1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:), dim=xhat, shift= 1)

               y1r = cshift(cshift(y3r(:,:,:,:,:,:), dim=xhat, shift= 2),
     &                   dim=yhat, shift=-3) 
               y1i = cshift(cshift(y3i(:,:,:,:,:,:), dim=xhat, shift= 2),
     &                   dim=yhat, shift=-3)

               call UUdag(t1r,t1i,x1r,x1i,y1r,y1i)

               x1r = cshift(cshift(x3r(:,:,:,:,:,:), dim=yhat, shift=-3), 
     &                dim=xhat, shift= -1)
               x1i = cshift(cshift(x3i(:,:,:,:,:,:), dim=yhat, shift=-3),
     &                dim=xhat, shift= -1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c       In the following multiplication we introduce a factor
c       C5 which will serve to create an improved action  
c
               t1r = C5*t1r
               t1i = C5*t1i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t2r,t2i,t1r,t1i,x1r,x1i)

               y1r = cshift(cshift(y3r(:,:,:,:,:,:), dim=xhat, shift= -1),
     &                   dim=yhat, shift=-3) 
               y1i = cshift(cshift(y3i(:,:,:,:,:,:), dim=xhat, shift= -1),
     &                   dim=yhat, shift=-3)

               t1r = 0.0d0
               t1i = 0.0d0

               call UU(t1r,t1i,t2r,t2i,y1r,y1i)

               x1r = cshift(ur(:,:,:,:,xhat,:,:), dim = xhat, shift = -1)
               x1i = cshift(ui(:,:,:,:,xhat,:,:), dim = xhat, shift = -1)

               call UU(stapler,staplei,t1r,t1i,x1r,x1i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
         
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for negative 3 x 3M loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C5 =', C5         
c                  write(7,*)

               end if
            end if                        ! closes the quad == 5 condition
c   
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Do positive half of the (3x3) staple, update link at right
c                                                        ___
c                                                       |   |
c                                                       |   |
c                                                       |__ |
            if (quad == 2 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0

               x1r = cshift(cshift(x3r(:,:,:,:,:,:), dim=yhat, shift=3),
     &              dim=xhat, shift= -2)
               x1i = cshift(cshift(x3i(:,:,:,:,:,:), dim=yhat, shift=3),
     &              dim=xhat, shift= -2)

               y1r = cshift(y3r(:,:,:,:,:,:), dim=xhat, shift=1)
               y1i = cshift(y3i(:,:,:,:,:,:), dim=xhat, shift=1)
 
               call UUdag(t1r,t1i,y1r,y1i,x1r,x1i)

               y1r = cshift(y3r(:,:,:,:,:,:), dim=xhat, shift=-2)
               y1i = cshift(y3i(:,:,:,:,:,:), dim=xhat, shift=-2)
         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c       In the following multiplication we introduce a factor
c       C5 which will serve to create an improved action  
c
               t1r = C5*t1r
               t1i = C5*t1i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UUdag(t2r,t2i,t1r,t1i,y1r,y1i)

               x1r = cshift(x2r(:,:,:,:,:,:), dim = xhat, shift=-2)
               x1i = cshift(x2i(:,:,:,:,:,:), dim = xhat, shift=-2)

               call UU(stapler,staplei,t2r,t2i,x1r,x1i)               

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
         
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for positive 3 x 3R loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C5 =', C5         
c                  write(7,*)

               end if
            end if                        ! closes the quad == 2,5 condition        
c
c ------------------------------------------------------------------------
c     Do negative half of the (3x3) staple, update link at right
c                                                        __
c                                                       |   |
c                                                       |   |
c                                                       |___|
            if (quad == 3 .or. quad == 5) then

               t1r = 0.0d0
               t1i = 0.0d0
               t2r = 0.0d0
               t2i = 0.0d0
         
               x1r = cshift(cshift(x3r(:,:,:,:,:,:), dim=yhat, shift=-3),
     &                   dim=xhat, shift= -2)
               x1i = cshift(cshift(x3i(:,:,:,:,:,:), dim=yhat, shift=-3),
     &                   dim=xhat, shift= -2)

               y1r = cshift(cshift(y3r(:,:,:,:,:,:), dim=xhat, shift= 1),
     &                   dim=yhat, shift=-3) 
               y1i = cshift(cshift(y3i(:,:,:,:,:,:), dim=xhat, shift= 1),
     &                   dim=yhat, shift=-3)

               call UdagUdag(t1r,t1i,y1r,y1i,x1r,x1i)

               y1r = cshift(cshift(y3r(:,:,:,:,:,:), dim=xhat, shift= -2),
     &                   dim=yhat, shift=-3) 
               y1i = cshift(cshift(y3i(:,:,:,:,:,:), dim=xhat, shift= -2),
     &                   dim=yhat, shift=-3)
        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c       In the following multiplication we introduce a factor
c       C5 which will serve to create an improved action  
c
               t1r = C5*t1r
               t1i = C5*t1i        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               call UU(t2r,t2i,t1r,t1i,y1r,y1i)

               x1r = cshift(x2r(:,:,:,:,:,:), dim = xhat, shift=-2)
               x1i = cshift(x2i(:,:,:,:,:,:), dim = xhat, shift=-2)
                
               call UU(stapler,staplei,t2r,t2i,x1r,x1i)

               if (report) then

                  testactionsq = 0.0d0

                  do ic=1,nc
                     do kc=1,nc

                        testactionsq(:,:,:,:)    = testactionsq(:,:,:,:)
     &                       + ( ur(:,:,:,:,xhat,ic,kc) * stapler(:,:,:,:,kc,ic)
     &                       -   ui(:,:,:,:,xhat,ic,kc) * staplei(:,:,:,:,kc,ic) )

                     end do
                  end do
         
                  testactionsq = sum(testactionsq) / ( nx*ny*nz*nt )
c                  write(7,*)  
c                  write(7,*)'the testactionsq for negative 3 x 3R loop, i , yhat'
c                  write(7,*) testactionsq(1,1,1,1),yhat
c                  write(7,*) 'C5 =', C5         
c                  write(7,*)

               end if
            end if              ! closes the quad == 3,5 condition
         end if                 ! closes if (loops == 0 .or. loops == 3 .or. loops == 5 .or. loops == 6)

         return

         end subroutine LotsaLoops

      END MODULE L_LotsaLoops
