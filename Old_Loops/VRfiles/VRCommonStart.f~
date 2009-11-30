      pi = 4.0d0 * atan(1.0d0)
      volume = nx * ny * nz * nt

      write(*,*)
      write(*,'(a)')'Please enter the configurations directory.'
      read (*,'(a132)') basedir
      write(*,'(a)') basedir(1:len_trim(basedir))

      write(*,*)
      write(*,'(a)')'Please enter the gauge field file name.'
      read (*,'(a132)') lastConfig
      write(*,'(a)') lastConfig(1:len_trim(lastConfig))

      write(*,*)
      write(*,'(a)')'Please enter the index of the first configuration.'
      read (*,*) incfg
      write(*,*) incfg

      write(*,*)
      write(*,'(a)')'Please enter the index of the last configuration.'
      read (*,*) ficfg
      write(*,*) ficfg

      write(*,*)
      write(*,'(a)')'Please enter directory to output to. (baryonCorrel/) .'
      write(*,'(a)')'Note: this subdirectory must exist.'
      read (*,'(a132)') prefix
      write(*,'(a)') prefix(1:len_trim(prefix))

      write(*,*)
      write(*,'(a)')'Which F_mu_nu do you desire? (3)'
      write(*,'(a)')'      1: Standard plaquette F_mu_nu'
      write(*,'(a)')'      2: (1x1), (1x2) F_mu_nu'
      write(*,'(a)')'      3: (1x1), (2x2), (3x3) F_mu_nu'
      write(*,'(a)')'      4: (1x1), (2x2), (1x2),(1X3) F_mu_nu (4LTQ #1)'
      write(*,'(a)')'      5: 5LI F_mu_nu'
      read (*,*) Qpaths
      write(*,*) Qpaths

      write(*,*)
      write(*,'(a)')'Please enter the amount of spatial smearing.'
      read (*,*) smear3d
      write(*,*) smear3d

      write(*,*)
      write(*,'(a)')'What would you like to correlate? (1)'
      write(*,'(a)')'      1: Action and topological charge.'
      write(*,'(a)')'      2: Electric and Magnetic fields.'
      read (*,*) Correlate
      write(*,*) Correlate
      
      write(*,*)'Smear Links For The FMuNu Subroutine? (.true. or .false.) (.true.)'
      read(*,'(l7)') smear_links
      write(*,'(l7)') smear_links
      write(*,*)
c'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialize the epsilon index for the calls from lt_oop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call setEpsilonIndex()        

      do cfg = incfg, ficfg

         CALL SYSTEM_CLOCK(start_count, count_rate)

         write(input,fmt='(a,a,i3.3)')trim(basedir),trim(lastconfig),cfg

         call ReadLinks(input,ur,ui,nfig,beta,nxold,nyold,nzold,ntold,lastPlaq,plaqbarAvg,uzero)

         if (nx.ne.nxold) then
            write(*,'(a)') 'Mismatch in nx.'
            goto 999
         else if (ny.ne.nyold) then
            write(*,'(a)') 'Mismatch in ny.'
            goto 999
         else if (nz.ne.nzold) then
            write(*,'(a)') 'Mismatch in nz.'
            goto 999
         else if (nt.ne.ntold) then
            write(*,'(a)') 'Mismatch in nt.'
            goto 999
         else if ( nfig.eq.0 ) then
            write(*,'(a)') 'Mismatch in average plaquette.'
            go to 999
         end if

         if ( smear_links ) then

            write (output,'(a,a,i3.3,a,i2.2,a,i2.2)')trim(prefix),trim(lastconfig),cfg,'.s',nsteps,'-s',smear3d

         else

            write(output,'(a,a,i3.3,a,i2.2)')trim(prefix),trim(lastconfig),cfg,'.s00-s',smear3d

         endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  The gauge field has been successfully read
c  Get the field strength tensor
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         utr = ur
         uti = ui
         if ( smear_links ) then
            do istep = 1,nsteps 
         
               call ape_smear(utr,uti,alpha)
               
            end do
         end if

         call tadpoleimp(utr,uti,uzero)      
         call fMuNu(utr,uti,Fr,Fi,uzero,Qpaths)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         if (Correlate == 1) then

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Calculation of topological charge density
c  Scale results to single instanton results
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            call TopQandReconAction(Fr,Fi,uzero,topQDensity,Q,reconAction,pluckbar)

            sons1 = beta * pluckbar * nx * ny * nz * nt * (mu*(mu-1)/2.0d0) / 
     &           ( 8.0d0*pi**2/(6.0d0/beta) )

            reconAction = reconAction * beta * (mu*(mu-1)/2.0d0) / ( 8.0d0*pi**2/(6.0d0/beta) )

            write(*,'(/,2(a,f15.8,/),/)') 
     &           'Reconstructed S/S_0 = ', sons1,
     &           'Topological Charge  = ', Q

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         else

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  We'll put the electric field in the reconAction and the
c                magnetic field in the topQdensity
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            call EleMagField(Fr,Fi,reconAction,topQdensity)

            write(*,'(/,2(a,e15.8,/),/)') 
     &           'Sum of Electric Field = ', SUM(reconAction), 
     &           'Sum of Magnetic Field = ', SUM(topQdensity)

         end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Get the Wilson Loops
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         utr = ur
         uti = ui

         do istep = 1,smear3d
            
            call ape_smear3D(utr,uti,alpha,that)

         end do

         if ( Correlate == 1) then
               
            fstr1='.action.correl.unf'
            fstr2='.topChg.correl.unf'
            
         else

            fstr1='.ele.correl.unf'
            fstr2='.mag.correl.unf'
            
         endif
