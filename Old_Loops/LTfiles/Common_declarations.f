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
