      MODULE L_epsilonIndex
c
c  creates the color permutation index vectors
c
        implicit none
        integer          :: la(36),lb(36),lc(36),lap(36),lbp(36),lcp(36)
        double precision :: fper(36)

        CONTAINS

          subroutine setEpsilonIndex()

            integer :: i,ii,j,jj

            DO jj=1,2
            IF (jj.EQ.1) j=0
            IF (jj.EQ.2) j=9
            DO i=1,3
              la(i+j)=1
              la(3+i+j)=2
              la(6+i+j)=3
              la(18+i+j)=2
              la(21+i+j)=1
              la(24+i+j)=3
              lb(i+j)=2
              lb(3+i+j)=3
              lb(6+i+j)=1
              lb(18+i+j)=1
              lb(21+i+j)=3
              lb(24+i+j)=2
              lc(i+j)=3
              lc(3+i+j)=1
              lc(6+i+j)=2
              lc(18+i+j)=3
              lc(21+i+j)=2
              lc(24+i+j)=1
            ENDDO
            ENDDO

            DO jj=1,2
              IF (jj.EQ.1) j=0
              IF (jj.EQ.2) j=18
              DO ii=1,3
                i=(ii-1)*3
                lap(1+i+j)=1
                lap(2+i+j)=2
                lap(3+i+j)=3
                lap(10+i+j)=2
                lap(11+i+j)=1
                lap(12+i+j)=3
                lbp(1+i+j)=2
                lbp(2+i+j)=3
                lbp(3+i+j)=1
                lbp(10+i+j)=1
                lbp(11+i+j)=3
                lbp(12+i+j)=2
                lcp(1+i+j)=3
                lcp(2+i+j)=1
                lcp(3+i+j)=2
                lcp(10+i+j)=3
                lcp(11+i+j)=2
                lcp(12+i+j)=1
              ENDDO
            ENDDO

            DO i=1,9
              fper(i)=1.D0
              fper(i+27)=1.D0
              fper(i+9)=-1.D0
              fper(i+18)=-1.D0
            ENDDO
            RETURN
          END subroutine setEpsilonIndex

      END MODULE L_epsilonIndex


