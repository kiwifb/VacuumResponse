!----------------------------------------------------------------------!
!                                                                      !
!     Initialize the intrinsic Fortran random number generator         !
!     using a single integer seed.                                     !
!                                                                      !
!                                                                      !
!     Fortran 90/95 has an intrinsic function RANDOM_NUMBER for        !
!     generating pseudo-random numbers. There is also an intrinsic     !
!     function RANDOM_SEED for initializing the random number          !
!     generator.                                                       !
!                                                                      !
!     The problem with the RANDOM_SEED intrinsic function is that      !
!     the number of seeds it takes as input parameters depends on      !
!     how the random number generator is implemented, so it will       !
!     vary with different Fortran compilers. It is up to the user      !
!     to check how many seeds are required and then to provide the     !
!     appropriate number of *random* seeds -- choosing correlated      !
!     seeds such as seed, seed+1, seed+2, etc will affect the          !
!     randomness properties of most random number generators.          !
!                                                                      !
!     Of course this means that you need to use a random number        !
!     generator in order to initialize the random number generator,    !
!     and that is what is provided here.                               !
!                                                                      !
!     The routine HPF_RANDOM_SEED takes a single seed as input,        !
!     checks the Fortran 90 intrinsic RANDOM_SEED to see how many      !
!     seeds are required, and then calls the routine bitrand_seeder,   !
!     which uses a simple linear congruential random number            !
!     generator to create the appropriate number of random seeds,      !
!     using random values for each bit of each seed. The seeds are     !
!     passed to RANDOM_SEED to initialize the RANDOM_NUMBER            !
!     intrinsic function.                                              !
!                                                                      !
!----------------------------------------------------------------------!
      MODULE HPFRANDOMSEED

      CONTAINS

      SUBROUTINE HPF_RANDOM_SEED(seed)
      IMPLICIT NONE
      INTEGER seed
      INTEGER N
      INTEGER, DIMENSION(:), ALLOCATABLE :: seeder_array

!--- Find out the size of the seed array used in the Fortran 90 
!--- random number generator
      CALL RANDOM_SEED(size=N)

!--- Initialize the seed array by setting each bit using a 
!--- simple linear congruential generator (in bitrand_seeder)
      ALLOCATE ( seeder_array(N) ) 
      CALL bitrand_seeder(seed,N,seeder_array)
      CALL RANDOM_SEED( put=seeder_array(1:N) )
      DEALLOCATE( seeder_array )

      END SUBROUTINE


!----------------------------------------------------------------------!
!                                                                      !
!     Initialize each bit of the 32-bit integer seed array using a     !
!     simple linear congruential generator.                            !
!                                                                      !
!----------------------------------------------------------------------!

      SUBROUTINE bitrand_seeder(seed,N,seed_array)
      IMPLICIT NONE
      INTEGER seed,N
      INTEGER seed_array(N)
      INTEGER i,j,k
      INTEGER s,t,randx
      DOUBLE PRECISION TWOTO31, TWOTONEG31
      PARAMETER (TWOTO31 = 2147483648.0)   ! 2^{31} 
      PARAMETER (TWOTONEG31=(1.0/TWOTO31)) 
      INTEGER A, C, M
      PARAMETER (A=1103515245, C=12345)
      DATA  M   /Z'7FFFFFFF'/

!--- Initialize the linear congruential generator RAND.
      randx = seed

      DO i = 1, N
        s = 0
        t = 1
        DO j = 0, 30
          randx = IAND((A * randx + C),M)  
          IF((randx * TWOTONEG31) .LT. 0.5) THEN
            s = s + t
          ENDIF
           t = 2 * t
!------- Throw away a few random numbers just to foil any 
!------- possible correlations that might affect RAND.
          DO k=1,5
            randx = IAND((A * randx + C),M)     
          ENDDO
        ENDDO
!------- Throw away a few random numbers just to foil any 
!------- possible correlations that might affect RAND.
        DO k=1,13
          randx = IAND((A * randx + C),M)     
        ENDDO

        seed_array(i) = s

      ENDDO         

      RETURN
      END SUBROUTINE

      END MODULE HPFRANDOMSEED
