PROGRAM sim

implicit none

 INTEGER :: m
 INTEGER, PARAMETER :: mmax = 1E6
 REAL(8), PARAMETER :: Pmin = 1.0D0
 REAL(8), PARAMETER :: Pmax = 1000.0D0
 REAL(8), PARAMETER :: alpha = -0.66666667D0
 REAL(8) :: logPmin, logPmax
 REAL(8) :: ran
 REAL(8), DIMENSION(mmax) :: Pp
 
 ! initialize
 call randomreal(.TRUE.,ran)
 logPmin = DLOG(Pmin)
 logPmax = DLOG(Pmax)
 
 DO WHILE ( m .LT. mmax )
 
   ! generate a random period between Pmin and Pmax
   call randomreal(.FALSE.,ran)
	 IF( alpha .EQ. -1.0D0 ) THEN
     period = DEXP( logPmin + (logPmax - logPmin)*ran )
	 ELSE
		 period = ran*Pmax**(1.0D0+alpha) + (1.0D0-ran)*Pmin**(1.0D0+alpha)
		 period = period**( 1.0D0/(1.0D0+alpha) )
	 END IF

   ! save the result
	 m = m + 1
	 Pp(m) = period
	 
 END DO
 
 ! export the results
 OPEN(unit=103,file='periods.dat')
 DO m=1,mmax
	 write(103,*) Pp(m)
 END DO
 CLOSE(103)
 
CONTAINS
	
  ! =======================================================
   SUBROUTINE randomreal(firsttime,r)

   INTEGER :: count
   LOGICAL :: firsttime
   REAL(8) :: r

   IF( firsttime ) THEN
     CALL SYSTEM_CLOCK(count)
     CALL srand(count)
     r = rand()
   ELSE
     r = rand()
   END IF

   END SUBROUTINE randomreal
  ! =======================================================

END PROGRAM sim