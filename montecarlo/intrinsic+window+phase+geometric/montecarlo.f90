PROGRAM sim

implicit none

 INTEGER :: m
 INTEGER, PARAMETER :: mmax = 1E6
 REAL(8), PARAMETER :: baseline = 27.4D0
 REAL(8), PARAMETER :: Pmin = 1.0D0
 REAL(8), PARAMETER :: Pmax = 1000.0D0
 REAL(8), PARAMETER :: alpha = -0.66666667D0
 REAL(8), PARAMETER :: rhoSun = 1411.0390435974123D0
 REAL(8), PARAMETER :: pi = 3.141592653589793D0
 REAL(8), PARAMETER :: third = 0.333333333D0
 REAL(8), PARAMETER :: Grv = 6.67408D-11
 REAL(8) :: logPmin, logPmax
 REAL(8) :: ran
 REAL(8) :: period, t0
 REAL(8) :: rhostar, inc, b, aR
 LOGICAL :: single, transit
 REAL(8), DIMENSION(mmax) :: Pp
 
 ! initialize
 call randomreal(.TRUE.,ran)
 logPmin = DLOG(Pmin)
 logPmax = DLOG(Pmax)
 rhostar = 1.0D0*rhoSun
 
 DO WHILE ( m .LT. mmax )
 
   ! generate a random period between Pmin and Pmax
   call randomreal(.FALSE.,ran)
	 IF( alpha .EQ. -1.0D0 ) THEN
     period = DEXP( logPmin + (logPmax - logPmin)*ran )
	 ELSE
		 period = ran*Pmax**(1.0D0+alpha) + (1.0D0-ran)*Pmin**(1.0D0+alpha)
		 period = period**( 1.0D0/(1.0D0+alpha) )
	 END IF
	 
   ! calculate a/R*
   aR = ( Grv*(period*86400.0D0)**2*rhoSun )/( 3.0D0*pi )
   aR = aR**third
 
	 ! generate a random phase
	 call randomreal(.FALSE.,ran)
	 t0 = period*ran
	 
   ! generate a random inclination
   call randomreal(.FALSE.,ran)
   inc = DACOS( 2.0D0*ran - 1.0D0 )
	 
   ! calculate impact parameter
   b = aR*DCOS(inc)
	 
   ! does it transit?
   IF( DABS(b) .LT. 1.0D0 ) THEN
	   transit = .TRUE.
   ELSE
	   transit = .FALSE.
   END IF
	 
   ! does it single transit in the baseline?
   IF( t0 .LT. baseline .AND. (t0+period) .GT. baseline .AND. (t0-period) .LT. 0.0D0 ) THEN
		 single = .TRUE.
	 ELSE
		 single = .FALSE.
	 END IF

   ! save the result
	 IF( single .AND. transit ) THEN
		 IF( t0 .GT. 4.95D0 .AND. t0 .LT. 5.05D0 ) THEN
		   m = m + 1
		   Pp(m) = period
		 END IF
	 END IF
	 
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