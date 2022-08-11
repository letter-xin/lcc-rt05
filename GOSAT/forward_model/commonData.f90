MODULE commonData

PUBLIC :: GasInfo, wvnm_resolution, SPLMI, SEVAL,SPLMI4, SEVAL4, LinTerp, LinTerpDouble,&
			 reallocate, reallocate2, BBFN, W4, VOIGT, E1, EN, eb, gasdev

TYPE GasInfo
  INTEGER :: Gas         ! 1:CO2, 2:H2O, 3:CO, 4:CH4, 5:C2H4, 6:NO, 7:OH
  REAL(8) :: P,T,x
END TYPE GasInfo

INTEGER,PARAMETER,PUBLIC :: CO2=1,H2O=2,CO=3,CH4=4,C2H4=5,NO=6,OH=7
REAL(8),PUBLIC :: wvnmst
REAL(8),PUBLIC :: wave_b = 200.d0  ! Beginning wavenumber
REAL(8),PUBLIC :: wave_e = 15000.d0  ! Ending wavenumber

CONTAINS

!------------------------------------------------------------------------------------
SUBROUTINE wvnm_resolution(Gas_Info,wave_switch,eta_b,eta_e,wvnmst,number)
! Based on the gas information, decide the wavenumber resolution
!  Input:
!    Gas_Info     --  information of the gas
!  Output:
!    wvnmst       --  wavenumber step for absco generation
!    number       --  total number of wavenumber locations for absco generation 
  IMPLICIT NONE
    TYPE(GasInfo),INTENT(IN) :: Gas_Info
    LOGICAL,INTENT(IN)  :: wave_switch
    REAL(8),INTENT(IN)  :: eta_b
    REAL(8),INTENT(IN)  :: eta_e
    REAL(8),INTENT(OUT) :: wvnmst
    INTEGER,INTENT(OUT) :: number

!!$ local variables
    REAL(8),PARAMETER :: EPS= 1.d-6
    REAL(8) :: P,x,T, w

    P= Gas_Info%P; x= Gas_Info%x; T= Gas_Info%T

IF(wave_switch) THEN
    IF (P> 1.1d0-EPS) THEN
	  IF(T> 330) THEN; w= 5.d-3
      ELSE; w = 1.d-2; ENDIF

	ELSEIF (P> 1.d0-EPS) THEN
	  IF(T> 290) THEN; w= 5.d-3
	  ELSE; w= 1.d-2; ENDIF

	ELSEIF (P> 0.9d0-EPS) THEN
	  IF(T> 250) THEN; w= 5.d-3
	  ELSE; w= 1.d-2; ENDIF

	ELSEIF (P> 0.8d0-EPS) THEN
	  IF(T> 220) THEN; w=5.d-3
	  ELSE; w=1.d-2; ENDIF				!1.d-2 or 5.d-3

	ELSEIF (P> 0.7d0-EPS) THEN
	  IF(T> 180) THEN; w=5.d-3
	  ELSE; w=1.d-2; ENDIF

	ELSEIF (P> 0.6d0-EPS) THEN
	  IF(T> 380) THEN; w=5.d-4
	  ELSEIF(T> 150) THEN; w= 5.d-3
	  ELSE; w= 1.d-2; ENDIF

	ELSEIF (P> 0.5d0-EPS) THEN
	  IF(T> 300) THEN; w=5.d-4
	  ELSEIF(T> 110) THEN; w= 5.d-3		!5.d-3 or 2.5d-3
	  ELSE; w= 1.d-2; ENDIF

	ELSEIF (P> 0.4d0-EPS) THEN
	  IF(T> 220) THEN; w=5.d-4
	  ELSEIF(T> 80) THEN; w= 5.d-3
	  ELSE; w= 1.d-2; ENDIF

	ELSEIF (P> 0.3d0-EPS) THEN
	  IF(T> 150) THEN; w=5.d-4
	  ELSEIF(T> 60) THEN; w= 5.d-3
	  ELSE; w= 1.d-2; ENDIF

	ELSEIF (P> 0.2d0-EPS) THEN
	  IF(T> 80) THEN; w=5.d-4
	  ELSEIF(T> 30) THEN; w= 5.d-3
	  ELSE; w= 1.d-2; ENDIF

	ELSEIF (P> 0.1d0-EPS) THEN
	  IF(T> 30) THEN; w=5.d-4
	  ELSEIF(T> 10) THEN; w= 5.d-3
	  ELSE; w= 1.d-2; ENDIF

	ELSEIF (P> 0.09d0-EPS) THEN
	  IF(T> 30) THEN; w=5.d-4
	  ELSEIF(T> 10) THEN; w= 5.d-3
	  ELSE; w= 1.d-2; ENDIF

	ELSEIF (P> 0.08d0-EPS) THEN
	  IF(T> 20) THEN; w=5.d-4
	  ELSEIF(T> 10) THEN; w= 5.d-3
	  ELSE; w= 1.d-2; ENDIF

	ELSEIF (P> 0.07d0-EPS) THEN
	  IF(T> 20) THEN; w=5.d-4
	  ELSE; w= 5.d-3; ENDIF

	ELSEIF (P> 0.05d0-EPS) THEN
	  IF(T> 10) THEN; w=5.d-4
	  ELSE; w= 5.d-3; ENDIF

	ELSEIF (P> 0.04d0-EPS) THEN
	  IF(T> 10) THEN; w=5.d-4
	  ELSE; w= 5.d-3; ENDIF

	ELSEIF (P> 0.02d0-EPS) THEN
	  w= 5.d-4

	ELSEIF (P> 0.01d0-EPS) THEN
	  IF(T> 50) THEN; w=1.d-4			!1.d-4 or 5.d-5
	  ELSE; w= 5.d-4; ENDIF				!5.d-4 or 2.5d-4

	ELSE
	w= 1.d-4
	ENDIF
ENDIF

    wvnmst= w
    number=(eta_e-eta_b)/wvnmst+EPS+1
  END SUBROUTINE wvnm_resolution
!!$------------------------------------------------------------------------

   SUBROUTINE LinTerp(Y1, Y2, Yout, XX, x)
   ! This routine performs a linear interpolation
   ! Input: 
   !        Y1     !! wvnmst must be the same for Y1 and Y2 so they are of the same length.
   !        Y2     !! Otherwise, you will be interpolating between abscos of different wavenumbers.
   !        XX     !! array of length 2 containing concentrations to interpolate over.
   !        x      !! Concentration to calc absco at.

   ! Output:
   !        Yout   !! The interpolated absco.

      IMPLICIT NONE
      REAL(4), INTENT(IN) :: Y1, Y2
      REAL(4), INTENT(OUT) :: Yout
      REAL(8), INTENT(IN) :: XX(2),x      ! X(1) = lower value; X(2) = upper value
      REAL(8) :: diff1, diff2

      diff1 = x - XX(1)
      diff2 = XX(2) - XX(1)
      Yout = Y1 + diff1/diff2*(Y2 - Y1)
   END SUBROUTINE

!!$------------------------------------------------------------------------

   SUBROUTINE LinTerpDouble(Y1, Y2, Yout, XX, x)
   ! This routine performs a linear interpolation
   ! Input: 
   !        Y1     !! wvnmst must be the same for Y1 and Y2 so they are of the same length.
   !        Y2     !! Otherwise, you will be interpolating between abscos of different wavenumbers.
   !        XX     !! array of length 2 containing concentrations to interpolate over.
   !        x      !! Concentration to calc absco at.

   ! Output:
   !        Yout   !! The interpolated absco.

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: Y1, Y2
      REAL(8), INTENT(OUT) :: Yout
      REAL(8), INTENT(IN) :: XX(2),x      ! X(1) = lower value; X(2) = upper value
      REAL(8) :: diff1, diff2

      diff1 = x - XX(1)
      diff2 = XX(2) - XX(1)
      Yout = Y1 + diff1/diff2*(Y2 - Y1)
   END SUBROUTINE

!--------------------------------------------------------------------------
SUBROUTINE SPLMI (N, X, Y, B, C, D)
   INTEGER :: N
   REAL(8) :: Y(N)
   DOUBLE PRECISION :: X(N), B(N), C(N), D(N),DELTA(N),H(N),AL(N),BE(N)
   !
   !  THE COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N ARE COMPUTED
   !  FOR A MONOTONICALLY VARYING CUBIC INTERPOLATING SPLINE
   !
   !    S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
   !
   !    FOR  X(I) .LE. X .LE. X(I+1)
   !    WITH Y(I+1).GE.Y(I) (ALL I) OR Y(I+1).LT.Y(I) (ALL I)
   !
   !  INPUT..
   !
   !    N = THE NUMBER OF DATA POINTS OR KNOTS (N.GE.2)
   !    X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
   !    Y = THE ORDINATES OF THE KNOTS
   !
   !  OUTPUT..
   !
   !    B, C, D  = ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.
   !
   !*************************************************************************
   !
   ! THEORY FROM 'MONOTONE PIECEWISE CUBIC INTERPOLATION',
   ! BY F.N. FRITSCH AND R.E. CARLSON IN SIAM J.NUMER.ANAL.,V.17,P.238
   !
   !*************************************************************************
   !
   !    Y(I) = S(X(I))
   !    B(I) = SP(X(I))
   !    C(I) = SPP(X(I))/2
   !    D(I) = SPPP(X(I))/6  (DERIVATIVE FROM THE RIGHT)
   !
   !  THE ACCOMPANYING FUNCTION SUBPROGRAM  SEVAL  CAN BE USED
   !  TO EVALUATE THE SPLINE.
   !
   !
      NM1 = N-1
      IF ( N .LT. 2 ) RETURN
      IF ( N .LT. 3 ) GO TO 100
!
! CALCULATE THE H(I) AND DELTA(I)
!
      DO 10 I=1,NM1
      AL(I)=0.
      BE(I)=0.
      H(I)=X(I+1)-X(I)
   10 DELTA(I)=(Y(I+1)-Y(I))/H(I)
!
! CALCULATE FIRST VALUES FOR AL AND BE BY 3-POINT DIFFERENCE
!
      IF(DELTA(1).EQ.0) GOTO 15
      AL(1)=((H(1)+H(2))**2*Y(2)-H(1)**2*Y(3)-H(2)*(2.*H(1)+H(2))&
             *Y(1))/(H(2)*(H(1)+H(2))*(Y(2)-Y(1)))
   15 DO 20 I=2,NM1
      IF(DELTA(I).EQ.0) GOTO 20
      AL(I)=(H(I-1)**2*Y(I+1)+(H(I)**2-H(I-1)**2)*Y(I)-H(I)**2*&
             Y(I-1))/(H(I-1)*(H(I)+H(I-1))*(Y(I+1)-Y(I)))
   20 CONTINUE
!
      NM2=N-2
      DO 30 I=1,NM2
      IF(DELTA(I).EQ.0.) GOTO 30
      BE(I)=(H(I)**2*Y(I+2)+(H(I+1)**2-H(I)**2)*Y(I+1)-H(I+1)**2*&
             Y(I))/(H(I+1)*(H(I)+H(I+1))*(Y(I+1)-Y(I)))
   30 CONTINUE
!
      IF(DELTA(N-1).EQ.0.) GOTO 35
      BE(N-1)=(H(N-2)*(2.*H(N-1)+H(N-2))*Y(N)-(H(N-1)+H(N-2))**2&
             *Y(N-1)+H(N-1)**2*Y(N-2))/(H(N-2)*(H(N-1)+H(N-2))&
             *(Y(N)-Y(N-1)))
!
! CORRECT VALUES FOR AL AND BE
!
   35 DO 40 I=1,NM1
      IF(AL(I)+BE(I).LE.2.) GOTO 40
      IF(2.*AL(I)+BE(I).LE.3.) GOTO 40
      IF(AL(I)+2.*BE(I).LE.3.) GOTO 40
      PHI=AL(I)-(2.*AL(I)+BE(I)-3.)**2/(AL(I)+BE(I)-2.)/3.
      IF(PHI.GE.0.) GOTO 40
      TI=3./SQRT(AL(I)**2+BE(I)**2)
      AL(I)=TI*AL(I)
      BE(I)=TI*BE(I)
   40 CONTINUE
!
! CALCULATE SPLINE COEFFICIENTS
!
      DO 50 I=1,NM1
      D(I)=(AL(I)+BE(I)-2.)*DELTA(I)/H(I)**2
      C(I)=(3.-2.*AL(I)-BE(I))*DELTA(I)/H(I)
      B(I)=AL(I)*DELTA(I)
      IF(B(I)*DELTA(I).GE.0.) GOTO 50
!        write(*,*)b(i),delta(i)
      B(I)=DELTA(I)
      C(I)=0.
      D(I)=0.
   50 CONTINUE
      B(N)=BE(N-1)*DELTA(N-1)
      C(N)=0.
      D(N)=0.
      RETURN
!
  100 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.
      D(1) = 0.
      B(2) = B(1)
      C(2) = 0.
      D(2) = 0.
      RETURN
END SUBROUTINE

!------------------------------------------------------------------

REAL(8) FUNCTION SEVAL(N, U, X, Y, B, C, D)
   INTEGER :: N
   REAL(8) :: Y(N)
   DOUBLE PRECISION :: U, X(N), B(N), C(N), D(N)
   !
   !  THIS SUBROUTINE EVALUATES THE CUBIC SPLINE FUNCTION
   !
   !    SEVAL = Y(I) + B(I)*(U-X(I)) + C(I)*(U-X(I))**2 + D(I)*(U-X(I))**3
   !
   !    WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE
   !
   !  IF  U .LT. X(1) THEN  I = 1  IS USED.
   !  IF  U .GE. X(N) THEN  I = N  IS USED.
   !
   !  INPUT..
   !
   !    N = THE NUMBER OF DATA POINTS
   !    U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
   !    X,Y = THE ARRAYS OF DATA ABSCISSAS AND ORDINATES
   !    B,C,D = ARRAYS OF SPLINE COEFFICIENTS COMPUTED BY SPLINE
   !
   !  IF  U  IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A
   !  BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.
   !
   INTEGER I, J, K
   REAL DX
   DATA I/1/
   IF ( I .GE. N ) I = 1
   IF ( U .LT. X(I) ) GO TO 10
   IF ( U .LE. X(I+1) ) GO TO 30

   !  BINARY SEARCH
   !
   10 I = 1
      J = N+1
   20 K = (I+J)/2
      IF ( U .LT. X(K) ) J = K
      IF ( U .GE. X(K) ) I = K
      IF ( J .GT. I+1 ) GO TO 20
   !
   !  EVALUATE SPLINE
   !
   30 DX = U - X(I)
      SEVAL = Y(I) + DX*(B(I) + DX*(C(I) + DX*D(I)))
      RETURN
END FUNCTION SEVAL

!--------------------------------------------------------------------------
SUBROUTINE SPLMI4 (N, X, Y, B, C, D)
   INTEGER :: N
   REAL(4) :: Y(N)
   DOUBLE PRECISION :: X(N), B(N), C(N), D(N),DELTA(N),H(N),AL(N),BE(N)
   !
   !  THE COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N ARE COMPUTED
   !  FOR A MONOTONICALLY VARYING CUBIC INTERPOLATING SPLINE
   !
   !    S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
   !
   !    FOR  X(I) .LE. X .LE. X(I+1)
   !    WITH Y(I+1).GE.Y(I) (ALL I) OR Y(I+1).LT.Y(I) (ALL I)
   !
   !  INPUT..
   !
   !    N = THE NUMBER OF DATA POINTS OR KNOTS (N.GE.2)
   !    X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
   !    Y = THE ORDINATES OF THE KNOTS
   !
   !  OUTPUT..
   !
   !    B, C, D  = ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.
   !
   !*************************************************************************
   !
   ! THEORY FROM 'MONOTONE PIECEWISE CUBIC INTERPOLATION',
   ! BY F.N. FRITSCH AND R.E. CARLSON IN SIAM J.NUMER.ANAL.,V.17,P.238
   !
   !*************************************************************************
   !
   !    Y(I) = S(X(I))
   !    B(I) = SP(X(I))
   !    C(I) = SPP(X(I))/2
   !    D(I) = SPPP(X(I))/6  (DERIVATIVE FROM THE RIGHT)
   !
   !  THE ACCOMPANYING FUNCTION SUBPROGRAM  SEVAL  CAN BE USED
   !  TO EVALUATE THE SPLINE.
   !
   !
      NM1 = N-1
      IF ( N .LT. 2 ) RETURN
      IF ( N .LT. 3 ) GO TO 100
!
! CALCULATE THE H(I) AND DELTA(I)
!
      DO 10 I=1,NM1
      AL(I)=0.
      BE(I)=0.
      H(I)=X(I+1)-X(I)
   10 DELTA(I)=(Y(I+1)-Y(I))/H(I)
!
! CALCULATE FIRST VALUES FOR AL AND BE BY 3-POINT DIFFERENCE
!
      IF(DELTA(1).EQ.0) GOTO 15
      AL(1)=((H(1)+H(2))**2*Y(2)-H(1)**2*Y(3)-H(2)*(2.*H(1)+H(2))&
             *Y(1))/(H(2)*(H(1)+H(2))*(Y(2)-Y(1)))
   15 DO 20 I=2,NM1
      IF(DELTA(I).EQ.0) GOTO 20
      AL(I)=(H(I-1)**2*Y(I+1)+(H(I)**2-H(I-1)**2)*Y(I)-H(I)**2*&
             Y(I-1))/(H(I-1)*(H(I)+H(I-1))*(Y(I+1)-Y(I)))
   20 CONTINUE
!
      NM2=N-2
      DO 30 I=1,NM2
      IF(DELTA(I).EQ.0.) GOTO 30
      BE(I)=(H(I)**2*Y(I+2)+(H(I+1)**2-H(I)**2)*Y(I+1)-H(I+1)**2*&
             Y(I))/(H(I+1)*(H(I)+H(I+1))*(Y(I+1)-Y(I)))
   30 CONTINUE
!
      IF(DELTA(N-1).EQ.0.) GOTO 35
      BE(N-1)=(H(N-2)*(2.*H(N-1)+H(N-2))*Y(N)-(H(N-1)+H(N-2))**2&
             *Y(N-1)+H(N-1)**2*Y(N-2))/(H(N-2)*(H(N-1)+H(N-2))&
             *(Y(N)-Y(N-1)))
!
! CORRECT VALUES FOR AL AND BE
!
   35 DO 40 I=1,NM1
      IF(AL(I)+BE(I).LE.2.) GOTO 40
      IF(2.*AL(I)+BE(I).LE.3.) GOTO 40
      IF(AL(I)+2.*BE(I).LE.3.) GOTO 40
      PHI=AL(I)-(2.*AL(I)+BE(I)-3.)**2/(AL(I)+BE(I)-2.)/3.
      IF(PHI.GE.0.) GOTO 40
      TI=3./SQRT(AL(I)**2+BE(I)**2)
      AL(I)=TI*AL(I)
      BE(I)=TI*BE(I)
   40 CONTINUE
!
! CALCULATE SPLINE COEFFICIENTS
!
      DO 50 I=1,NM1
      D(I)=(AL(I)+BE(I)-2.)*DELTA(I)/H(I)**2
      C(I)=(3.-2.*AL(I)-BE(I))*DELTA(I)/H(I)
      B(I)=AL(I)*DELTA(I)
      IF(B(I)*DELTA(I).GE.0.) GOTO 50
!        write(*,*)b(i),delta(i)
      B(I)=DELTA(I)
      C(I)=0.
      D(I)=0.
   50 CONTINUE
      B(N)=BE(N-1)*DELTA(N-1)
      C(N)=0.
      D(N)=0.
      RETURN
!
  100 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.
      D(1) = 0.
      B(2) = B(1)
      C(2) = 0.
      D(2) = 0.
      RETURN
END SUBROUTINE

!------------------------------------------------------------------

REAL(8) FUNCTION SEVAL4(N, U, X, Y, B, C, D)
   INTEGER :: N
   REAL(4) :: Y(N)
   DOUBLE PRECISION :: U, X(N), B(N), C(N), D(N)
   !
   !  THIS SUBROUTINE EVALUATES THE CUBIC SPLINE FUNCTION
   !
   !    SEVAL = Y(I) + B(I)*(U-X(I)) + C(I)*(U-X(I))**2 + D(I)*(U-X(I))**3
   !
   !    WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE
   !
   !  IF  U .LT. X(1) THEN  I = 1  IS USED.
   !  IF  U .GE. X(N) THEN  I = N  IS USED.
   !
   !  INPUT..
   !
   !    N = THE NUMBER OF DATA POINTS
   !    U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
   !    X,Y = THE ARRAYS OF DATA ABSCISSAS AND ORDINATES
   !    B,C,D = ARRAYS OF SPLINE COEFFICIENTS COMPUTED BY SPLINE
   !
   !  IF  U  IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A
   !  BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.
   !
   INTEGER I, J, K
   REAL DX
   DATA I/1/
   IF ( I .GE. N ) I = 1
   IF ( U .LT. X(I) ) GO TO 10
   IF ( U .LE. X(I+1) ) GO TO 30

   !  BINARY SEARCH
   !
   10 I = 1
      J = N+1
   20 K = (I+J)/2
      IF ( U .LT. X(K) ) J = K
      IF ( U .GE. X(K) ) I = K
      IF ( J .GT. I+1 ) GO TO 20
   !
   !  EVALUATE SPLINE
   !
   30 DX = U - X(I)
      SEVAL4 = Y(I) + DX*(B(I) + DX*(C(I) + DX*D(I)))
      RETURN
END FUNCTION SEVAL4

!------------------------------------------------------
FUNCTION reallocate(p, n)               ! reallocate REAL
   REAL(8), POINTER, DIMENSION(:) :: p, reallocate
   INTEGER, intent(in) :: n
   INTEGER :: nold, ierr
   ALLOCATE(reallocate(1:n), STAT=ierr)
   IF(ierr /= 0) STOP "allocate error"
   IF(.NOT. ASSOCIATED(p)) RETURN
   nold = MIN(SIZE(p), n)
   reallocate(1:nold) = p(1:nold)
   reallocate(nold+1:) = 0.d0
   DEALLOCATE(p)
END FUNCTION reallocate
!------------------------------------------------------
FUNCTION reallocate2(p,m,n)               ! reallocate REAL
   REAL(8), POINTER, DIMENSION(:,:) :: p, reallocate2
   INTEGER, intent(in) :: m,n
   INTEGER :: nold, ierr
   ALLOCATE(reallocate2(m,n), STAT=ierr)
   IF(ierr /= 0) STOP "allocate error"
   IF(.NOT. ASSOCIATED(p)) RETURN
   nold = MIN(SIZE(p,2), n)
   reallocate2(:,1:nold) = p(:,1:nold)
   reallocate2(:,nold+1:) = 0.d0
   DEALLOCATE(p)
END FUNCTION reallocate2
!------------------------------------------------------
  REAL(8) FUNCTION BBFN(Q)
!!$     ********************************************************************
!!$     *  This subroutine calculates the fractional blackbody             *
!!$     *  emissive power f(n*lambda*T), where X=n*lambda*T in (micro-m*K) *
!!$     ********************************************************************
    REAL(8) :: PI,CC,C2,EPS,V,EX,M,EM,VM,BM,Q
    REAL(8),PARAMETER :: x1= 21.d0, x2= 75671.d0

    IF (Q <x1) THEN        ! if Q is too small, return 0.0
       BBFN=0.d0; RETURN
    ELSEIF (Q>x2) THEN     ! if Q is too large, return 1.0
       BBFN= 1.d0; RETURN
    ENDIF

    PI=4.D0*DATAN(1.D0)
    CC=1.5D1/PI**4
    C2=1.4388D4
    EPS=1.D-16

    V=C2/Q
    EX=DEXP(V)

    M=0
    BBFN=0.D0
    EM=1.D0
5   M=M+1
    VM=M*V
    BM=(6.D0+VM*(6.D0+VM*(3.D0+VM)))/M**4
    EM=EM/EX
    BBFN=BBFN+BM*EM
    IF(VM**3*EM.GT.EPS) GOTO 5
    BBFN=CC*BBFN
    RETURN
  END FUNCTION BBFN


!************************************************************************
  DOUBLE COMPLEX FUNCTION W4(Z)
! COMPUTES THE COMPLEX PROBABILITY FUNCTION W(Z)=EXP(-Z**2)*ERFC(-I*Z)
! IN THE UPPER HALF-PLANE Z=X+I*Y (I.E. FOR Y>=0)
! MAXIMUM RELATIVE ERROR OF BOTH REAL AND IMAGINARY PARTS IS <1*10**(-4)
    IMPLICIT NONE
      DOUBLE COMPLEX :: Z,T,U
      REAL(8) :: X,Y,S
      X=REAL(Z)
      Y=AIMAG(Z)
      T=CMPLX(Y,-X)
      S=ABS(X)+Y
      IF(S.GE.15.) THEN       !  ***  REGION I
        W4=T*.5641896/(.5+T*T)
      ELSEIF (S.GE.5.5) THEN  !  ***  REGION II
        U=T*T
        W4=T*(1.410474+U*.5641896)/(.75+U*(3.+U))
      ELSEIF (Y.GE..195*ABS(X)-.176) THEN  !  ***  REGION III
        W4=(16.4955+T*(20.20933+T*(11.96482+T*(3.778987+T*.5642236))))/ &
           (16.4955+T*(38.82363+T*(39.27121+T*(21.69274+T*(6.699398+T)))))
      ELSE                    !  ***  REGION IV
        U=T*T
        W4=CDEXP(U)-T*(36183.31-U*(3321.9905-U*(1540.787-U*(219.0313 &
           -U*(35.76683-U*(1.320522-U*.56419)))))) &
           /(32066.6-U*(24322.84-U*(9022.228-U*(2186.181-U*(364.2191 &
           -U*(61.57037-U*(1.841439-U)))))))
      ENDIF
      RETURN
    END FUNCTION W4

    SUBROUTINE VOIGT(S,BL,BD,DETA,KETA)
! COMPUTES ABSORPTION COEFFICIENT FOR A VOIGT LINE USING THE SUBROUTINE
! BY J. HUMLICEK, J. QUANT. SPECTR. RAD. TRANSFER, V.27, PP.437-444, 1982
!
! INPUT:
!   S     LINE INTENSITY (CM-2)
!   BL    LORENTZ HALFWIDTH AT HALFHEIGHT (CM-1)
!   BD    DOPPLER HALFWIDTH AT HALFHEIGHT (CM-1)
!   DETA  WAVENUMBERS AWAY FROM LINE CENTER (CM-1)
! OUTPUT:
!   KETA  SPECTRAL ABSORPTION COEFFICIENT AT A WAVENUMBER DETA REMOVED 
!         FROM LINE CENTER (CM-1)
!
    IMPLICIT NONE
      REAL(8) :: S, BL,BD,DETA,KETA,BD2
      DOUBLE COMPLEX :: Z
      REAL(8),PARAMETER:: RTLN2=0.8325546               ! SQRT(ALOG(2))
      REAL(8),PARAMETER:: RTPI=1.772454                 ! SQRT(PI)
      IF(BD.LT.1.E-10) BD=1.E-10    ! Avoid division by zero
! BD2 IS THE DOPPLER HALFWIDTH AT 1/e HEIGHT (CM-1)
      BD2=BD/RTLN2                  ! CORRECTED BD
      Z=CMPLX(DETA/BD2,BL/BD2)
      KETA=S/(RTPI*BD2)*REAL(W4(Z),8)
    END SUBROUTINE VOIGT

!-----------------------------------------------------------------------
 REAL(8) FUNCTION E1(X)
! called by EN(n,x) 
  IMPLICIT NONE
  REAL(8),INTENT(in) :: X
  REAL(8) :: A0,A1,A2,A3,A4,A5,B1,B2,B3,B4
  IF(X .LE. 1.0d0 ) THEN
     A0=-.57721566;A1= .99999193;A2=-.24991055
     A3= .05519968;A4=-.00976004;A5= .00107857
     E1=A0+X*(A1+X*(A2+X*(A3+X*(A4+X*A5))))-DLOG(X+1.D-8)
  ELSE
     A1= 8.5733287401;A2=18.0590169730;A3= 8.6347608925;A4=  .2677737343
     B1= 9.5733223454;B2=25.6329561486;B3=21.0996530827;B4= 3.9584969228
     E1=(X*(A3+X*(A2+X*(A1+X)))+A4)/(X*(B3+X*(B2+X*(B1+X)))+B4)*DEXP(-X)/X
  ENDIF
  RETURN
 END FUNCTION E1
!-----------------------------------------------------------------------
 REAL(8) FUNCTION EN(N,X)
! exponential integral of order n, eq. 13.31 in Modest book.
  IMPLICIT NONE
  INTEGER,INTENT(in) :: N
  REAL(8),INTENT(in) :: X
  INTEGER :: i
  REAL(8) :: EX
  EN=E1(X)
  IF(N .GT. 1 ) THEN
     EX=DEXP(-X)
     DO I=2,N
        EN=(EX-X*EN)/(I-1.D0)
     END DO
  ENDIF
  RETURN
 END FUNCTION EN
!-----------------------------------------------------------------------
 REAL(8) FUNCTION eb(T,wn)
! T in K and wn is wavenumber in 1/cm
! eb is the spectral blackbody emissive power (eq. 1.14 of Modest book)
! eb/pi is the spectral blackbody intensity, or Planck function

 IMPLICIT NONE
 REAL(8),INTENT(in) :: T, wn
 REAL(8),PARAMETER :: C1=3.7419d-12,C2=1.4388d0
 eb=C1*wn**3/(EXP(C2/T*wn)-1.d0)
 RETURN
 END FUNCTION eb  
!-----------------------------------------------------------------------
SUBROUTINE subtract(x, fun1, fun2, diff, lsquare)
! This subroutine calculates the absolute and relative errors
! between two functions and outputs it all to a file that is 
! also an input.
IMPLICIT NONE

!Declare constants and variables
REAL(8), INTENT(IN) :: x(:)		! Abscissa array for both fun1 and fun2
REAL(8), INTENT(IN) :: fun1(:)	! Ordinates to be compared with fun2
REAL(8), INTENT(IN) :: fun2(:)	! Ordinates to be compared with fun1
REAL(8), INTENT(OUT) :: diff(:)	! Ordinates to be compared with fun1
REAL(8), INTENT(OUT) :: lsquare	! Ordinates to be compared with fun1
INTEGER :: i, j, ierr	! Loop index and I/O variable

! If the arrays are not the same size, we cannot properly subtract them.
IF(SIZE(fun1) /= SIZE(fun2)) THEN
	WRITE(*,*) 'Error: Size mismatch. SIZE(fun1) /= SIZE(fun2)!'
	RETURN
ENDIF

lsquare = 0.d0
! Calculate errors and write to file
j = 1
DO i = 1,SIZE(fun1)
	IF(x(i) > 2070.d0 .AND. x(i) < 2200.d0) THEN
!	IF((x(i) > 3585.d0 .AND. x(i) < 3610.d0) .OR. (x(i) > 3617.d0 .AND. x(i) < 3635.d0) &
!		.OR. (x(i) > 3690.d0 .AND. x(i) < 3712.d0) .OR. (x(i) > 3717.d0 .AND. x(i) < 3740.d0)) THEN
		diff(i) = fun2(i) - fun1(i)
		j = j+1
	ELSE
		diff(i) = 0.d0
	ENDIF
	lsquare = lsquare + diff(i)**2.d0
END DO
lsquare = SQRT(lsquare/REAL(j-1))
END SUBROUTINE subtract
!-----------------------------------------------------------------------------
FUNCTION gasdev(idum)
IMPLICIT NONE
INTEGER :: idum
REAL(8) :: gasdev
INTEGER :: iset
REAL(8) :: fac,gset,rsq,v1,v2
SAVE    :: iset,gset
DATA iset/0/
IF(iset==0) THEN
1   v1 = 2.*ran1(idum)-1.
    v2 = 2.*ran1(idum)-1.
    rsq = v1**2+v2**2
    IF(rsq>=1.d0 .OR. rsq==0.d0) GOTO 1
    fac = SQRT(-2.*LOG(rsq)/rsq)
    gset = v1*fac
    gasdev = v2*fac
    iset = 1
ELSE
    gasdev = gset
    iset = 0
ENDIF
RETURN
CONTAINS
	FUNCTION ran1(idum)
	IMPLICIT NONE
	INTEGER :: idum,NDIV
	REAL(8) :: ran1,AM,EPS,RNMX
	INTEGER, PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32
	INTEGER :: j,k,iv(NTAB),iy
	SAVE iv,iy
	DATA iv /NTAB*0/, iy /0/
	
	AM = 1./IM;NDIV=1+(IM-1)/NTAB;EPS=1.2e-7;RNMX=1.-EPS
	IF(idum<=0 .OR. iy==0) THEN
	    idum = max(-idum,1)
	    DO 11 j = NTAB+8,1,-1
	        k = idum/IQ
	        idum = IA*(idum-k*IQ)-IR*k
	        IF (idum<0) idum=idum+IM
	        IF (j<=NTAB) iv(j)=idum
	11  CONTINUE
	    iy = iv(1)
	ENDIF
	k = idum/IQ
	idum = IA*(idum-k*IQ)-IR*k
	IF (idum<0) idum=idum+IM
	j = 1+iy/NDIV
	iy = iv(j)
	iv(j) = idum
	ran1 = min(AM*iy,RNMX)
	RETURN
	END FUNCTION
END FUNCTION

END MODULE commonData

