!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program reads in tabulated absorption coefficients and interpolates among
! them to calculate the absorption coefficient at a desired temperature (T) and
! concentratrion (x). The total pressure is not adjustable in this code; rather,
! the absorption coefficient database for a given total pressure should be
! created prior to using this code so that only interpolations in T and x are
! needed.

! After interpolating the absorption coefficient, this program calculates the
! emitted and/or transmitted normal intensity from a 1-D medium. It then
! convolves the intensity with the FTIR Instrument Line Function (ILF) at a 
! desired "resolution".

!  REVISIONS:
!
!     Date      Programmer      Changes
!     ====      ==========      =======
!	13Dec2010    T.A. Reeder     Updated the program to read newest version of
!                               HITEMP (2010).
!
!	09Mar2011	T.A. Reeder		  Implemented FFTW into code to calculate 
!                               convolution of ILF with calculated intensity.
!
!	21Mar2011    T.A. Reeder	  Added bilinear interpolation for calculating
!                               absorption coefficient from tabled data.
!
!	04Jun2011    T.A. Reeder     Implemented Anquan's unformatted absco DB and
!                               spline routines.
!
!	27Jun2011	 T.A. Reeder	  Made this code into a module to be called by
!										  inverse calculations.


MODULE forward_calc
USE commonData
INCLUDE '/public/home/lcc-rt/opt/fftw/api/fftw3.f'

! Other parameters
REAL(4),PUBLIC,ALLOCATABLE :: absc2(:)
REAL(4),PUBLIC,ALLOCATABLE :: absc_co2(:,:,:)
REAL(4),PUBLIC,ALLOCATABLE :: absc_h2o(:,:,:)
REAL(4),PUBLIC,ALLOCATABLE :: absc_co(:,:,:)
REAL(8),PRIVATE,PARAMETER :: wvnm_b = 200.d0	
REAL(8),PRIVATE,PARAMETER :: wvnm_e = 15000.d0
REAL(8),PRIVATE,PARAMETER :: wvnm_b_real = 5600.d0	
REAL(8),PRIVATE,PARAMETER :: wvnm_e_real = 6600.d0
INTEGER,PRIVATE :: number


CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE readDB(P,nP,nGb,nGe)
IMPLICIT NONE

! Passed variables
INTEGER,INTENT(IN) :: nP						! Number of DB Pressure
REAL(8),INTENT(IN) :: P(nP)						! DB pressure
INTEGER,INTENT(IN) :: nGb  			      		! beginning gas integer (1 = CO2, 2 = H2O, 3 = CO)
INTEGER,INTENT(IN) :: nGe						! ending gas integer

! FFTW parameters
INTEGER,PARAMETER :: nT = 8						! number of temperatures
REAL(8) :: Temp(nT)								! DB Temperatures
REAL(8) :: x = 0.0000							! DB mole fractions
INTEGER :: iT,iP,iG								! loop variables
REAL(8) :: Tmin=50.d0							! lower limit of temperature range
REAL(8) :: Tmax=400.d0							! upper limit of temperature range
INTEGER,PARAMETER :: habsc = 10             	! File handles
CHARACTER(80) :: abscDB                     	! File names
INTEGER :: dummy,ierr,nTdummy,number_dummy  
REAL(4) :: dummy1_r4, dummy2_r4

! Other parameters
TYPE(GasInfo) :: Gas_Info						! Gas information(gas, P, T, x)
REAL(8) :: t_start, t_finish
INTEGER	:: i
INTEGER	:: k
REAL(8),DIMENSION(:),ALLOCATABLE :: wvnm
!!!!!!!!!!!!!!!!!!!!!!!!!! End variable declaration !!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*) 'You are reading the DB for x = ', x

! Calculate temperature vector and determine bounds on T
DO iT = 1,SIZE(Temp)
   Temp(iT)=Tmin + (Tmax-Tmin)*(iT-1)/(nT-1)
END DO

DO iG = nGb,nGe		! This loop is for gas species
   IF(iG == 1) THEN
      Gas_Info = GasInfo(CO2,P(1),Temp(1),x)
      ALLOCATE(absc_co2(int((wvnm_e_real - wvnm_b_real)*10000.d0+100.d0),nP,SIZE(Temp))) ! = int((wvnm_e_real - wvnm_b_real)*10000.d0+100.d0)
      WRITE(*,*) 'Reading absco for CO2'
   ELSEIF(iG == 2) THEN
      Gas_Info = GasInfo(H2O,P(1),Temp(1),x)
      ALLOCATE(absc_h2o(7400001,nP,SIZE(Temp)))
      WRITE(*,*) 'Reading absco for H2O'
   ELSEIF(iG == 3) THEN
      Gas_Info = GasInfo(CO,P(1),Temp(1),x)
   	  ALLOCATE(absc_co(7400001,nP,SIZE(Temp)))
      WRITE(*,*) 'Reading absco for CO'
   ENDIF

   Gas_Info%x = x       ! Define Gas_Info pressure

   DO iP = 1,nP			! This loop is for concentration
      ! Write filepaths
      IF(Gas_Info%Gas==CO2) THEN
         WRITE(abscDB,101) nT,P(iP),x
			Gas_Info%P = P(iP)
      ELSEIF(Gas_Info%Gas==H2O) THEN
         WRITE(abscDB,201) nT,P(iP),x
			Gas_Info%P = P(iP)
      ELSEIF(Gas_Info%Gas==CO) THEN
         WRITE(abscDB,301) nT,P(iP),x
			Gas_Info%P = P(iP)
      ENDIF
      IF (abscDB(58:58)==' ') abscDB(58:58)='0'

!/public/home/lcc-rt/opt/abscDB_at/abscCO2.08.P.
!/public/home/lcc-rt04/abscDB_12/abscDB_at_9/abscCO2.08.P.
		! FORMAT statements
		100 FORMAT (1X, A4, F8.2, A12, I15)
		!101 FORMAT('/opt/abscDB/abscCO2.',I2.2,'.P.',F4.1,'m.',F4.2,'.sp.vt.dat')
		101 FORMAT('/public/home/lcc-rt04/DB/abscDB_at_9/abscCO2.',I2.2,'.P.',F5.3,'m.',F6.4,'.sp.vt.dat')
		201 FORMAT('/public/home/lcc-rt/opt/abscDB_at/abscH2O.',I2.2,'.P.',F4.2,'m.',F4.2,'.sp.vt.dat')
		301 FORMAT('/public/home/lcc-rt/opt/abscDB_at/abscCO_.',I2.2,'.P.',F4.2,'m.',F4.2,'.sp.vt.dat')
		! Read in file header. The file header contains information about the file, such as P, x, wvnm_b,
		! wvnm_e, and SIZE(Temp)
		
		WRITE(*,110) 'Reading absorbtion coefficient for P = ', P(iP), 'bar'
		110 FORMAT (1X, A39, F5.3, A5)
		
      OPEN (habsc, FILE = abscDB, FORM = 'UNFORMATTED', STATUS = 'OLD', ACTION = 'READ', IOSTAT = ierr)
      IF (ierr == 0) THEN
         READ(habsc) dummy1_r4,dummy2_r4	!Pressure Concentration
         READ(habsc) dummy1_r4,dummy2_r4	!Wvnm_b Wvnm_e
         READ(habsc) nTdummy
      ELSE
         WRITE(*,*) 'Error: Trouble opening '//abscDB
      ENDIF
		
		! Read in data for each temperature. Since the wavenumber spacing is not uniform, 'wvnm_resolution'
		! must be called each time through the loop. The header information for each temperature includes
		! Temp(iT), 'wvnmst', and 'number'. If desired, you can print these values.
      DO iT = 1, SIZE(Temp)
		Gas_Info%T = Temp(iT)	! Define Gas_Info temperature
		CALL wvnm_resolution(Gas_Info,.TRUE.,wvnm_b,wvnm_e,wvnmst,number)			! Determine wvnm stepsize
		ALLOCATE(absc2(number))
		WRITE(*,100)'T = ', Temp(iT), 'number = ', number
		! Read absco data from file
		READ(habsc) dummy1_r4, dummy2_r4	!Temperature Wvnmst
		READ(habsc) number_dummy			!Number
		!write(*,*) dummy1_r4, dummy2_r4
		!write(*,*) number_dummy
		!write(*,*) size(absc2),number
		READ(habsc) (absc2(i), i = 1,number)			!Read absorbtion coefficient

OPEN(UNIT = 1,FILE='testdata/DBNAN.dat',STATUS ='REPLACE',ACTION='WRITE', IOSTAT = ierr)
		
		k = 1
		DO i = 1,number

IF(isnan(absc2(i))) THEN
	WRITE(1,*) i*wvnmst
ENDIF

			!absc2(i) = absc2(i)*x(ix)*P*0.976d0		!Convert from pressure-based to linear.
			!absc2(i) = absc2(i)*P(iP)						!Convert from pressure-based to linear.    P/1.01325
				IF(Gas_Info%Gas == CO2) THEN
					
					IF (((i * wvnmst + wvnm_b) >= wvnm_b_real) .AND. ((i * wvnmst + wvnm_b) < wvnm_e_real)) THEN
						absc_co2(k,iP,iT) = absc2(i) 		! Assign to appropriate array. These values are now stored
						k = k + 1
					ENDIF
				ELSEIF(Gas_Info%Gas == H2O) THEN		! in global variables that can be called from any
					absc_h2o(i,iP,iT) = absc2(i)		! function or subroutine, etc.
				ELSEIF(Gas_Info%Gas == CO) THEN		
					absc_co(i,iP,iT) = absc2(i)
				ENDIF
		END DO !i
		DEALLOCATE(absc2)
      ENDDO !iT
      CLOSE(habsc)
   END DO !iP
   
END DO !iG
CLOSE(1)
END SUBROUTINE readDB


SUBROUTINE interp(P,T, xCO2, xH2O, xCO, wave_array, wavesize, absc, dab_dx)
IMPLICIT NONE
!!! CHANGABLE PARAMETERS!! MAKE SURE THESE ARE CORRECT BEFORE PROCEEDING!!!
LOGICAL :: switchEmit = .TRUE.  ! TRUE=emission; FALSE=transmission
LOGICAL :: WriteDataSwitch = .FALSE. ! Write intensity data at end of code
LOGICAL :: wave_switch = .FALSE. ! wvnm stepsize. false = OFF(const); true = ON

! Filepath beginnings for saving emitted/transmitted intensities (optional).

!!!!!!!!!!!!!!!!!!!!!!!!! END CHANGEABLE PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!

! Passed variables
REAL(8), INTENT(IN) :: P						! Total pressure
REAL(8), INTENT(IN) :: T, xCO2, xH2O, xCO		! Temp and concentrations
INTEGER, INTENT(IN) :: wavesize					! size of wave_array and absc
REAL(8), INTENT(IN) :: wave_array(wavesize)		! array of wavenumbers at which there are data
REAL(4), INTENT(OUT) :: absc(wavesize)			! Interpolated absorption coefficient
REAL(8), INTENT(OUT) :: dab_dx(wavesize)		! derivative of absco wrt part press

! Absorption Coefficient Interpolation Parameters
INTEGER,PARAMETER :: nGas = 1				! number of gases to loop over
INTEGER,PARAMETER :: nT = 8					! number of temperatures
INTEGER,PARAMETER :: np = 22				! number of mole fractions
REAL(8) :: DB_P(np)=(/0.01d0,0.02d0,0.03d0,0.04d0,0.05d0,0.06d0,0.07d0,0.08d0,0.09d0,&
			0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.0d0,1.1d0,1.2d0,1.3d0/)	! DB pressures
!REAL(8) :: DB_P(np)=(/0.8d0,0.9d0/)
REAL(8) Temp(nT)							! DB Temperatures
REAL(8) :: x = 0.0000						! DB mole fraction
INTEGER :: iP,iT,ix,iG						! loop variables
INTEGER,PARAMETER :: N = 4					! # of spline knots
REAL(8) :: TT(N)                            ! spline ordinates
INTEGER :: NTT(N), nxP(2), diff,diffT		! spline ordinate indices
REAL(8) :: B(N), C(N), D(N)                 ! spline coefficients
REAL(8) :: pwg                              ! = xCO2,xH2O,xCO...for current gas
REAL(8) :: xP(2),xT(2)						! upper and lower limits for x
INTEGER :: r1,r2							! Bound indices for wavenumber interval of interest
REAL(8) :: Tmin=50.d0						! lower limit of temperature range
REAL(8) :: Tmax=400.d0						! upper limit of temperature range
REAL(4),DIMENSION(:,:,:),ALLOCATABLE :: absc_cur ! output of lin. interp.
REAL(4) :: absc_n(wavesize,4)				! Intermediate value of absco
REAL(8) :: dab_dx_T(wavesize,4)				! Intermediate value of dab_dx
! Other parameters
TYPE(GasInfo) :: Gas_Info					! Gas information(gas, T, P, x)
REAL(8) :: tstart, tfinish					! For timing
REAL(8) :: t_start, t_finish				! For timing
INTEGER :: i,j,ierr,mi						! Loop indices & I/O success
!REAL(8),DIMENSION(:),ALLOCATABLE :: wvnm	! wavenumber
REAL(8) :: wvnm_act(2)						! used for interpolating wvnm
!!!!!!!!!!!!!!!!!!!!!!!!!! End variable declaration !!!!!!!!!!!!!!!!!!!!!!!!!!!

!================== get absco for single gas / mixture ========================
! Determine wavenumber stepsize. Typically, wvnmst will be fixed at 0.01; if
! you want to use the complete absco database, turn wave_switch ON. NOTE: IF 
! wave_switch is ON, you will also need to adjust the filepaths to the DB.

! Open output data file. ENSURE THIS IS CORRECT...WILL BE REPLACED
!OPEN(UNIT = 16, FILE = 'inverseDiagnostica.dat', &
!                                STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
!IF(ierr /= 0) THEN
!    WRITE(*,*) 'Error! Could not open output file!'
!    STOP
!ENDIF


! IF(wave_switch) CALL wvnm_resolution(Gas_Info,wave_switch,wvnmst,number)
absc_n=0.d0									! Initializing absc_n
dab_dx_T=0.d0								! Initializing dab_dx_T
absc=0.d0									! Initializing absc
dab_dx=0.d0									! Initializing dab_dx
110 FORMAT (1X, A40, F6.0, A3, F6.0)
! Calculate temperature vector and determine bounds on T
DO iT = 1,SIZE(Temp)
   Temp(iT)=Tmin + (Tmax-Tmin)*(iT-1)/(nT-1)
END DO
!WRITE(*,*)'*******************'
! Determine the indices for spline interpolation.
IF(T <= 50.d0) THEN
	DO j = 1,4
		TT(j) = 50.d0
		NTT(j) = 1
	END DO
ELSEIF(T < 100) THEN
	DO j = 1,4
		TT(j) = Temp(j)
		NTT(j) = j
	END DO
ELSEIF(T < 350) THEN
	diff = T/50
	DO j = 1,4
		TT(j) = Temp(diff-2+j)
		NTT(j) = diff-2+j
	END DO
ELSEIF(T < 400) THEN
	diff = T/50
	DO j = 1,4
		TT(j) = Temp(diff-3+j)
		NTT(j) = diff-3+j
	END DO
ELSEIF(T >= 400.D0) THEN
	DO j = 1,4
		TT(j) = 400.d0
		NTT(j) = 8
	END DO
ELSE
	WRITE(*,*) 'Error! Problem determining bounds on T.'
	STOP
ENDIF
!WRITE(*,*)''
DO iG = 1,nGas
   IF(iG == 1) THEN
		IF(xCO2 == 0) CYCLE
		Gas_Info = GasInfo(CO2,P,Temp(1),xCO2)
!		pwg = P
   ELSEIF(iG == 2) THEN
		IF(xH2O == 0) CYCLE
		Gas_Info = GasInfo(H2O,P,Temp(1),xH2O)
!		pwg = P
   ELSEIF(iG == 3) THEN
		IF(xCO == 0) CYCLE
		Gas_Info = GasInfo(CO,P,Temp(1),xCO)
!		pwg = P
   ENDIF
  
! Determine bounds on concentration
	Gas_Info%x = x
!write(*,*)pwg
	IF(P < 0.1d0) THEN
		diff = INT(P/0.01)
	ELSEIF (0.1d0 <= P <= 1.3d0) THEN
		diff = INT(P/0.1) + 9
	ELSE
		diff = 20
	ENDIF
!	diff = 1
	IF((P < 0.01d0 ).OR. (P > 1.3d0)) THEN
		xP(1) = DB_P(diff+1)
		xP(2) = DB_P(diff+2)
		nxP(1) = diff+1
		nxP(2) = diff+2
	ELSE
		xP(1) = DB_P(diff)
		xP(2) = DB_P(diff+1)
		nxP(1) = diff
		nxP(2) = diff+1
	ENDIF

	diffT = (T-100.D0)/50
	IF((T > 50.d0).and.(T < 399.9D0)) THEN
		xT(1) = Temp(diffT+1)
   		xT(2) = Temp(diffT+2)
	ELSE
		xT(1) = Temp(diffT+2)
		xT(2) = Temp(diffT+3)
	ENDIF

   ALLOCATE(absc_cur(wavesize,SIZE(xP,1),SIZE(TT)))

	! Linear interpolation between wavenumber
	DO iT = 1,SIZE(TT)
		DO iP = 1,SIZE(xP,1)	
			Gas_Info%P = xP(iP)
			Gas_Info%T = TT(iT)
OPEN(UNIT = 2,FILE='testdata/DBwvnm.dat',STATUS ='REPLACE',ACTION='WRITE', IOSTAT = ierr)
WRITE(2,*)xP(iP),TT(iT)
			! Determine wavenumbers for wavenumber interpolation
			CALL wvnm_resolution(Gas_Info,.TRUE., wvnm_b, wvnm_e, wvnmst, number)
			CALL CPU_TIME(t_finish)

			DO i = 1, SIZE(wave_array,1)
				wvnm_act(1) = INT(wave_array(i)/wvnmst)*wvnmst
				wvnm_act(2) = INT((wave_array(i)+1.d-8 + wvnmst)/wvnmst)*wvnmst
				r1 = INT((wvnm_act(1) - wvnm_b_real)/wvnmst + 1 + 1.d-6)
				r2 = INT((wvnm_act(2) - wvnm_b_real)/wvnmst + 1 + 1.d-6)

IF(wave_array(i)>=199.9d0 .AND. wave_array(i)<=205.d0) THEN
WRITE(2,*) wave_array(i), wvnm_act(1),wvnm_act(2)
WRITE(2,*) wvnmst,r1,r2
WRITE(2,*) absc_co2(r1,nxP(iP),NTT(iT)),absc_co2(r2,nxP(iP),NTT(iT))
ENDIF

				IF(Gas_Info%Gas == CO2) THEN
					CALL LinTerp(absc_co2(r1,nxP(iP),NTT(iT)),&
		 			absc_co2(r2,nxP(iP),NTT(iT)),absc_cur(i,iP,iT),&
					wvnm_act(:), wave_array(i))
					absc_cur(i,iP,iT) = absc_cur(i,iP,iT)*xCO2*P
				ELSEIF(Gas_Info%Gas == H2O) THEN
	            	CALL LinTerp(absc_h2o(r1,nxP(iP),NTT(iT)),&
					absc_h2o(r2,nxP(iP),NTT(iT)), absc_cur(i,iP,iT),&
					wvnm_act(:), wave_array(i))
				ELSEIF(Gas_Info%Gas == CO) THEN
					CALL LinTerp(absc_co(r1,nxP(iP),NTT(iT)),&
					absc_co(r2,nxP(iP),NTT(iT)), absc_cur(i,iP,iT),&
					wvnm_act(:), wave_array(i))
				ENDIF
			END DO !i
!			DEALLOCATE(wvnm)
		END DO !iP
	END DO !iT
	
	! Linearly interpolate between pressure
	DO iT = 1,SIZE(TT)
		DO i = 1,SIZE(wave_array,1)
			CALL LinTerp(absc_cur(i,1,iT), absc_cur(i,2,iT), absc_n(i,iT),xP(:),P)
			dab_dx_T(i,iT) = (absc_cur(i,2,iT) - absc_n(i,iT))/(xP(2) - P)
!			absc_n(i,iT) = absc_n(i,iT)*pwg
!			dab_dx_T(i,iT) = dab_dx_T(i,iT)*pwg
			absc_n(i,iT) = absc_n(i,iT)
			dab_dx_T(i,iT) = dab_dx_T(i,iT)
		END DO
	END DO

	!
	DO i = 1,SIZE(absc)
		IF(T <= 50.d0 .OR. T >= 400.d0) THEN
			absc(i) = absc_n(i,1)
			CALL SPLMI(SIZE(TT), TT, dab_dx_T(i,:), B, C, D)
			dab_dx_T(i,2) = dab_dx_T(i,1)
		ELSE
			CALL SPLMI4(SIZE(TT), TT, absc_n(i,:), B, C, D)
			absc(i) = absc(i) + SEVAL4(SIZE(TT), T, TT, absc_n(i,:),B,C,D)
			IF(T <= 50) THEN
				dab_dx(i)=0.D0
			ELSEIF(T < 100) THEN
				dab_dx(i)=(absc_n(i,2)-absc(i))/(xT(2)-T)
			ELSEIF(T < 350) THEN
				dab_dx(i)=(absc_n(i,3)-absc(i))/(xT(2)-T)
			ELSEIF(T < 400) THEN
				dab_dx(i)=(absc_n(i,4)-absc(i))/(xT(2)-T)
			ELSEIF(T >= 400) THEN
				dab_dx(i)=0.D0
			ENDIF
		CALL SPLMI(SIZE(TT), TT, dab_dx_T(i,:), B, C, D)
		dab_dx_T(i,2) = dab_dx_T(i,2)+SEVAL(SIZE(TT),T,TT,dab_dx_T(i,:),B,C,D)
		ENDIF
	END DO
!	write(*,*) absc(1)
!	DEALLOCATE(absc_cur)
END DO ! nGas
END SUBROUTINE interp

SUBROUTINE soot_absc(wave_array, wavesize, absc)
	! Passed variables
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: wavesize				! size of wave_array and absc
	REAL(8), INTENT(IN) :: wave_array(wavesize) ! array of wavenumbers at which there are data
	REAL(4), INTENT(OUT) :: absc(wavesize) 		! Interpolated absorption coefficient		
	REAL(8) :: n(wavesize), k(wavesize), lamd(wavesize)	
	REAL(8),parameter		:: pi=3.14159265d0
	INTEGER:: i
	DO i=1,wavesize
	lamd(i)=10000.d0/wave_array(i)
	n(i)=1.811d0+0.1263d0*(log(lamd(i)))+0.027d0*(log(lamd(i)))**2+0.0417d0*(log(lamd(i)))**3
	k(i)=0.5821d0+0.1213d0*(log(lamd(i)))+0.2309d0*(log(lamd(i)))**2-0.01d0*(log(lamd(i)))**3
	absc(i)=36.d0*pi*n(i)*k(i)/((n(i)**2-k(i)**2+2.d0)**2+4.d0*n(i)**2*k(i)**2)/lamd(i)*10000.D0
	END DO  
END SUBROUTINE soot_absc

!------------------------------------------------------------------------------
!SUBROUTINE wvnm_resolution(Gas_Info,wave_switch,wvnmst,number)
!!!$ Based on the gas information, decide the wavenumber resolution
!!!$  Input:
!!!$    Gas_Info     --  information of the gas
!!!$  Output:
!!!$    wvnmst       --  wavenumber step for absco generation
!!!$    number       --  total number of wavenumber locations for absco generation 
!  IMPLICIT NONE
!    TYPE(GasInfo),INTENT(IN) :: Gas_Info
!    LOGICAL,INTENT(IN) :: wave_switch
!    REAL(8),INTENT(OUT)      :: wvnmst
!    INTEGER,INTENT(OUT)      :: number
!
!!!$ local variables
!    REAL(8),PARAMETER :: EPS= 1.d-6
!    REAL(8) :: P,x,T, w
!
!    P= Gas_Info%P; x= Gas_Info%x; T= Gas_Info%T
!
!IF(wave_switch) THEN
!    IF (P>5.0d0) THEN
!      w= 0.01d0
!
!    ELSEIF (P> 1.d0-EPS) THEN
!      IF (x<0.25d0) THEN
!        IF (T>1900.d0) THEN; w= 4.d-3
!        ELSEIF (T> 1100.d0) THEN; w= 5.d-3
!        ELSE; w= 1.d-2; ENDIF
!      ELSEIF (x<0.5d0) THEN
!        IF (T>2300.d0) THEN; w= 5.d-3
!        ELSE; w= 1.d-2; ENDIF
!      ELSE; w= 1.d-2; ENDIF
!
!    ELSEIF (P> 0.7d0-EPS) THEN
!      IF (x<0.25d0) THEN
!        IF (T>1500.d0) THEN; w= 2.d-3
!        ELSEIF (T>1300.d0) THEN; w= 4.d-3
!        ELSEIF (T>900.d0) THEN; w= 5.d-3
!        ELSE; w= 1.d-2; ENDIF
!      ELSEIF (x<0.5d0) THEN
!        IF (T> 1600.d0) THEN; w= 5.d-3
!        ELSE; w= 1.d-2; ENDIF
!      ELSE; w= 1.d-2; ENDIF
!    ENDIF
!ELSE
!	 w = 1.d-2
!ENDIF
!
!    wvnmst= w
!    number=(wvnm_e-wvnm_b)/wvnmst+EPS+1
!  END SUBROUTINE wvnm_resolution
!!$------------------------------------------------------------------------
!
!   SUBROUTINE LinTerp(Y1, Y2, Yout, XX, x)
!   ! This routine performs a linear interpolation
!   ! Input: 
!   !        Y1     !! wvnmst must be the same for Y1 and Y2 so they are of the same length.
!   !        Y2     !! Otherwise, you will be interpolating between abscos of different wavenumbers.
!   !        XX     !! array of length 2 containing concentrations to interpolate over.
!   !        x      !! Concentration to calc absco at.
!
!   ! Output:
!   !        Yout   !! The interpolated absco.
!
!      IMPLICIT NONE
!      REAL(4), INTENT(IN) :: Y1, Y2
!      REAL(4), INTENT(OUT) :: Yout
!      REAL(8), INTENT(IN) :: XX(2),x      ! X(1) = lower value; X(2) = upper value
!      REAL(8) :: diff1, diff2
!
!      diff1 = x - XX(1)
!      diff2 = XX(2) - XX(1)
!      Yout = Y1 + diff1/diff2*(Y2 - Y1)
!   END SUBROUTINE
!
!!!$------------------------------------------------------------------------
!
!   SUBROUTINE LinTerpDouble(Y1, Y2, Yout, XX, x)
!   ! This routine performs a linear interpolation
!   ! Input: 
!   !        Y1     !! wvnmst must be the same for Y1 and Y2 so they are of the same length.
!   !        Y2     !! Otherwise, you will be interpolating between abscos of different wavenumbers.
!   !        XX     !! array of length 2 containing concentrations to interpolate over.
!   !        x      !! Concentration to calc absco at.
!
!   ! Output:
!   !        Yout   !! The interpolated absco.
!
!      IMPLICIT NONE
!      REAL(8), INTENT(IN) :: Y1, Y2
!      REAL(8), INTENT(OUT) :: Yout
!      REAL(8), INTENT(IN) :: XX(2),x      ! X(1) = lower value; X(2) = upper value
!      REAL(8) :: diff1, diff2
!
!      diff1 = x - XX(1)
!      diff2 = XX(2) - XX(1)
!      Yout = Y1 + diff1/diff2*(Y2 - Y1)
!   END SUBROUTINE
!
!!--------------------------------------------------------------------------
!
!SUBROUTINE SPLMI (N, X, Y, B, C, D)
!   INTEGER :: N
!   REAL(4) :: Y(N)
!   DOUBLE PRECISION :: X(N), B(N), C(N), D(N),DELTA(N),H(N),AL(N),BE(N)
!   !
!   !  THE COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N ARE COMPUTED
!   !  FOR A MONOTONICALLY VARYING CUBIC INTERPOLATING SPLINE
!   !
!   !    S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
!   !
!   !    FOR  X(I) .LE. X .LE. X(I+1)
!   !    WITH Y(I+1).GE.Y(I) (ALL I) OR Y(I+1).LT.Y(I) (ALL I)
!   !
!   !  INPUT..
!   !
!   !    N = THE NUMBER OF DATA POINTS OR KNOTS (N.GE.2)
!   !    X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
!   !    Y = THE ORDINATES OF THE KNOTS
!   !
!   !  OUTPUT..
!   !
!   !    B, C, D  = ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.
!   !
!   !*************************************************************************
!   !
!   ! THEORY FROM 'MONOTONE PIECEWISE CUBIC INTERPOLATION',
!   ! BY F.N. FRITSCH AND R.E. CARLSON IN SIAM J.NUMER.ANAL.,V.17,P.238
!   !
!   !*************************************************************************
!   !
!   !    Y(I) = S(X(I))
!   !    B(I) = SP(X(I))
!   !    C(I) = SPP(X(I))/2
!   !    D(I) = SPPP(X(I))/6  (DERIVATIVE FROM THE RIGHT)
!   !
!   !  THE ACCOMPANYING FUNCTION SUBPROGRAM  SEVAL  CAN BE USED
!   !  TO EVALUATE THE SPLINE.
!   !
!   !
!      NM1 = N-1
!      IF ( N .LT. 2 ) RETURN
!      IF ( N .LT. 3 ) GO TO 100
!!
!! CALCULATE THE H(I) AND DELTA(I)
!!
!      DO 10 I=1,NM1
!      AL(I)=0.
!      BE(I)=0.
!      H(I)=X(I+1)-X(I)
!   10 DELTA(I)=(Y(I+1)-Y(I))/H(I)
!!
!! CALCULATE FIRST VALUES FOR AL AND BE BY 3-POINT DIFFERENCE
!!
!      IF(DELTA(1).EQ.0) GOTO 15
!      AL(1)=((H(1)+H(2))**2*Y(2)-H(1)**2*Y(3)-H(2)*(2.*H(1)+H(2))&
!             *Y(1))/(H(2)*(H(1)+H(2))*(Y(2)-Y(1)))
!   15 DO 20 I=2,NM1
!      IF(DELTA(I).EQ.0) GOTO 20
!      AL(I)=(H(I-1)**2*Y(I+1)+(H(I)**2-H(I-1)**2)*Y(I)-H(I)**2*&
!             Y(I-1))/(H(I-1)*(H(I)+H(I-1))*(Y(I+1)-Y(I)))
!   20 CONTINUE
!!
!      NM2=N-2
!      DO 30 I=1,NM2
!      IF(DELTA(I).EQ.0.) GOTO 30
!      BE(I)=(H(I)**2*Y(I+2)+(H(I+1)**2-H(I)**2)*Y(I+1)-H(I+1)**2*&
!             Y(I))/(H(I+1)*(H(I)+H(I+1))*(Y(I+1)-Y(I)))
!   30 CONTINUE
!!
!      IF(DELTA(N-1).EQ.0.) GOTO 35
!      BE(N-1)=(H(N-2)*(2.*H(N-1)+H(N-2))*Y(N)-(H(N-1)+H(N-2))**2&
!             *Y(N-1)+H(N-1)**2*Y(N-2))/(H(N-2)*(H(N-1)+H(N-2))&
!             *(Y(N)-Y(N-1)))
!!
!! CORRECT VALUES FOR AL AND BE
!!
!   35 DO 40 I=1,NM1
!      IF(AL(I)+BE(I).LE.2.) GOTO 40
!      IF(2.*AL(I)+BE(I).LE.3.) GOTO 40
!      IF(AL(I)+2.*BE(I).LE.3.) GOTO 40
!      PHI=AL(I)-(2.*AL(I)+BE(I)-3.)**2/(AL(I)+BE(I)-2.)/3.
!      IF(PHI.GE.0.) GOTO 40
!      TI=3./SQRT(AL(I)**2+BE(I)**2)
!      AL(I)=TI*AL(I)
!      BE(I)=TI*BE(I)
!   40 CONTINUE
!!
!! CALCULATE SPLINE COEFFICIENTS
!!
!      DO 50 I=1,NM1
!      D(I)=(AL(I)+BE(I)-2.)*DELTA(I)/H(I)**2
!      C(I)=(3.-2.*AL(I)-BE(I))*DELTA(I)/H(I)
!      B(I)=AL(I)*DELTA(I)
!      IF(B(I)*DELTA(I).GE.0.) GOTO 50
!!        write(*,*)b(i),delta(i)
!      B(I)=DELTA(I)
!      C(I)=0.
!      D(I)=0.
!   50 CONTINUE
!      B(N)=BE(N-1)*DELTA(N-1)
!      C(N)=0.
!      D(N)=0.
!      RETURN
!!
!  100 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
!      C(1) = 0.
!      D(1) = 0.
!      B(2) = B(1)
!      C(2) = 0.
!      D(2) = 0.
!      RETURN
!END SUBROUTINE
!
!!------------------------------------------------------------------
!
!REAL(8) FUNCTION SEVAL(N, U, X, Y, B, C, D)
!   INTEGER :: N
!   REAL(4) :: Y(N)
!   DOUBLE PRECISION :: U, X(N), B(N), C(N), D(N)
!   !
!   !  THIS SUBROUTINE EVALUATES THE CUBIC SPLINE FUNCTION
!   !
!   !    SEVAL = Y(I) + B(I)*(U-X(I)) + C(I)*(U-X(I))**2 + D(I)*(U-X(I))**3
!   !
!   !    WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE
!   !
!   !  IF  U .LT. X(1) THEN  I = 1  IS USED.
!   !  IF  U .GE. X(N) THEN  I = N  IS USED.
!   !
!   !  INPUT..
!   !
!   !    N = THE NUMBER OF DATA POINTS
!   !    U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
!   !    X,Y = THE ARRAYS OF DATA ABSCISSAS AND ORDINATES
!   !    B,C,D = ARRAYS OF SPLINE COEFFICIENTS COMPUTED BY SPLINE
!   !
!   !  IF  U  IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A
!   !  BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.
!   !
!   INTEGER I, J, K
!   REAL DX
!   DATA I/1/
!   IF ( I .GE. N ) I = 1
!   IF ( U .LT. X(I) ) GO TO 10
!   IF ( U .LE. X(I+1) ) GO TO 30
!
!   !  BINARY SEARCH
!   !
!   10 I = 1
!      J = N+1
!   20 K = (I+J)/2
!      IF ( U .LT. X(K) ) J = K
!      IF ( U .GE. X(K) ) I = K
!      IF ( J .GT. I+1 ) GO TO 20
!   !
!   !  EVALUATE SPLINE
!   !
!   30 DX = U - X(I)
!      SEVAL = Y(I) + DX*(B(I) + DX*(C(I) + DX*D(I)))
!      RETURN
!END FUNCTION SEVAL
!
!!------------------------------------------------------
!FUNCTION reallocate(p, n)               ! reallocate REAL
!	REAL(8), POINTER, DIMENSION(:) :: p, reallocate
!	INTEGER, intent(in) :: n
!	INTEGER :: nold, ierr
!	ALLOCATE(reallocate(1:n), STAT=ierr)
!	IF(ierr /= 0) STOP "allocate error"
!	IF(.NOT. ASSOCIATED(p)) RETURN
!	nold = MIN(SIZE(p), n)
!	reallocate(1:nold) = p(1:nold)
!	reallocate(nold+1:) = 0.d0
!	DEALLOCATE(p) 
!END FUNCTION reallocate
!!------------------------------------------------------
END MODULE forward_calc
