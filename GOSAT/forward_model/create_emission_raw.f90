	PROGRAM testsinc8
	USE forward_calc
	USE absc_kg_gen
	USE commonData
	IMPLICIT NONE
	! 	Date:			Author:			Notes:
	!	=====			=======			======
	!	16Nov2011	T.A. Reeder		This code calculates absorption coefficent either
	!								by reading and interpolating the absorption
	!								coefficent DB or by direct calculation. It then 
	!								calculates the transmitttance and uses the free
	!								software FFTW to convolve it with an FTIR
	!								instrument line function/shape (ILF or ILS). Both
	!								actual and ideal ILFs are used.
	!
	!								The relative and absolute errors between measured
	!								and calculated data are written to a file at the
	!								end of the program.
	!
	!								It is important to know that FFTW requires both
	!								input arrays to have the same length and the same
	!								spacing between data. If otherwise, you will
	!								either (a) get an error or (b) obtain an incorrect
	!								convolution.
	!-------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER,PARAMETER 				::	species=1												!	CO2==1,H2O==2,CO==3
	INTEGER,PARAMETER				::	nP = 22													!	Number of DB Pressure
	REAL(8),PARAMETER 				:: 	DB_P(np)=(/0.01d0,0.02d0,0.03d0,0.04d0,0.05d0,0.06d0,0.07d0,0.08d0,0.09d0,&
										0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.0d0,1.1d0,1.2d0,1.3d0/)	! 	DB Pressure [bar]
!	REAL(8),PARAMETER				::	DB_P(np)=(/0.8d0,0.9d0/)
	INTEGER,PARAMETER				::  nsb = 28												! number of parameters
	REAL(8),PARAMETER				::  L = 9000000.0d0											! optical path length [cm] L1=10.d0.l2=100.d0
	REAL(8)							::  P(nsb),T(nsb),lx(nsb),dx(nsb)
	REAL(8),ALLOCATABLE				::  seed(:,:), seed_ALB(:,:)
	INTEGER							:: 	n_profile
	REAL(8)							::  fvCO2(nsb),fvH2O(nsb),fvCO(nsb)
	REAL(8),PARAMETER 				:: 	pi=3.14159265d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------
	LOGICAL 						:: Pressure_switch = .FALSE.						! .FALSE. = DB; .TRUE. = calculated
	LOGICAL 						:: convmethod = .TRUE.								! Brute force (false) of FFTW(true)?
	REAL(8)							:: xCO2 											! CO2 concentration [-]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------

	REAL(8), PARAMETER  			:: ils_step_l = 0.00407029d0						! Step size to interpolate ils/absco 
	REAL(8)							:: ils_step											! Step size to interpolate ils/absco 

	!! you have to change 'wvnm_b_real' and 'wvnm_e_real' in 'forward_calc.f90' too
	REAL(8) 						:: wave_bb = 5781.13d0	!259.3941586793d0
	REAL(8) 						:: wave_ee = 6600.0d0						!7997.8258441756d0 
	REAL(8) :: wave_43_b = 6100.33d0
	REAL(8) :: wave_43_e = 6319.78d0
	
INTEGER							:: idum1
	REAL(8) 						:: resolution											! Return value ftir

	!-------------------------------------------------------------------------------------------
	REAL(8) 						:: ils_wave_ref											! used to properly space ils data
	REAL(8) 						:: wave_ref												! used to properly space tran_final
	REAL(8) 						:: wave_step											! wvnm(2) - wvnm(1)
	REAL(8) 						:: ils_interp(2)										! bordering values passed to LinTerp
	REAL(8) 						:: tran_interp(2)										! bordering values passed to LinTerp
	REAL(4), ALLOCATABLE 			:: absco(:,:),absco_in(:)				! interpolated and spaced abs.coeff.
	REAL(8), ALLOCATABLE 			:: dab_dx_in(:)											! interpolated and spaced abs.coeff.
	REAL(8), ALLOCATABLE 			:: wvnm(:)												! wavenumbers read from FTIR data
	REAL(8), ALLOCATABLE 			:: ils_wave(:)											! wavenumbers read from ils spectrum
	REAL(8), ALLOCATABLE 			:: absco_wave(:)										! wavenumbers calculated for absco
	REAL(8), ALLOCATABLE 			:: ils(:)												! ils values from ils spectrum
	REAL(8), POINTER 				:: ils_spaced(:)										! interpolated and spaced ils
	REAL(8), POINTER 				:: ils_spaced1(:)
	REAL(8), POINTER 				:: transmit(:)											! calculated transmittance
	REAL(8), ALLOCATABLE			:: rad(:) 
	REAL(8), ALLOCATABLE			:: label_rad(:) 
	REAL(8), ALLOCATABLE			:: labels_OT(:) 
	REAL(8), ALLOCATABLE			:: ALBEDO(:),ALBEDO_ILSF(:) 
	REAL(8), ALLOCATABLE	 		:: opticalThickness(:)									! calculated opticalThickness
	REAL(8), ALLOCATABLE 			:: tran_interim(:)										! Normalized tran with wings truncated
	REAL(8), ALLOCATABLE 			:: tran_final(:)										! Final transmittance...output to file
	REAL(8), POINTER 				:: tran_compare(:,:)									! Return value transmittance
	REAL(8), POINTER 				:: wvnm_compare(:,:)									! Return value wavenumbers
	INTEGER 						:: i, j, k, q, m, s,ne,ii,conv							! loop indices
	INTEGER 						:: sample,sample_b,sample_e,n_profile_b,n_profile_e	! loop indices
	INTEGER 						:: ierr													! I/O variable; 0=success
	INTEGER 						:: wavesize												! size of ftir data
	INTEGER 						:: ils_size												! size of ils file
	INTEGER 						:: PROFILE_size											! size of PROFILE file
	INTEGER 						:: ALBEDO_size											! size of ALBEDO file
	INTEGER 						:: absco_size											! size of initial absco
	CHARACTER(100)					:: suffix
	! FFTW parameters
	INTEGER 						:: convsize												! size of arrays into FFTW
	REAL(8), ALLOCATABLE 			:: tran_conv(:)											! Output from FFTW
	COMPLEX(8), ALLOCATABLE 		:: ils_ft(:)											! Fourier transform of ils
	COMPLEX(8), ALLOCATABLE 		:: tran_ft(:)											! FT of transmittance
	COMPLEX(8), ALLOCATABLE 		:: conv_in(:)											! conv_in(i)=ils_ft(i)*tran_ft(i)
	INTEGER(8) 									:: plan1, plan2								! FFTW plan integers
	! Other parameters and variables
	INTEGER							:: tcount,trate,tmax
	CHARACTER(100) 					:: dummy
	CHARACTER(150) 					:: outputfile 
	CHARACTER(150) 					:: ilsdatafile 
	REAL(8) 							:: T_low,slop,T_slop,xCO2_slop,theta_slop,fai_slop,ALB_slop
	REAL(8) 							:: xCO2_low
	REAL(8), external:: ran1
	INTEGER :: i_43_b,i_43_e,i_27_b,i_27_e

	REAL(8)::SS(162919,3)
	REAL(8)::SO(162919)

	REAL(8)	::Rayleigh(2000000,28),ALB(3083,26),ALB_ref(20,1)
	REAL(8)	::n(2000000)
	REAL(8) ::rn(2000000)
	REAL(8)	::Ps, Ts, Ns, theta, fai 

	TYPE(GasInfo) 							:: Gas_Info										! Derived data type; gas properties

	CALL SYSTEM_CLOCK(tcount,trate,tmax)
	write(*,*)tcount,trate,tmax

	idum1 = -tcount
	OPEN(UNIT=12, FILE= 'flame/ALBEDO_reference.dat', STATUS='OLD',ACTION='READ')
	IF(ierr /= 0) THEN
		WRITE(*,*) 'Error! Could not open output file!'
	ENDIF
	DO i = 1,20
		READ(12,*)(ALB_ref(i,1))
	END DO
	
	OPEN(UNIT=13, FILE= 'flame/ALBEDO_090728.dat', STATUS='OLD',ACTION='READ')	!ALBEDO_090701_120731 ALBEDO_090728	
	OPEN(UNIT=14, FILE= 'flame/PROFILE_090728.dat', STATUS='OLD',ACTION='READ') !PROFILE_090701_120731	PROFILE_090728

!!!!! Read PROFILE data to determine length of file
	i = 1
	DO
		READ(13,FMT = '(A28)',IOSTAT = ierr) dummy
		IF(ierr == -1) THEN ! Reached end of file, exit loop
		    ALBEDO_size = i-1
		    EXIT
		ENDIF
		i = i+1
	END DO
	ALLOCATE(seed_ALB(ALBEDO_size,26))
	Rewind(13)
	DO i = 1,ALBEDO_size
		READ(13,*)(seed_ALB(i,j),j=1,26)
	END DO
	i = 1
	DO
		READ(14,FMT = '(A28)',IOSTAT = ierr) dummy
		IF(ierr == -1) THEN ! Reached end of file, exit loop
		    PROFILE_size = i-1
		    EXIT
		ENDIF
		i = i+1
	END DO
	ALLOCATE(seed(PROFILE_size,28))
	Rewind(14)
	DO i = 1,PROFILE_size
		READ(14,*)(seed(i,j),j=1,28)
	END DO
	WRITE(*,*)'PROFILE_size',PROFILE_size
	WRITE(*,*)'ALBEDO_size',ALBEDO_size
	
	suffix = '090728'
	n_profile_b = 1
	n_profile_e = ALBEDO_size
	sample_b = 1
	sample_e = 1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	OPEN(UNIT = 18, FILE= 'testdata/labels_FIT_'//trim(suffix)//'.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
	IF(ierr /= 0) THEN
		WRITE(*,*) 'Error! Could not open output file!'
	ENDIF

	OPEN(UNIT = 19, FILE= 'testdata/labels_log_'//trim(suffix)//'.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
	IF(ierr /= 0) THEN
		WRITE(*,*) 'Error! Could not open output file!'
	ENDIF

	OPEN(UNIT=50, FILE= 'flame/irradabs1560.dat', STATUS='OLD',ACTION='READ')
	DO i=1,162919
		READ(50,*)(SS(i,j),j=1,3)
		SO(162919-i+1) = SS(i,1)*SS(i,2)/SS(i,3)/10000.d0
	END DO
	write(*,*) 'Solar Radiation sinfo : ',size(SO), SO(1),SO(2), SO(73706)

	OPEN(UNIT = 20, FILE= 'testdata/labels_log_OT_ILSF_'//trim(suffix)//'.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
	IF(ierr /= 0) THEN
		WRITE(*,*) 'Error! Could not open output file!'
	ENDIF

	OPEN(UNIT = 21, FILE= 'testdata/labels_T_'//trim(suffix)//'.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
	IF(ierr /= 0) THEN
		WRITE(*,*) 'Error! Could not open output file!'
	ENDIF

	OPEN(UNIT = 22, FILE= 'testdata/labels_x_'//trim(suffix)//'.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr) !labels_CO2
	IF(ierr /= 0) THEN
		WRITE(*,*) 'Error! Could not open output file!'
	ENDIF
	OPEN(UNIT = 24, FILE= 'testdata/OT_ILSF_'//trim(suffix)//'.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)  !labels_CO
	IF(ierr /= 0) THEN
		WRITE(*,*) 'Error! Could not open output file!'
	ENDIF
	OPEN(UNIT = 28, FILE= 'testdata/ALBEDO_ILSF_'//trim(suffix)//'.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
	IF(ierr /= 0) THEN
		WRITE(*,*) 'Error! Could not open output file!'
	ENDIF
	OPEN(UNIT = 29, FILE= 'testdata/Theta&Fai_'//trim(suffix)//'.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
	IF(ierr /= 0) THEN
		WRITE(*,*) 'Error! Could not open output file!'
	ENDIF
	
	ils_step = ils_step_l
	resolution=0.1995d0				!real(0.05d0*resolu)
	write(*,*)'reslution=',resolution
	write(*,*)'ils_step=',ils_step
	WRITE(ilsdatafile,100)
	100 FORMAT('./flame/ILSF_6230.dat')

	wavesize = INT((wave_ee - wave_bb)/resolution + 1 + 1.E-6)
	WRITE(*,*) 'waesize =' ,wavesize
	ALLOCATE(wvnm(wavesize))
	wvnm(1)=wave_bb
	DO i=1,wavesize-1
		wvnm(i+1) = wvnm(i) + resolution
	END DO

	wave_b = wvnm(1)
	wave_e = wvnm(wavesize)

	absco_size = INT((wvnm(SIZE(wvnm))-wvnm(1))/ils_step) + 1 + 1.E-6
	write(*,*)'absco_size =',absco_size
	ALLOCATE(absco_wave(absco_size))

	absco_wave(1) = wvnm(1)
	DO i=1,absco_size-1
		absco_wave(i+1) = absco_wave(i) + ils_step
	END DO

	OPEN(UNIT = 40, FILE = ilsdatafile, STATUS = 'OLD', ACTION = 'READ')	
		! Read ILF data to determine length of file
		i = 1
		DO
			READ(40,FMT = '(A28)',IOSTAT = ierr) dummy
			IF(ierr == -1) THEN ! Reached end of file, exit loop
			    ils_size = i-1
			    EXIT
			ENDIF
			i = i+1
		END DO

	IF( ALLOCATED(ils) )  DEALLOCATE(ils )
	write(*,*)'ils_size = ', ils_size
		! Allocate ils_wave, ils

		ALLOCATE(ils(ils_size))
		ALLOCATE(ils_wave(ils_size))

		! Rewind the file and read in ils_wave and ils
		REWIND(UNIT = 40)
		DO i = 1,ils_size
			READ(40,*,IOSTAT = ierr) ils_wave(i), ils(i)
		END DO
		CLOSE(40)

!------------------------------------------------
	! write(*,*)'DB_Pressure: ',DB_P
	CALL readDB(DB_P, nP, 1, 1)				! Output: absc_co2(i,iP,iT)
	
	DO n_profile = n_profile_b , n_profile_e
	write(*,*) 'n_profile= ', n_profile

    DO sample = sample_b , sample_e
	write(*,*) 'sample = ', sample
	
	T_slop = 0.d0
	xCO2_slop = 0.d0
	theta_slop = 0.d0
	fai_slop = 0.d0
	ALB_slop = 0.d0
	
	! T_slop = 8.d0*(ran1(idum1)-0.5d0)				!2 or 4
	! xCO2_slop = 50.d0*(ran1(idum1)-0.5d0)*1.d-6! 2 or 40
	! theta_slop = 0.2d0*(ran1(idum1)-0.5d0)
	! fai_slop = 0.2d0*(ran1(idum1)-0.5d0)
	! ALB_slop = 0.4d0*(ran1(idum1)-0.5d0) !0.02 or 0.04
	
	DO i=1,nsb
	  P(i) = seed(n_profile*4-2,i)/1000.d0
	  T(i) = seed(n_profile*4-1,i) + T_slop
	  fvCO2(i) = seed(n_profile*4,i)*1.d-6 + xCO2_slop
	  IF (i == 1) THEN
		dx(i) = seed(n_profile*4-3,i)*1.d5	!cm
	  ELSE
		dx(i) = (seed(n_profile*4-3,i) - seed(n_profile*4-3,i-1))*1.d5
	  END IF
	END DO

	ALLOCATE(absco(absco_size,nsb),absco_in(absco_size))
	ALLOCATE(dab_dx_in(absco_size))
	DO i=1,nsb
		CALL interp(P(i),T(i),fvCO2(i), fvH2O(i), fvCO(i), absco_wave, absco_size, absco_in(:), dab_dx_in(:))
		absco(:,i) = absco_in(:)
	END DO
	
	ALLOCATE(ils_spaced(INT((ils_wave(ils_size) - ils_wave(1))/ils_step)))
		! Now interpolate ils to be at 'ils_step' spacing
	i=1
	j=1
	k=1
	ils_spaced(1) = ils(1)
	ils_wave_ref = ils_wave(1)
	DO
		DO j = i,ils_size-1
			IF(ils_wave(j) - ils_wave_ref > ils_step) THEN		
				ils_interp(1) = ils(j-1)
				ils_interp(2) = ils(j)
				CALL LinTerpDouble(ils_interp(1), ils_interp(2), ils_spaced(k+1),&
												 ils_wave(j-1:j), ils_wave_ref + ils_step)
				ils_wave_ref = ils_wave_ref + ils_step
				EXIT
			ENDIF
		END DO
		IF(j >= ils_size .OR. k+1 == SIZE(ils_spaced)) EXIT
		i = j
		k = k+1
	END DO
	WRITE(*,*) 'Calculating transmittance....'
!------------------------------------------------OpricalThickness_Radiance
	theta = (seed_ALB(n_profile,1)*(1+theta_slop))*pi/180.d0	!solarzenith
	fai = (seed_ALB(n_profile,2)*(1+fai_slop))*pi/180.d0	!satelitezenith
	DO i = 5,24
		ALB(n_profile,i) = seed_ALB(n_profile,i) + seed_ALB(n_profile,i)*ALB_slop
	END DO
	
	Ns = 2.54743d19
	Ps = 1013.25d0
	Ts = 288.15d0
	ALLOCATE(label_rad(absco_size))
	ALLOCATE(opticalThickness(absco_size))
	ALLOCATE(rad(absco_size))
	ALLOCATE(ALBEDO(absco_size))

	k = 1
	q = 1
	DO i = 1,absco_size
		n(i) = ( 5791814.d0 / ( 238.0185d0 - (absco_wave(i)*1d-4)**2 )   +  167929.d0/(57.362d0-(absco_wave(i)*1d-4)**2) )*1d-8 + 1.d0
		rn(i) = 1.007482d-2 + 0.7990914d0/(47.48717d0 - (absco_wave(i)*1d-4)**2)
		Rayleigh(i,1) = 1.d0/(Ns* seed(4*n_profile-2,1)/Ps *Ts/seed(4*n_profile-1,1)) *24.d0* pi**3 *(n(i)**2-1.d0)**2 /&
			&(10000.d0/absco_wave(i))**4/Ns**2/(n(i)**2+2.d0)**2 * (6.d0+3.d0*rn(i)) / (6.d0-7.d0*rn(i))
		opticalThickness(i)=dx(1)*(absco(i,1) + Rayleigh(i,1))
		DO j=2,nsb
			Rayleigh(i,j) =1.d0/ (Ns* seed(4*n_profile-2,j)/Ps *Ts/seed(4*n_profile-1,j)) *24.d0* pi**3 *(n(i)**2-1.d0)**2 /&
				&(1d4/absco_wave(i))**4/(n(i)**2+2.d0)**2 * (6.d0+3.d0*rn(i))/(6.d0-7.d0*rn(i))
			opticalThickness(i)	= opticalThickness(i)+dx(j)* (absco(i,j) + Rayleigh(i,j))
		END DO
		IF ((i*ils_step + wave_bb) >= wave_43_b .AND.  (i * ils_step+wave_bb) <= wave_43_e) THEN
			DO m = q,19
				IF((i*ils_step + wave_bb) >= ALB_ref(m,1) .AND. (i*ils_step+wave_bb) < ALB_ref(m+1,1)) THEN
                                                                                Call LinTerpDouble(ALB(n_profile,m+4), ALB(n_profile,m+5), ALBEDO(i), ALB_ref(m:m+1,1), i*ils_step + wave_bb)
					rad(i) = cos(theta)*exp(-1.d0*opticalThickness(i)/cos(theta))*exp(-1.d0*opticalThickness(i)/cos(fai))*SO(k)/pi * ALBEDO(i)
					IF(m /= q) THEN
						q = q + 1
					ENDIF
					EXIT
				ELSEIF((i*ils_step+wave_bb) < ALB_ref(1,1)) THEN
					rad(i) = cos(theta)*exp(-1.d0*opticalThickness(i)/cos(theta))*exp(-1.d0*opticalThickness(i)/cos(fai))*SO(k)/pi * ALB(n_profile,5)
					ALBEDO(i) = ALB(n_profile,5)
					EXIT
				ELSEIF((i*ils_step+wave_bb) >= ALB_ref(20,1)) THEN
					rad(i) = cos(theta)*exp(-1.d0*opticalThickness(i)/cos(theta))*exp(-1.d0*opticalThickness(i)/cos(fai))*SO(k)/pi * ALB(n_profile,24)
					ALBEDO(i) = ALB(n_profile,24)
					EXIT
				ENDIF
			END DO
			k = k + 1
		ELSEIF((i*ils_step+wave_bb) > wave_43_e) THEN
			rad(i) = 0.d0!cos(theta)*exp(-1.d0*opticalThickness(i)/cos(theta))*exp(-1.d0*opticalThickness(i)/cos(fai))* 3.7418d-16*absco_wave(i)**3/(exp(1.43880d0*absco_wave(i)/5762.d0)-1.d0)*1.0d4/pi * ALB(n_profile,24)
			ALBEDO(i) = ALB(n_profile,24)
		ELSEIF((i*ils_step + wave_bb) < wave_43_b) THEN
			rad(i) = 0.d0!cos(theta)*exp(-1.d0*opticalThickness(i)/cos(theta))*exp(-1.d0*opticalThickness(i)/cos(fai))* 3.7418d-16*absco_wave(i)**3/(exp(1.43880d0*absco_wave(i)/5762.d0)-1.d0)*1.0d4/pi * ALB(n_profile,5)
			ALBEDO(i) = ALB(n_profile,5)
		ENDIF

		label_rad(i) = log( cos(theta)*ALBEDO(i)*exp(-1.d0*opticalThickness(i)/cos(theta))*exp(-1.d0*opticalThickness(i)/cos(fai)) )
		IF(label_rad(i) < -99999) THEN
			label_rad(i) = -999999999.d0
		ENDIF
	END DO ! i

WRITE(*,*) k
!------------------------------------------------Convolution
convsize = absco_size + SIZE(ils_spaced)
	
DO conv = 1,4
	ALLOCATE(transmit(absco_size))
	ALLOCATE(ils_spaced1(INT((ils_wave(ils_size) - ils_wave(1))/ils_step)))
	ALLOCATE(ils_ft(convsize/2 + 1), tran_ft(convsize/2 + 1),conv_in(convsize/2 + 1))
	ALLOCATE(tran_interim(convsize),tran_conv(convsize))
	ALLOCATE(tran_compare(2,wavesize),wvnm_compare(2,wavesize))
	ALLOCATE(tran_final(wavesize))
	if (conv == 1) THEN
		transmit = rad
	ELSEIF(conv == 2) THEN
		transmit = log(ALBEDO)
		DO i = 1,absco_size
			IF(transmit(i) < -99999) THEN
				transmit(i) = -999999999.d0
			ENDIF
		END DO
	ELSEIF(conv == 3) THEN
		transmit = label_rad
	ELSEIF(conv == 4) THEN
		transmit = opticalThickness
	ENDIF
	ils_spaced1 = ils_spaced

	! This loop goes around only twice. The first time is to calculate the
	! transmittance with the actual ILF, the second time with the ideal ILF.
	  actorideal:DO s=2,2	
	! Allocate ils_spaced

	IF(convmethod) THEN
		! Need to size transmit and ils for FFTW. For a circular convolution, 
		! we want convsize = SIZE(transmit) + SIZE(ils_spaced)

		transmit => reallocate(transmit,convsize)
		ils_spaced1 => reallocate(ils_spaced1,convsize)
		!write(*,*)SIZE(ils_spaced),absco_size,convsize

		! Call FFTW to get Fourier transform for sinc2(s) and intensity(iT,i)
		! Visit http://fftw.org/fftw3_doc/ for information regarding FFTW.

		CALL dfftw_plan_dft_r2c_1d(plan1,convsize,ils_spaced1,ils_ft,FFTW_ESTIMATE)
		CALL dfftw_plan_dft_r2c_1d(plan2,convsize,transmit,tran_ft,FFTW_ESTIMATE)

		CALL dfftw_execute_dft_r2c(plan1, ils_spaced1, ils_ft)
		CALL dfftw_execute_dft_r2c(plan2, transmit, tran_ft)

		CALL dfftw_destroy_plan(plan1)
		CALL dfftw_destroy_plan(plan2)

		! In frequency domain, convolution is simply a multiplication. This is why
		! arrays must be of the same size!
		DO i = 1,convsize/2 + 1
			 conv_in(i) = ils_ft(i)*tran_ft(i)/convsize*ils_step
		END DO

		! Now call inverse FT to get back to 'time' domain. 'tran_conv' is the
		! convolution of transmittance with ILF.
		CALL dfftw_plan_dft_c2r_1d(plan1,convsize,conv_in,tran_conv,FFTW_ESTIMATE)
		CALL dfftw_execute_dft_c2r(plan1, conv_in, tran_conv)
		CALL dfftw_destroy_plan(plan1)
	ELSE
		convsize = absco_size + SIZE(ils_spaced)
		transmit => reallocate(transmit,convsize)
		ils_spaced1 => reallocate(ils_spaced1,convsize)

		DO i = 1,convsize
			tran_conv(i) = 0.d0
			DO j = 1,convsize
				IF(i-j+1 > 0) THEN
					tran_conv(i) = tran_conv(i) + transmit(i-j+1)*ils_spaced1(j)
				ELSE
					EXIT
				ENDIF
			END DO
		END DO
	ENDIF	

		! Now that tran_conv has been calculated, cut off unecessary signal.
		! Note: NEED TO UPDATE IF STATEMENT BELOW...USE SOMETHING MORE ROBUST THAN 
		! THE VALUE 2350.d0

		j = 1
		DO i = 1,convsize
			IF(i > ((convsize + absco_size)/2) .OR.&
									 i < ((convsize - absco_size)/2)) THEN
				CYCLE
			ELSE
				tran_interim(j) = tran_conv(i)
				j = j+1
			ENDIF
		END DO

!------------------------------------------------Interpolation

		! We will now interpolate tran_interim so the spacing is the same as the 
		! FTIR data.

		i=2
		j=1
		k=1
		tran_final(1) = tran_interim(1)
		wave_ref = wvnm(1) 		! The first wavenumber of ils_spaced is same as ils
		wave_step = wvnm(2) - wvnm(1)
		DO
			 DO j = i,absco_size
			    IF(absco_wave(j) - wave_ref > (wave_step)) THEN
			       tran_interp(1) = tran_interim(j-1)			
			       tran_interp(2) = tran_interim(j)
			       CALL LinTerpDouble(tran_interp(1), tran_interp(2),&
								 tran_final(k+1), absco_wave(j-1:j), wave_ref + wave_step)
			       wave_ref = wave_ref + wave_step

			       EXIT
			    ENDIF
			 END DO
			 IF(j >= absco_size .OR. k == wavesize-1) EXIT
			 i = j
			 k = k+1
		END DO
		j = 1
		DO i = 1,wavesize
			IF(wvnm(i) >= wave_bb .AND. wvnm(i) <= wave_ee) THEN
				tran_compare(s,j) = tran_final(i)
				wvnm_compare(s,j) = wvnm(i)
				j = j+1
			ELSE
				CYCLE
			ENDIF
		END DO
		! Deallocate variables
	END DO actorideal

!------------------------------------------------Output
	i_43_b = INT((6180.d0 - wvnm_compare(2,1))/resolution) + 1 + 1.E-6
	i_43_e = INT((6280.d0 - wvnm_compare(2,1))/resolution) + 1 + 1.E-6

	If (conv == 1) THEN
!!!!! (18:labels_FIT)
		WRITE(18,"(*(E15.9,:,','))") (tran_compare(2,i) , i= i_43_b , i_43_e)
	ELSEIF(conv == 2) THEN
		ALLOCATE(ALBEDO_ILSF(wavesize))
		ALBEDO_ILSF = tran_compare(2,:)
	ELSEIF(conv == 3) THEN
!!!!! (19:labels_log)	(20:labels_log_OT_ILSF)
		WRITE(19,"(*(E15.9,:,','))") (tran_compare(2,i), i= i_43_b , i_43_e)
		ALLOCATE(labels_OT(i_43_e - i_43_b + 1))
		q = 12
		DO i = i_43_b,i_43_e
			DO m = q,16
				IF((i*resolution+wvnm_compare(2,1)) >= ALB_ref(m,1) .AND. (i*resolution+wvnm_compare(2,1)) < ALB_ref(m+1,1)) THEN
					labels_OT(i+1-i_43_b) = ( log(cos(theta)) + ALBEDO_ILSF(i) - tran_compare(2,i) ) / (cos(theta)+cos(fai))*cos(theta)*cos(fai)
					IF(m /= q) THEN
						q = q + 1
				ENDIF
					EXIT
				ELSEIF((i*resolution+wvnm_compare(2,1)) < ALB_ref(12,1)) THEN
					labels_OT(i+1-i_43_b) = ( log(cos(theta)) + ALBEDO_ILSF(i) - tran_compare(2,i) ) / (cos(theta)+cos(fai))*cos(theta)*cos(fai)
					EXIT
				ELSEIF((i*resolution+wvnm_compare(2,1)) >= ALB_ref(16,1)) THEN
					labels_OT(i+1-i_43_b) = ( log(cos(theta)) + ALBEDO_ILSF(i) - tran_compare(2,i) ) / (cos(theta)+cos(fai))*cos(theta)*cos(fai)
					EXIT
				ENDIF
			END DO
		END DO
		WRITE(20,"(*(E15.9,:,','))") (labels_OT(i), i = 1 , size(labels_OT))
		DEALLOCATE(labels_OT)

	ELSEIF(conv == 4) THEN
!!!!! (24:OT_ILSF)
		WRITE(24,"(*(E15.9,:,','))") (tran_compare(2,i), i= i_43_b , i_43_e)
	ENDIF

	WRITE(*,*) i_43_b , i_43_e
	WRITE(*,*) "Begin, Begin+1, End, Resolution: ",wvnm_compare(2,i_43_b),wvnm_compare(2,i_43_b+1),wvnm_compare(2,i_43_e),resolution

	DEALLOCATE(ils_spaced1,transmit)
	DEALLOCATE(tran_interim, tran_conv,ils_ft, tran_ft, conv_in)
	DEALLOCATE(tran_compare,wvnm_compare)
	DEALLOCATE(tran_final)

END DO !!conv		

	WRITE(21,"(*(F8.4,:,','))") (T(i), i=1,nsb)
	WRITE(22,"(*(F12.10,:,','))") (fvCO2(i), i=1,nsb)
	WRITE(28,"(*(F9.6,:,','))") (ALBEDO_ILSF(i) , i=i_43_b,i_43_e) 
	WRITE(29,"(F8.6,',')") (theta,fai)

	DEALLOCATE(ils_spaced)
	DEALLOCATE(absco,absco_in)
	DEALLOCATE(dab_dx_in)
	DEALLOCATE(rad,opticalThickness,label_rad)
	DEALLOCATE(ALBEDO,ALBEDO_ILSF)

	END DO !SAMPLE
	END DO !PROFILE
	DEALLOCATE(seed,seed_ALB)
	DEALLOCATE(absco_wave,wvnm)
	DEALLOCATE(ils_wave)

	CLOSE(13)
	CLOSE(14)
	CLOSE(18)
	CLOSE(19)
	CLOSE(20)
	CLOSE(21)
	CLOSE(22)
	CLOSE(24)
	CLOSE(28)
	CLOSE(29)
	END PROGRAM testsinc8


    FUNCTION gasdev(idum)
! Gaussian random number with 0 mean and unity standard variation
	implicit none
    INTEGER		:: idum
    REAL(8)		:: gasdev
!    USES ran1
    INTEGER		:: iset
    REAL(8)		:: fac,gset,rsq,v1,v2,ran1
    SAVE		:: iset,gset
    DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END FUNCTION

	FUNCTION ran1(idum)
 	implicit none
    INTEGER		:: idum,NDIV
    REAL(8)		:: ran1,AM,EPS,RNMX
    INTEGER, PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32
    INTEGER		:: j,k,iv(NTAB),iy
    SAVE iv,iy
    DATA iv /NTAB*0/, iy /0/
	AM = 1./IM;NDIV=1+(IM-1)/NTAB;EPS=1.2e-7;RNMX=1.-EPS
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END FUNCTION
