INCLUDE './absc_kg_gen.h'
MODULE absc_kg_gen
USE commonData

! Variables and Constants
  INTEGER,PARAMETER,PRIVATE :: rows=1!14241164 ! Max number of lines in DB for H2O
! INTEGER,PARAMETER,PRIVATE :: rows=14241164 ! Max number of lines in DB for CO2
  REAL(8),ALLOCATABLE,PRIVATE :: DATA(:,:)
  INTEGER,PRIVATE           :: lines

  PUBLIC :: GetData,absc_setup,absc_gen

  LOGICAL,PUBLIC :: Pressure_based= .false.             ! linear or pressure-based Planck mean absco

CONTAINS

!-----------------------------------------------------------
SUBROUTINE maincalc(Gas_Info,eta_b,eta_e,step,absc)
IMPLICIT NONE

REAL(8), INTENT(IN) :: eta_b, eta_e
REAL(8), INTENT(IN) :: step
INTEGER :: number,i,j,it
TYPE(GasInfo) :: Gas_Info 
REAL(4),DIMENSION(:),allocatable :: absc
REAL(8),DIMENSION(:),ALLOCATABLE :: absc1

! Define length of absco array
number=(eta_e-eta_b)/step+2+1.d-6
ALLOCATE(absc(number), absc1(number))

absc=0.d0   ! Initialize absc
absc1 = 0.d0

WRITE(*,*) 'Calculating absorption coefficient'
   IF (Gas_Info%x > 0.d0) THEN
      CALL absc_setup(Gas_Info)
      CALL absc_gen(Gas_Info, absc1(:), step, .TRUE.)
   ELSE
      RETURN
   ENDIF
absc = absc1

END SUBROUTINE maincalc
!------------------------------------------------------
  SUBROUTINE absc_setup(Gas_Info)
!!$ Prepare for the absorption coefficient generation, reading spectral lines
!!$ from the spectroscopic database. If 
!!$ Input:
!!$   GasInfo    --  information of the gas
    IMPLICIT NONE
    TYPE(GasInfo),INTENT(IN) :: Gas_Info
    INTEGER  :: ifg,ierr=0
    ifg = Gas_Info%Gas
    ALLOCATE(DATA(rows,6),STAT= ierr)
    IF (ierr>0) WRITE(*,*) 'Insufficient Memory(Cannot allocate data).'
    CALL Getdata(Gas_Info,lines)
  END SUBROUTINE absc_setup

!-------------------------------------------------------
  SUBROUTINE absc_free
!!$ free memory after the absorption coefficient generation
    IMPLICIT NONE
    INTEGER :: ierr=0
    DEALLOCATE(DATA,STAT= ierr)
    IF(ierr>0) WRITE(*,*)'Unsuccessful memory free(data).'
  END SUBROUTINE absc_free

!!$------------------------------------------------------------------------

  SUBROUTINE absc_gen(Gas_Info,absc,wvnmst,voigt_pfl)
!!$ Calculate the absorption coefficients at given (P,T,x) for the full spectrum
!!$   Input:
!!$     Gas_Info     --  information of gas
!!$     wvnmst       --  wavenumber step, output of subroutine 'wvnm_resolution'
!!$     voigt_pfl    --  (optional) If 'True', Voigt profile calculation is used for every line
!!$                      If 'False', it is decided automatically to use either Voigt or Lorentz profile
!!$                      Default value is 'False'
!!$   Output:
!!$     absc         --  absorption coefficients calculated for the full spectrum
    IMPLICIT NONE
    TYPE(GasInfo),INTENT(IN)  :: Gas_Info
    REAL(8),INTENT(OUT)       :: absc(:)
    REAL(8),INTENT(IN)        :: wvnmst
    LOGICAL,INTENT(IN),OPTIONAL :: voigt_pfl
!!$ local variables
    REAL(8),PARAMETER :: hhh=6.626076d-34,ccc=2.997925d10 ! hhh = Planck's constant (Js); ccc = speed of light in vaccuum (cm/s)
    REAL(8),PARAMETER :: kkk=1.380658d-23,nnn=6.022136d23 ! kkk = boltzman constant (J/K); nnn = Avogadro's # (mol^-1)
    REAL(8),PARAMETER :: PI=3.14159265357989d0
    REAL(8),PARAMETER :: c1=3.7419d-12,c2=1.4388d0,sigma=5.67d-12
    REAL(8)           :: hck,hckt,hckt0,hcktt0,patm,T0,klmin
    REAL(8)           :: V,V0,VR0,VV0,RR0,dk,deta,mass,bd0,bd
    REAL(8)           :: intensity,line_width,klmax,TT0,cr,crx,dummy
    REAL(8)           :: wvnm
    INTEGER           :: ifg,i,icl,l,number
    LOGICAL           :: voigt_profile
    REAL(8)           :: T,P,xmfr

    IF(PRESENT(voigt_pfl)) THEN
       voigt_profile= voigt_pfl
    ELSE
       voigt_profile= .FALSE.
    ENDIF
    T= Gas_Info%T; P= Gas_Info%P; xmfr= Gas_Info%x; ifg= Gas_Info%Gas
    number= (wave_e-wave_b)/wvnmst+1+1.d-6

    IF(ifg==CO2) THEN
      T0= 296.0d0
       mass= 44.d-3/nnn
    ELSE IF(ifg==H2O) THEN
       T0= 296.0d0
       mass= 18.d-3/nnn
    ELSE IF (ifg==CO) THEN
       T0= 296.0d0
       mass= 28.d-3/nnn
    ELSE IF(ifg==CH4) THEN
       T0= 296.0d0
       mass= 16.d-3/nnn
    ELSE IF(ifg==C2H4) THEN
       T0= 296.0d0
       mass= 28.05d-3/nnn
    ELSE IF(ifg==NO) THEN
       T0 = 296.0d0
       mass = 30.d-3/nnn
    ELSE IF(ifg==OH) THEN
       T0 = 296.0d0
       mass = 17.d-3/nnn
    END IF
    bd0=1.d2* DSQRT(2.d0*kkk*T*LOG(2.d0)/mass)/ccc
    hck=hhh*ccc/kkk
    hckt0=hck/T0
    patm= P/1.01325

!!$ absco is additive
    hckt=hck/T
    hcktt0=hckt0-hckt
    cr=nnn/8.314d1/T
    crx=cr*xmfr*Patm
    IF(pressure_based) THEN
       klmin= 1.d-9
    ELSE
       klmin= 1.d-9*xmfr*P
    ENDIF

    IF(ifg==1) THEN
       !!$ VIBRATIONAL PARTITION FUNCTION FOR CO2
       V0=(1d0-DEXP(-666d0*hckt0))**2*(1d0-DEXP(-2396d0*hckt0))*(1d0-DEXP(-1351d0*hckt0))
       V=(1d0-DEXP(-666d0*hckt))**2*(1d0-DEXP(-2396d0*hckt))*(1d0-DEXP(-1351d0*hckt))
       VV0=V/V0
       !!$ ROTATIONAL PARTITION FUNCTION
       TT0=T0/T
       RR0=TT0
    ELSE IF(ifg==2) THEN
       !!$ VIBRATIONAL PARTITION FUNCTION FOR H2O
       V0=(1d0-DEXP(-3652d0*hckt0))*(1d0-DEXP(-1595d0*hckt0))*(1d0-DEXP(-3756d0*hckt0))
       V=(1d0-DEXP(-3652d0*hckt))*(1d0-DEXP(-1595d0*hckt))*(1d0-DEXP(-3756d0*hckt))
       VV0=V/V0
       !!$ ROTATIONAL PARTITION FUNCTION
       TT0=T0/T
       RR0=TT0**1.5
    ELSE IF (ifg==3) THEN
       !!$ VIBRATIONAL PARTITION FUNCTION FOR CO
       V0=1d0-DEXP(-2143d0*hckt0)
       V=1d0-DEXP(-2143d0*hckt)
       VV0=V/V0
       !!$ ROTATIONAL PARTITION FUNCTION
       TT0=T0/T
       RR0=TT0    
    ELSE IF(ifg==4) THEN
       !!$ VIBRATIONAL PARTITION FUNCTION FOR CH4
       V0=(1d0-DEXP(-2917d0*hckt0))*(1d0-DEXP(-1534d0*hckt0))**2*   &
            (1d0-DEXP(-3019d0*hckt0))**3*(1d0-DEXP(-1306d0*hckt0))**3
       V=(1d0-DEXP(-2917d0*hckt))*(1d0-DEXP(-1534d0*hckt))**2*   &
            (1d0-DEXP(-3019d0*hckt))**3*(1d0-DEXP(-1306d0*hckt))**3
       VV0=V/V0
       !!$ ROTATIONAL PARTITION FUNCTION
       TT0=T0/T
       RR0=TT0**1.5
    ELSE IF(ifg==5) THEN
       !!$ VIBRATIONAL PARTITION FUNCTION FOR C2H4
       V0=(1d0-DEXP(-3026d0*hckt0))*(1d0-DEXP(-1623d0*hckt0))*(1d0-DEXP(-1342d0*hckt0))* &
            (1d0-DEXP(-1023d0*hckt0))*(1d0-DEXP(-3103d0*hckt0))*(1d0-DEXP(-1236d0*hckt0))* &
            (1d0-DEXP(-949d0*hckt0))*(1d0-DEXP(-943d0*hckt0))*(1d0-DEXP(-3106d0*hckt0))* &
            (1d0-DEXP(-826d0*hckt0))*(1d0-DEXP(-2989d0*hckt0))*(1d0-DEXP(-1444d0*hckt0))
       V= (1d0-DEXP(-3026d0*hckt))*(1d0-DEXP(-1623d0*hckt))*(1d0-DEXP(-1342d0*hckt))* &
            (1d0-DEXP(-1023d0*hckt))*(1d0-DEXP(-3103d0*hckt))*(1d0-DEXP(-1236d0*hckt))* &
            (1d0-DEXP(-949d0*hckt))*(1d0-DEXP(-943d0*hckt))*(1d0-DEXP(-3106d0*hckt))*	&
            (1d0-DEXP(-826d0*hckt))*(1d0-DEXP(-2989d0*hckt))*(1d0-DEXP(-1444d0*hckt))
       VV0=V/V0
       !!$ ROTATIONAL PARTITION FUNCTION
       TT0=T0/T
       RR0=TT0**1.5
    ELSEIF(ifg == 6) THEN
       !!$ VIBRATIONAL PARTITION FUNCTION FOR NO
       V0=1d0-DEXP(-1876d0*hckt0)
       V=1d0-DEXP(-1876d0*hckt)
       VV0=V/V0
       !!$ ROTATIONAL PARTITION FUNCTION
       TT0=T0/T
       RR0=TT0
    ELSEIF(ifg == 7) THEN
       !!$ VIBRATIONAL PARTITION FUNCTION FOR OH
       V0=1d0-DEXP(-3569d0*hckt0)
       V=1d0-DEXP(-3569d0*hckt)
       VV0=V/V0
       !!$ ROTATIONAL PARTITION FUNCTION
       TT0=T0/T
       RR0=TT0
    ENDIF

    VR0=VV0*RR0
!!$ Scan over lines
    DO l=1,lines
       line_width=(DATA(l,4)*xmfr+DATA(l,3)*(1.d0-xmfr))*patm*TT0**DATA(l,6)
!!$ molecule-based intensity
       intensity=DATA(l,2)*VR0*DEXP(hcktt0*DATA(l,5))*(1d0-DEXP(-DATA(l,1)*hckt))/&
            (1d0-DEXP(-DATA(l,1)*hckt0))
       IF (Pressure_Based) THEN
!!$ pressure-based intensity
          intensity=intensity*cr
       ELSE
!!$ linear intensity
          intensity= intensity*crx
       ENDIF

!!$ absorption coefficient at line center
       klmax=intensity/(PI*line_width)
!!$ Doppler line half-width
       bd= DATA(l,1)*bd0
       IF(klmax<klmin) CYCLE
!!$ Find wavenumber (subscript) closest to line center
       icl=(DATA(l,1)-wave_b)/wvnmst+1
       IF (icl<1) icl=1
       IF (icl>number) icl=number  ! if the line is outside the range start it's contribution at the first or last wavenumber 
!!$ Scan over adjacent wavenumbers to see whether line makes contribution to absco
       DO i=icl,number
          wvnm= wave_b+(i-1)*wvnmst
          deta=(DATA(l,1)-wvnm)/line_width
          IF (voigt_profile.OR.(line_width/bd < 10.0)) THEN
             CALL VOIGT(intensity,line_width,bd,DATA(l,1)-wvnm,dk)
          ELSE
             dk=klmax/(1.d0+deta*deta)
          ENDIF
          absc(i)=absc(i)+dk
          IF(dk<klmin) EXIT
       ENDDO !i
       DO i=icl-1,1,-1
          wvnm= wave_b+(i-1)*wvnmst
          deta=(DATA(l,1)-wvnm)/line_width
          IF (voigt_profile.OR.(line_width/bd < 10.0)) THEN
             CALL VOIGT(intensity,line_width,bd,DATA(l,1)-wvnm,dk)
          ELSE
             dk=klmax/(1.d0+deta*deta)
          ENDIF
          absc(i)=absc(i)+dk
          IF(dk<klmin) EXIT
       ENDDO !i
    ENDDO ! l

  END SUBROUTINE absc_gen

!!$*************************************************************************
  SUBROUTINE Getdata(Gas_Info,i)
!!$ Read the spectral lines from HITEMP for H2O
!!$ Read the spectral lines from HITRAN for CH4,C2H4 and CO 
!!$   Output:     i  --  number of lines read

    USE hitemp_hitran_cdsd
    IMPLICIT NONE

    TYPE(GasInfo),INTENT(IN)  :: Gas_Info
    INTEGER :: i,j,c,h,ifg,ier,dummy1, no_files, begoffile, endoffile
    INTEGER, PARAMETER :: lu=269
    CHARACTER :: dummy2*10
    CHARACTER(40), ALLOCATABLE :: strDataFileName(:)
    CHARACTER(80) :: strFilePath
    
    ifg = Gas_Info%gas
    IF (ifg == CO2) THEN
        c = 3                   ! Change value of c according to database (1 = hitran2008,
        IF (c == 1) THEN        ! 2 = cdsd2008, 3 = hitemp2010)
           strFilePath = hitran_CO2
           no_files = 1
           ALLOCATE(strDataFileName(no_files))

           strDataFileName(1) = "02_hit08.par"                        
        ELSEIF (c == 2) THEN
           strFilePath = cdsd_CO2
           no_files = 20
           ALLOCATE(strDataFileName(no_files))
           
           strDataFileName(1) = "cdsd_1000_01"
           strDataFileName(2) = "cdsd_1000_02"
           strDataFileName(3) = "cdsd_1000_03"
           strDataFileName(4) = "cdsd_1000_04"
           strDataFileName(5) = "cdsd_1000_05"
           strDataFileName(6) = "cdsd_1000_06"
           strDataFileName(7) = "cdsd_1000_07"
           strDataFileName(8) = "cdsd_1000_08"
           strDataFileName(9) = "cdsd_1000_09"
           strDataFileName(10) = "cdsd_1000_10"
           strDataFileName(11) = "cdsd_1000_11"
           strDataFileName(12) = "cdsd_1000_12"
           strDataFileName(13) = "cdsd_1000_13"
           strDataFileName(14) = "cdsd_1000_14"
           strDataFileName(15) = "cdsd_1000_15"
           strDataFileName(16) = "cdsd_1000_16"
           strDataFileName(17) = "cdsd_1000_17"
           strDataFileName(18) = "cdsd_1000_18"
           strDataFileName(19) = "cdsd_1000_19"
           strDataFileName(20) = "cdsd_1000_20"

        ELSEIF (c == 3) THEN
           strFilePath = hitemp_CO2
           no_files = 20
           ALLOCATE(strDataFileName(no_files))
   
           strDataFileName(1) = "02_0-500_HITEMP2010.par"
           strDataFileName(2) = "02_500-625_HITEMP2010.par"
           strDataFileName(3) = "02_625-750_HITEMP2010.par"
           strDataFileName(4) = "02_750-1000_HITEMP2010.par"
           strDataFileName(5) = "02_1000-1500_HITEMP2010.par"
           strDataFileName(6) = "02_1500-2000_HITEMP2010.par"
           strDataFileName(7) = "02_2000-2125_HITEMP2010.par"
           strDataFileName(8) = "02_2125-2250_HITEMP2010.par"
           strDataFileName(9) = "02_2250-2500_HITEMP2010.par"
           strDataFileName(10) = "02_2500-3000_HITEMP2010.par"
           strDataFileName(11) = "02_3000-3250_HITEMP2010.par"
           strDataFileName(12) = "02_3250-3500_HITEMP2010.par"
           strDataFileName(13) = "02_3500-3750_HITEMP2010.par"
           strDataFileName(14) = "02_3750-4000_HITEMP2010.par"
           strDataFileName(15) = "02_4000-4500_HITEMP2010.par"
           strDataFileName(16) = "02_4500-5000_HITEMP2010.par"
           strDataFileName(17) = "02_5000-5500_HITEMP2010.par"
           strDataFileName(18) = "02_5500-6000_HITEMP2010.par"
           strDataFileName(19) = "02_6000-6500_HITEMP2010.par"
           strDataFileName(20) = "02_6500-12785_HITEMP2010.par"
        ENDIF

    ELSE IF (ifg == H2O) THEN
        h = 2
        IF(h == 1) THEN
           strFilePath = hitran_H2O
           no_files = 1
           ALLOCATE(strDataFileName(no_files))
           strDataFileName(1) = "01_hit08.par"

        ELSEIF(h == 2) THEN
        
           strFilePath = hitemp_H2O
   
           no_files = 34
   
           ALLOCATE(strDataFileName(no_files))
   
           strDataFileName(1) = "01_0-50_HITEMP2010.par"
           strDataFileName(2) = "01_50-150_HITEMP2010.par"
           strDataFileName(3) = "01_150-250_HITEMP2010.par"
           strDataFileName(4) = "01_250-350_HITEMP2010.par"
           strDataFileName(5) = "01_350-500_HITEMP2010.par"
           strDataFileName(6) = "01_500-600_HITEMP2010.par"
           strDataFileName(7) = "01_600-700_HITEMP2010.par"
           strDataFileName(8) = "01_700-800_HITEMP2010.par"
           strDataFileName(9) = "01_800-900_HITEMP2010.par"
           strDataFileName(10) = "01_900-1000_HITEMP2010.par"
           strDataFileName(11) = "01_1000-1150_HITEMP2010.par"
           strDataFileName(12) = "01_1150-1300_HITEMP2010.par"
           strDataFileName(13) = "01_1300-1500_HITEMP2010.par"
           strDataFileName(14) = "01_1500-1750_HITEMP2010.par"
           strDataFileName(15) = "01_1750-2000_HITEMP2010.par"
           strDataFileName(16) = "01_2000-2250_HITEMP2010.par"
           strDataFileName(17) = "01_2250-2500_HITEMP2010.par"
           strDataFileName(18) = "01_2500-2750_HITEMP2010.par"
           strDataFileName(19) = "01_2750-3000_HITEMP2010.par"
           strDataFileName(20) = "01_3000-3250_HITEMP2010.par"
           strDataFileName(21) = "01_3250-3500_HITEMP2010.par"
           strDataFileName(22) = "01_3500-4150_HITEMP2010.par"
           strDataFileName(23) = "01_4150-4500_HITEMP2010.par"
           strDataFileName(24) = "01_4500-5000_HITEMP2010.par"
           strDataFileName(25) = "01_5000-5500_HITEMP2010.par"
           strDataFileName(26) = "01_5500-6000_HITEMP2010.par"
           strDataFileName(27) = "01_6000-6500_HITEMP2010.par"
           strDataFileName(28) = "01_6500-7000_HITEMP2010.par"
           strDataFileName(29) = "01_7000-7500_HITEMP2010.par"
           strDataFileName(30) = "01_7500-8000_HITEMP2010.par"
           strDataFileName(31) = "01_8000-8500_HITEMP2010.par"
           strDataFileName(32) = "01_8500-9000_HITEMP2010.par"
           strDataFileName(33) = "01_9000-11000_HITEMP2010.par"
           strDataFileName(34) = "01_11000-30000_HITEMP2010.par"
        ENDIF

    ELSE IF (ifg==CO) THEN
 
       strFilePath = hitemp_CO
       no_files = 1
       ALLOCATE(strDataFileName(no_files))
       strDataFileName(1) = "05_HITEMP2010new.par"
!		strDataFileName(1) = "05_hit04.par"

    ELSE IF (ifg==CH4) THEN

       strFilePath = hitran_CH4
       no_files = 1
       ALLOCATE(strDataFileName(no_files))
       strDataFileName(1) = "06_hit08.par"

    ELSE IF (ifg==C2H4) THEN

       strFilePath = hitran_C2H4
       no_files = 1
       ALLOCATE(strDataFileName(no_files))
       strDataFileName(1) = "38_hit08.par"
 
    ELSE IF (ifg==NO) THEN

       strFilePath = hitran_NO
       no_files = 1
       ALLOCATE(strDataFileName(no_files))
       strDataFileName(1) = "08_hit08.par"

    ELSE IF (ifg==OH) THEN

       strFilePath = hitran_NO
       no_files = 1
       ALLOCATE(strDataFileName(no_files))
       strDataFileName(1) = "13_hit08.par"

    END IF

    IF(ier==1) THEN
       IF (ifg==1) THEN
          WRITE(*,*) '!!!Fatal Error: CO2 HITEMP database does not exist!!!'
       ELSE IF (ifg==2) THEN
          WRITE(*,*) '!!!Fatal Error: H2O HITEMP database does not exist!!!'
       ELSE IF (ifg==3) THEN
          WRITE(*,*) '!!!Fatal Error: CO HITEMP database does not exist!!!'
       ELSE IF (ifg==4) THEN
          WRITE(*,*) '!!!Fatal Error: CH4 HITRAN database does not exist!!!'
       ELSE IF (ifg==5) THEN
          WRITE(*,*) '!!!Fatal Error: C2H4 HITRAN database does not exist!!!'
       ELSE IF (ifg==6) THEN
          WRITE(*,*) '!!!Fatal Error: NO HITRAN database does not exist!!!'
       ELSE IF (ifg==7) THEN
          WRITE(*,*) '!!!Fatal Error: OH HITRAN database does not exist!!!'
       END IF
       STOP
    ENDIF
 

    ! data(i,1) = wavenumber
    ! data(i,2) = intensity
    ! data(i,3) = b_air
    ! data(i,4) = b_self
    ! data(i,5) = E''
    ! data(i,6) = exponent for b

    IF (ifg == CO2 .OR. ifg == H2O) THEN
        ! Find beginning file. This prevents the code from needing to cycle through 
        ! unnecessary files.

        beginning: DO i = 1, no_files
           OPEN(lu, FILE = TRIM(ADJUSTL(ADJUSTR(strFilePath)//ADJUSTL(strDataFileName(i)))), &
                STATUS = 'OLD', ACTION = 'READ',IOSTAT = ier)

           READ(lu,FMT=160,IOSTAT = ier) dummy1,DATA(1,1)

           IF (DATA(1,1) > wave_b-250.d0 .AND. i>1) THEN
               begoffile = i-1
               EXIT
           ELSE
               begoffile = i
           ENDIF 
           CLOSE(lu)
        END DO beginning
        160 FORMAT(I3, F12.6)
 
        ! Find ending file. This prevents the code from cycling through unnecessary files.
        ending: DO i = 1, no_files
           OPEN(lu, FILE = TRIM(ADJUSTL(ADJUSTR(strFilePath)//ADJUSTL(strDataFileName(i)))), &
                STATUS = 'OLD', ACTION = 'READ',IOSTAT = ier)

           READ(lu,FMT=160,IOSTAT = ier) dummy1,DATA(1,1)

           IF (DATA(1,1) > wave_e + 250.d0 .AND. i>1) THEN
               endoffile = i-1
               EXIT
           ELSE
               endoffile = i
           ENDIF
           CLOSE(lu)
        END DO ending

    ELSE
        begoffile = 1
        endoffile = 1
    ENDIF

! Read Database
i = 1
    gas: DO j = begoffile, endoffile

       OPEN(lu, FILE = TRIM(ADJUSTL(ADJUSTR(strFilePath)//ADJUSTL(strDataFileName(j)))), &
            STATUS = 'OLD', ACTION = 'READ',IOSTAT = ier)

       write(*,*)TRIM(ADJUSTL(ADJUSTR(strFilePath)//ADJUSTL(strDataFileName(j))))

       wave: DO
           READ(lu,FMT=150,IOSTAT = ier) dummy1,DATA(i,1),DATA(i,2),dummy2,DATA(i,3),    &
                DATA(i,4),DATA(i,5),DATA(i,6)

           IF(ier < 0) THEN  ! Reached end of file
              IF(i>1) i = i-1
              close(lu)
              EXIT
           ENDIF

           IF(DATA(i,1)<wave_b-250.d0) CYCLE

           IF(DATA(i,1)>wave_e+250.d0) THEN
              i = i-1
              close(lu)
              EXIT
           ENDIF

           iF(DATA(i,5)<0.d0) THEN
              DATA(i,5)=0.d0
           ENDIF

           IF(DATA(i,4)<1.e-5) DATA(i,4)=5.d0*DATA(i,3)

!              SELECT CASE(dummy1)
!                 CASE(21)
!                    DATA(i,2) = DATA(i,2)*0.984204d0
!                 CASE(22)
!                    DATA(i,2) = DATA(i,2)*0.0d0
!                 CASE(23)
!                    DATA(i,2) = DATA(i,2)*0.0d0
!                 CASE(24)
!                    DATA(i,2) = DATA(i,2)*0.0d0
!                 CASE(25)
!                    DATA(i,2) = DATA(i,2)*0.0d0
!                 CASE(26)
!                    DATA(i,2) = DATA(i,2)*0.0d0
!                 CASE(27)
!                    DATA(i,2) = DATA(i,2)*0.0d0
!              END SELECT
           i=i+1
       END DO wave

    ! ES10.3 used to be E10.3 in prior code and F5.3 was F5.4
    150 FORMAT(I3,F12.6,ES10.3,A10,F5.4,F5.3,F10.4,F4.2)
    END DO gas
    RETURN

  END SUBROUTINE Getdata
!------------------------------------------------------
!  REAL(8) FUNCTION BBFN(Q)
!!!$     ********************************************************************
!!!$     *  This subroutine calculates the fractional blackbody             *
!!!$     *  emissive power f(n*lambda*T), where X=n*lambda*T in (micro-m*K) *
!!!$     ********************************************************************
!    REAL(8) :: PI,CC,C2,EPS,V,EX,M,EM,VM,BM,Q
!    REAL(8),PARAMETER :: x1= 21.d0, x2= 75671.d0
!
!    IF (Q <x1) THEN        ! if Q is too small, return 0.0
!       BBFN=0.d0; RETURN
!    ELSEIF (Q>x2) THEN     ! if Q is too large, return 1.0
!       BBFN= 1.d0; RETURN
!    ENDIF
!
!    PI=4.D0*DATAN(1.D0)
!    CC=1.5D1/PI**4
!    C2=1.4388D4
!    EPS=1.D-16
!
!    V=C2/Q
!    EX=DEXP(V)
!
!    M=0
!    BBFN=0.D0
!    EM=1.D0
!5   M=M+1
!    VM=M*V
!    BM=(6.D0+VM*(6.D0+VM*(3.D0+VM)))/M**4
!    EM=EM/EX
!    BBFN=BBFN+BM*EM
!    IF(VM**3*EM.GT.EPS) GOTO 5
!    BBFN=CC*BBFN
!    RETURN
!  END FUNCTION BBFN
!
!
!!************************************************************************
!  DOUBLE COMPLEX FUNCTION W4(Z)
!! COMPUTES THE COMPLEX PROBABILITY FUNCTION W(Z)=EXP(-Z**2)*ERFC(-I*Z)
!! IN THE UPPER HALF-PLANE Z=X+I*Y (I.E. FOR Y>=0)
!! MAXIMUM RELATIVE ERROR OF BOTH REAL AND IMAGINARY PARTS IS <1*10**(-4)
!    IMPLICIT NONE
!      DOUBLE COMPLEX :: Z,T,U
!      REAL(8) :: X,Y,S
!      X=REAL(Z)
!      Y=AIMAG(Z)
!      T=CMPLX(Y,-X)
!      S=ABS(X)+Y
!      IF(S.GE.15.) THEN       !  ***  REGION I
!        W4=T*.5641896/(.5+T*T)
!      ELSEIF (S.GE.5.5) THEN  !  ***  REGION II
!        U=T*T
!        W4=T*(1.410474+U*.5641896)/(.75+U*(3.+U))
!      ELSEIF (Y.GE..195*ABS(X)-.176) THEN  !  ***  REGION III
!        W4=(16.4955+T*(20.20933+T*(11.96482+T*(3.778987+T*.5642236))))/ &
!           (16.4955+T*(38.82363+T*(39.27121+T*(21.69274+T*(6.699398+T)))))
!      ELSE                    !  ***  REGION IV
!        U=T*T
!        W4=CDEXP(U)-T*(36183.31-U*(3321.9905-U*(1540.787-U*(219.0313 &
!           -U*(35.76683-U*(1.320522-U*.56419)))))) &
!           /(32066.6-U*(24322.84-U*(9022.228-U*(2186.181-U*(364.2191 &
!           -U*(61.57037-U*(1.841439-U)))))))
!      ENDIF
!      RETURN
!    END FUNCTION W4
!
!    SUBROUTINE VOIGT(S,BL,BD,DETA,KETA)
!! COMPUTES ABSORPTION COEFFICIENT FOR A VOIGT LINE USING THE SUBROUTINE
!! BY J. HUMLICEK, J. QUANT. SPECTR. RAD. TRANSFER, V.27, PP.437-444, 1982
!!
!! INPUT:
!!   S     LINE INTENSITY (CM-2)
!!   BL    LORENTZ HALFWIDTH AT HALFHEIGHT (CM-1)
!!   BD    DOPPLER HALFWIDTH AT HALFHEIGHT (CM-1)
!!   DETA  WAVENUMBERS AWAY FROM LINE CENTER (CM-1)
!! OUTPUT:
!!   KETA  SPECTRAL ABSORPTION COEFFICIENT AT A WAVENUMBER DETA REMOVED 
!!         FROM LINE CENTER (CM-1)
!!
!    IMPLICIT NONE
!      REAL(8) :: S, BL,BD,DETA,KETA,BD2
!      DOUBLE COMPLEX :: Z
!      REAL(8),PARAMETER:: RTLN2=0.8325546               ! SQRT(ALOG(2))
!      REAL(8),PARAMETER:: RTPI=1.772454                 ! SQRT(PI)
!      IF(BD.LT.1.E-10) BD=1.E-10    ! Avoid division by zero
!! BD2 IS THE DOPPLER HALFWIDTH AT 1/e HEIGHT (CM-1)
!      BD2=BD/RTLN2                  ! CORRECTED BD
!   
!      Z=CMPLX(DETA/BD2,BL/BD2)
!      KETA=S/(RTPI*BD2)*REAL(W4(Z),8)
!    END SUBROUTINE VOIGT
!
!!-----------------------------------------------------------------------
! REAL(8) FUNCTION E1(X)
!! called by EN(n,x) 
!  IMPLICIT NONE
!  REAL(8),INTENT(in) :: X
!  REAL(8) :: A0,A1,A2,A3,A4,A5,B1,B2,B3,B4
!  IF(X .LE. 1.0d0 ) THEN
!     A0=-.57721566;A1= .99999193;A2=-.24991055
!     A3= .05519968;A4=-.00976004;A5= .00107857
!     E1=A0+X*(A1+X*(A2+X*(A3+X*(A4+X*A5))))-DLOG(X+1.D-8)
!  ELSE
!     A1= 8.5733287401;A2=18.0590169730;A3= 8.6347608925;A4=  .2677737343
!     B1= 9.5733223454;B2=25.6329561486;B3=21.0996530827;B4= 3.9584969228
!     E1=(X*(A3+X*(A2+X*(A1+X)))+A4)/(X*(B3+X*(B2+X*(B1+X)))+B4)*DEXP(-X)/X
!  ENDIF
!  RETURN
! END FUNCTION E1
!!-----------------------------------------------------------------------
! REAL(8) FUNCTION EN(N,X)
!! exponential integral of order n, eq. 13.31 in Modest book.
!  IMPLICIT NONE
!  INTEGER,INTENT(in) :: N
!  REAL(8),INTENT(in) :: X
!  INTEGER :: i
!  REAL(8) :: EX
!  EN=E1(X)
!  IF(N .GT. 1 ) THEN
!     EX=DEXP(-X)
!     DO I=2,N
!        EN=(EX-X*EN)/(I-1.D0)
!     END DO
!  ENDIF
!  RETURN
! END FUNCTION EN
!!-----------------------------------------------------------------------
!
! REAL(8) FUNCTION eb(T,wn)
!! T in K and wn is wavenumber in 1/cm
!! eb*C1 is the spectral blackbody emissive power (eq. 1.14 of Modest book)
!! eb*C1/pi is the spectral blackbody intensity, or Planck function
!
! IMPLICIT NONE
! REAL(8),INTENT(in) :: T, wn
! REAL(8),PARAMETER :: C1=3.7419d-12,C2=1.4388d0
!
! eb=wn**3/(EXP(C2/T*wn)-1.d0)
! RETURN
! END FUNCTION eb
!!-----------------------------------------------------------------------

END MODULE absc_kg_gen
