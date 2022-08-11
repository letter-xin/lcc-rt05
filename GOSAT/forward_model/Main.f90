PROGRAM mlpabscotest
use commonData
USE MLP
IMPLICIT NONE

INTEGER,PARAMETER::nstate=28!atmosphere layer number
REAL(8):: P(nstate),T(nstate),dx(nstate),xCO2(nstate)
REAL(8):: angel(2)
REAL(8):: zenith, azimuth
INTEGER:: i, j, k, q, m, s,ne,ii! loop indices
REAL(8), PARAMETER:: resolution = 0.01d0
REAL(8):: wave_bb = 6180.d0!259.3941586793d0
REAL(8):: wave_ee = 6280.0d0!7997.8258441756d0 
INTEGER:: wavesize
INTEGER:: absco_size
REAL(8), ALLOCATABLE:: wvnm(:)
REAL(8), ALLOCATABLE:: absco_wave(:)
REAL(8), ALLOCATABLE:: absco(:,:)
REAL(8), ALLOCATABLE:: opticalThickness(:)
REAL(8), ALLOCATABLE:: tran(:)
REAL(8), ALLOCATABLE:: tran_theta(:)
REAL(8), ALLOCATABLE:: rad(:)
REAL(8), ALLOCATABLE:: alb(:)
REAL(8), ALLOCATABLE:: SOO(:)
REAL(8),PARAMETER::pi=3.14159265d0
INTEGER::ierr
REAL(8)::SS(209550,3)
REAL(8)::Rayleigh(2000000,nstate)
REAL(8)::n(2000000)
REAL(8)::rn(2000000)
REAL(8)::SO(209550,2)
REAL(8)::albedo(24,2)

REAL(8)	::Ps, Ts, Ns   !consider Rayleigh

integer,parameter :: Nq=10001
real(8),dimension(Nq) :: absco_in
TYPE(GasMixInfo) :: locMix_f        ! Local mixture information


OPEN(UNIT=20, FILE= 'testdata/solar.dat', STATUS='OLD',ACTION='READ')
DO i=1,209550
    READ(20,*)(SS(i,j),j=1,3)
    SO(i,1) = SS(i,3)
    SO(i,2) = SS(i,1)*SS(i,2)/SS(i,3)/10000.d0/pi
END DO
OPEN(UNIT = 21, FILE= 'testdata/labels_T.dat', STATUS='OLD',ACTION='READ', IOSTAT = ierr)
IF(ierr /= 0) THEN
    WRITE(*,*) 'Error! Could not open output file!'
ENDIF
READ(21,*)(T(i),i=1,nstate)

OPEN(UNIT = 22, FILE= 'testdata/labels_P.dat', STATUS='OLD',ACTION='READ', IOSTAT = ierr)
IF(ierr /= 0) THEN
    WRITE(*,*) 'Error! Could not open output file!'
ENDIF
READ(22,*)(P(i),i=1,nstate)
OPEN(UNIT = 23, FILE= 'testdata/labels_x.dat', STATUS='OLD',ACTION='READ', IOSTAT = ierr)
IF(ierr /= 0) THEN
    WRITE(*,*) 'Error! Could not open output file!'
ENDIF
READ(23,*)(xCO2(i),i=1,nstate)
OPEN(UNIT = 24, FILE= 'testdata/labels_dx.dat', STATUS='OLD',ACTION='READ', IOSTAT = ierr)
IF(ierr /= 0) THEN
    WRITE(*,*) 'Error! Could not open output file!'
ENDIF
READ(24,*)(dx(i),i=1,nstate)
OPEN(UNIT = 25, FILE= 'testdata/labels_angel.dat', STATUS='OLD',ACTION='READ', IOSTAT = ierr)
IF(ierr /= 0) THEN
    WRITE(*,*) 'Error! Could not open output file!'
ENDIF
READ(25,*)(angel(j),j=1,2)
OPEN(UNIT = 26, FILE= 'testdata/labels_albedo.dat', STATUS='OLD',ACTION='READ', IOSTAT = ierr)
IF(ierr /= 0) THEN
    WRITE(*,*) 'Error! Could not open output file!'
ENDIF
DO i=1,24
    READ(26,*)(albedo(i,j),j=1,2)
END DO
OPEN(UNIT = 31, FILE= 'testdata/out_absco.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
IF(ierr /= 0) THEN
    WRITE(*,*) 'Error! Could not open output file!'
ENDIF
OPEN(UNIT = 32, FILE= 'testdata/out_trans.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
IF(ierr /= 0) THEN
    WRITE(*,*) 'Error! Could not open output file!'
ENDIF
OPEN(UNIT = 33, FILE= 'testdata/out_rad.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
IF(ierr /= 0) THEN
    WRITE(*,*) 'Error! Could not open output file!'
ENDIF
OPEN(UNIT = 34, FILE= 'testdata/out_SOO.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
IF(ierr /= 0) THEN
    WRITE(*,*) 'Error! Could not open output file!'
ENDIF

OPEN(UNIT = 50, FILE= 'testdata/out_alb.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
IF(ierr /= 0) THEN
    WRITE(*,*) 'Error! Could not open output file!'
ENDIF



write(*,*)'reslution=',resolution
wavesize = INT((wave_ee-wave_bb)/resolution + 1 + 1.E-6)
WRITE(*,*) 'waesize =' ,wavesize
ALLOCATE(wvnm(wavesize))
wvnm(1)=wave_bb
DO i=1,wavesize-1
    wvnm(i+1) = wvnm(i) + resolution
END DO

absco_size = wavesize
ALLOCATE(absco_wave(absco_size))
absco_wave(:) = wvnm(:)

!!!----------------absco_calculate----------------
WRITE(*,*) 'Calculating abscorption coefficient....'
ALLOCATE(absco(absco_size,nstate))
P(:)=P(:)/1000.d0
xCO2(:)=xCO2(:)/1000000.d0

DO i = 1,nstate
IF(P(i)-0.01d0<=0.0d0) THEN
    P(i)=0.010d0
ENDIF
END DO
DO i = 1,nstate-1
!locMix_f%P=(P(i)+P(i+1))/2
!locMix_f%T=(T(i)+T(i+1))/2
locMix_f%P=P(i)
locMix_f%T=T(i)
if (locMix_f%P<=0.065) then      
CALL get_k_mlp_low(locMix_f,absco_in(:))
else
CALL get_k_mlp_high(locMix_f,absco_in(:))
end if
absco_in(:)=exp(absco_in(:))
absco(:,i)=absco_in(:)
END DO
!locMix_f%P=(P(nstate)+1.0d0)/2
locMix_f%P=P(nstate)
locMix_f%T=T(nstate)
CALL get_k_mlp_high(locMix_f,absco_in(:))
absco_in(:)=exp(absco_in(:))
absco(:,nstate)=absco_in(:)



Do j=1,nstate
WRITE(31,"(*(E11.4,:,','))") (absco(i,j) , i= 1 , wavesize)
end do
!!!----------------transmittance----------------
WRITE(*,*) 'Calculating transmittance....'
Ns = 2.54743d19    !consider Rayleigh
Ps = 1013.25d0     !consider Rayleigh
Ts = 288.15d0      !consider Rayleigh
ALLOCATE(opticalThickness(absco_size))
ALLOCATE(tran(absco_size))
ALLOCATE(tran_theta(absco_size))
ALLOCATE(SOO(absco_size))
ALLOCATE(alb(absco_size))
ALLOCATE(rad(absco_size))
zenith=angel(1)/180*pi
azimuth=angel(2)/180*pi
DO i = 1,absco_size
    n(i) = ( 5791814.d0 / ( 238.0185d0 - (absco_wave(i)*1d-4)**2 )   +  167929.d0/(57.362d0-(absco_wave(i)*1d-4)**2) )*1d-8 + 1.d0
    rn(i) = 1.007482d-2 + 0.7990914d0/(47.48717d0 - (absco_wave(i)*1d-4)**2)
    Rayleigh(i,1) = 1.d0/(Ns* P(1)*1000/Ps *Ts/T(1)) *24.d0* pi**3 *(n(i)**2-1.d0)**2 /&
        &(10000.d0/absco_wave(i))**4/Ns**2/(n(i)**2+2.d0)**2 * (6.d0+3.d0*rn(i)) / (6.d0-7.d0*rn(i))
    opticalThickness(i)=1000d0*100d0*(dx(1)-dx(2))*(absco(i,1) + Rayleigh(i,1))*xCO2(1)  !(absco(i,1) + Rayleigh(i,1))
    DO j=2,nstate-1
        Rayleigh(i,j) =1.d0/ (Ns* P(j)*1000/Ps *Ts/T(j)) *24.d0* pi**3 *(n(i)**2-1.d0)**2 /&
            &(1d4/absco_wave(i))**4/(n(i)**2+2.d0)**2 * (6.d0+3.d0*rn(i))/(6.d0-7.d0*rn(i))
        opticalThickness(i)=opticalThickness(i)+1000d0*100d0*(dx(j)-dx(j+1))*(absco(i,j) + Rayleigh(i,j))*xCO2(j)  !(absco(i,j) + Rayleigh(i,j))
    END DO
Rayleigh(i,nstate) =1.d0/ (Ns* P(nstate)*1000/Ps *Ts/T(nstate)) *24.d0* pi**3 *(n(i)**2-1.d0)**2 /&
            &(1d4/absco_wave(i))**4/(n(i)**2+2.d0)**2 * (6.d0+3.d0*rn(i))/(6.d0-7.d0*rn(i))
    opticalThickness(i)=opticalThickness(i)+1000d0*100d0*dx(nstate)*(absco(i,nstate)+Rayleigh(i,nstate))*xCO2(nstate)
END DO
tran(:) = exp(-1.d0*opticalThickness(:))
opticalThickness(:)=opticalThickness(:)*(1/cos(zenith)+1/cos(azimuth))
tran_theta(:)=exp(-1.d0*opticalThickness(:))
k=1
ii=1
DO i = 1,absco_size
    Do while (absco_wave(i)-SO(k,1)>0.00407d0)
        k=k+1
    end do
    Do while (absco_wave(i)>=albedo(ii,1))
        ii=ii+1
    end do
    Call LinTerpDouble(albedo(ii-1,2), albedo(ii,2), alb(i), albedo(ii-1:ii,1), absco_wave(i))
    Call LinTerpDouble(SO(k,2), SO(k+1,2), SOO(i), SO(k:k+1,1), absco_wave(i))
    rad(i)=SOO(i)*tran_theta(i)*cos(zenith)*alb(i)
end do
WRITE(*,*) cos(zenith)
WRITE(*,*) cos(azimuth)

!!!----------------output----------------

!WRITE(21,*) (T(i), i=1,nstate)
!WRITE(22,*) (P(i), i=1,nstate)
!WRITE(23,*) (xCO2(i), i=1,nstate)
WRITE(32,*) (tran(i), i=1,wavesize)
!WRITE(33,*) (rad(i), i=14,wavesize)
Do i = 14,absco_size
WRITE(33,*) absco_wave(i),rad(i)
END DO
WRITE(34,*) (SOO(i), i=1,wavesize)
WRITE(50,*) (alb(i), i=1,wavesize)

CLOSE(20)
CLOSE(21)
CLOSE(22)
CLOSE(23)
CLOSE(24)
CLOSE(25)
CLOSE(26)
CLOSE(31)
CLOSE(32)
CLOSE(33)
CLOSE(34)


END PROGRAM mlpabscotest
