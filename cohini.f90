MODULE cohini
USE param
USE functions

CONTAINS

SUBROUTINE INITCOH

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:41:21
!************************************************************************

!    *MAINPROG*   COHERENS PROGRAM (MAIN PROGRAM)

!       AUTHOR - PATRICK LUYTEN, KAJO MIRBACH, KEVIN RUDDICK, CLAIRE SMITH,
!                PAUL TETT, KAREN WILD-ALLEN

!       LAST UPDATE - 1 Sep 1999 @(COHERENS)mainprog.f 8.4

!       DESCRIPTION - READ THE con-FILE WITH CONTROL PARAMETERS AND SWITCHES
!                   - OPEN FILES FOR DATA INPUT
!                   - INITIALISE THE PROGRAM
!                   - UPDATE PROGRAM VARIABLES IN THE FORM OF A TIME LOOP
!                   - WRITE FINAL CONDITIONS FILE
!                   - EXIT THE PROGRAM

!       REFERENCE - COHERENS User Documentation

!       CALLING PROGRAM -

!       EXTERNALS - ANALYS, BBCIN, BCSIN, BC2IN, BIOLGY, BSTRES, CBCIN,
!                   CONTNY, CRRNT1C, CRRNT2, CRRNT3C, CRRNT3P, DENSTY,
!                   DENSTY1, ERRMOD, HBCIN, HEAT, HEDDY, INITC, INTEGR,
!                   MASS, METIN, NEWTIM, OUTPA, OUTPAR, OUTPT, PBCIN,
!                   PRINT, RDCON, SALT, SEARHO, SEDEUL, SEDLAG, TRANSV,
!                   VEDDY1, VEDDY2, WAVIN, WAVNUM, WCALC

!************************************************************************
!     1. INITIALISATION OF THE MODEL
!     -----------------------------

CALL preproc

!     1.3 TIME-STEPS FOR INCLUSION OF VARIOUS PROCESSES
!     -------------------------------------------------

del3 = ic3d*delt
delm = icmet*delt

!     1.4 CALL INITIALISATION SUBROUTINE
!     ----------------------------------

CALL initc


!     1.6 INITIALISE BOTTOM STRESS AND DENSITY RELATED PARAMETERS
!     -----------------------------------------------------------

CALL bcsin
CALL bstres
CALL searho

do k = 1,nz
do j = 1,nr
do i = 1,nc
  uqden0(k,j,i) = 0.0
  vqden0(k,j,i) = 0.0
end do
end do
end do

do k = 1,nz
do j = 1,nr
  uqden0(k,j,nc+1) = 0.0
end do
do i = 1,nc
vqden0(k,nr+1,i) = 0.0
end do
end do

CALL densty

do k = 1,nz
do j = 1,nr
do i = 1,nc
  uqden0(k,j,i) = uqden(k,j,i)
  vqden0(k,j,i) = vqden(k,j,i)
end do
end do
end do

do k = 1,nz
do j = 1,nr
  uqden0(k,j,nc+1) = uqden(k,j,nc+1)
end do
do i = 1,nc
vqden0(k,nr+1,i) = vqden(k,nr+1,i)
end do
end do

CALL wrcon
CALL metin

END SUBROUTINE INITCOH

!********************************

SUBROUTINE preproc

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:42:37

!***********************************************************************

!    *PREPROC*     PREPROCESSOR PROGRAM (MAIN PROGRAM)

!       AUTHOR -  PATRICK LUYTEN AND KEVIN RUDDICK

!       LAST UPDATE - 1 Sep 1999       @(COHERENS)preproc.f    8.4F

!       DESCRIPTION - SET UP DEFAULT VALUES (DEFAULT, INCTUR)
!                   - DEFINE MODEL PARAMETERS (DEFCON)
!                   - READ OUTPUT SPECIFIERS (READIN)
!                   - DEFINE MODEL GRID (DEFGRID, DEFGRID1)
!                   - DEFINE 2-D OPEN BOUNDARY CONDITIONS AND DATA
!                     (DEFOB2D, WROBDAT)
!                   - DEFINE 3-D OPEN BOUNDARY CONDITIONS AND DATA
!                     (DEFOB3D, WROBDAT)
!                   - DEFINE INITIAL CONDITIONS (DEFICS)
!                   - DEFINE METEOROLOGICAL AND WAVE DATA (WRFCDAT)
!                   - WRITE INPUT FILES FOR MAIN PROGRAM: MODEL PARAMETERS
!                     (WRCON), MODEL GRID (WRGRD), OPEN BOUNDARY
!                     CONDITIONS/DATA (BCSOUT,WROBDAT), INITIAL CONDITIONS
!                     (OUTPA), METEOROLOGICAL/WAVE DATA (WRFCDAT)
!                   - VERIFIES MODEL PARAMATERS AND ARRAYS (ERRMOD)

!       REFERENCE -  COHERENS User Documentation

!       EXTERNALS - BCSOUT, DEFAULT, DEFCON, DEFGRID, DEFGRID1, DEFICS,
!                   DEFOB2D, DEFOB3D, ERRICS, ERRMOD, ERROR, INCTUR, INICON,
!                   NDIFF, OUTPA, READIN, SEARHO, STAUCH, TRANSH, TRANSV,
!                   WRCON, WRFCDAT, WRGRD, WROBDAT

!***********************************************************************
!*    LOCAL VARIABLES

REAL :: trefc(nz),srefc(nz)
LOGICAL :: erflag
INTEGER :: i, j, k, n, nsecs, kml, ktop, kbot, ncc, kint
INTEGER :: n1,n2,n3,n4,nsum,ist,ii,jj,kpos
REAL :: hc,dh,aa,xzero,yzero,bb,yb,xx,yy,ff,fac,yre
REAL :: randm,xxpos,yypos,zzpos,zdist,depp
REAL :: TINI(8000),tsum,TINI2(8000),SINI2(8000)
REAL :: SINI(8000),ssum
REAL :: CINI1(8000),c1sum, CINI2(8000),c2sum, CINI3(8000),c3sum
REAL :: CINI4(8000),c4sum
REAL :: dep2(nr,nc)

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *NSECS*    INTEGER   NUMBERS OF SECONDS BETWEEN START AND END DATE

!-----------------------------------------------------------------------

!     1. INTIALISE CONTROL PARAMETERS, WRITE .con FILE
!     ------------------------------------------------

!     ---set default values
!      WRITE (*,*) 'Setting default values'
CALL default

!     ---user-defined values for .con file
! WRITE (*,*) 'User-defined model parameters'
CALL defcon

!     ---initialise date/time
nt = 0
hour = 0.0

!     2. BATHYMETRY
!     -------------
!
open(unit=96,file='topo1.dat',form='formatted',status='unknown')
DO j = 1,nr
read(96,*)(dep(j,i),i=1,nc)
END DO
close(unit=96)

do i = 1,nc
do j = 1,nr
  dep(j,i) = min(dep(j,i),1000.0)
end do
END DO

do i = 1,nc
  dep(nr,i) = 0.0
end do

open(unit=96,file='dat/topo.dat',form='formatted',status='unknown')
DO j = 1,nr
write(96,'(250F16.8)')(dep(j,i),i=1,nc)
END DO
close(unit=96)

!     2. GENERAL SIGMA COORDINATES

! standard sigma coordinates

DO i = 1,nc
 DO j = 1,nr
  DO k = 1,nz+1
    gz0(k,j,i) = real(k-1)/real(nz)
  END DO
 END DO
END DO

! modified sigma coordinates

DO i = 1,nc
  DO j = 1,nr
    ddd = dep(j,i)
    if (ddd > 600.) then
      do k = 1,6
        gz0(k,j,i) = real(k-1)*10./ddd
      end do
      do k = 7,nz+1
        gz0(k,j,i) = gz0(6,j,i)+real(k-6)/real(55)*(ddd-50.)/ddd
      end do
     end if
  END DO
 END DO


!     ---uniform default values for roughness length

DO  i=1,nc
  DO  j=1,nr
    cdz0(j,i) = cdz0un
  END DO
END DO
!     ---universal and model constants
CALL inicon

!     2. GENERATE GRID, WRITE .grd FILE
!     ---------------------------------

!     ---grid parameters and arrays
WRITE (*,*) 'Generating model grid'
CALL defgrid

!     ---initialise grid spacing, geogr. positions and Coriolis
CALL transh

CALL transv


!     ---write space/time resolution to stdout
WRITE (*,*) 'NC=',nc
WRITE (*,*) 'NR=',nr
WRITE (*,*) 'NZ=',nz
!     5. GENERATE INITIAL CONDITIONS, WRITE .ics FILE
!     -----------------------------------------------

! *initial T and S stratification in cartesian coordinates with 0.5 m depth increment from top to bottom 
! *initial T and S stratification in cartesian coordinates with 0.5 m depth increment from top to bottom 
DO k = 1,4000
  TINI(k) = 0.0 !-REAL(k-1)/1000.0*5.0
  SINI(k) = +REAL(k-1)/1000.0*1.0
END DO

!     ---initialise salinity and temperature
DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
      kbot = int(h2atc(j,i)*(1.-gz0(k,j,i)))
      ktop = int(h2atc(j,i)*(1.-gz0(k+1,j,i)))+1
      tsum = 0.
      ssum = 0.
      ncc = 0
      do kint = ktop,kbot
       tsum = tsum+TINI(kint)
       ssum = ssum+SINI(kint)
       ncc = ncc+1
      end do
      tsum = tsum/real(ncc)
      ssum = ssum/real(ncc)
      s(k,j,i) = ssum
      ssref(k,j,i) = ssum
      tr(k,j,i) = tsum
      t(k,j,i) = tsum
      ttref(k,j,i) = tsum
      END DO
    END IF
  END DO
END DO

!     ---initialise density related parameters

  DO  i=1,nc
    DO  j=1,nr
      DO  k=1,nz
        ro(k,j,i) = 0.0
        buoy(k,j,i) = 0.0
        IF (nwd(j,i) == 1) THEN
          sbet(k,j,i) = sbetun
          tbet(k,j,i) = tbetun
        ELSE
          sbet(k,j,i) = 0.0
          tbet(k,j,i) = 0.0
        END IF
      END DO
    END DO
  END DO

! PROFILE NUMBERS ARRAYS 3D Boundary conditions
!     -------------------------
!     ---U-open boundaries
      DO II=1,NOBU
        IVPOBU(II) = II
      END DO
!    
!     ---V-open boundaries
       DO JJ=1,NOBV
        IVPOBV(JJ) = NOBU+JJ
       END DO
!

call rhoini

CALL searho

tref = 20.0
sref = 34.0

OPEN(33,file="TS.txt",form="formatted")
 DO k = 1,4000
  write(33,*)TINI(K)+tref,SINI(K)+sref
 END DO
CLOSE(33) 


!     ---intialise bottom/surface stress and turbulence arrays (default)
WRITE (*,*) 'Initialise turbulence arrays'
CALL inctur

! ---initialize euler concentration fields

DO k = 1,8000
CINI1(k) = TINI(k)
CINI2(k) = TINI(k)
CINI3(k) = 0.
CINI4(k) = 0.
END DO

!DO k = 201,300
!CINI2(k) = 1.
!END DO
!DO k = 301,400
!CINI3(k) = 1.
!END DO
!DO k = 401,1000
!CINI4(k) = 1.
!END DO

!     ---initialise concentration fields
DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
      kbot = int(2*h2atc(j,i)*(1.-gz0(k,j,i)))
      ktop = int(2*h2atc(j,i)*(1.-gz0(k+1,j,i)))+1
      c1sum = 0.
      c2sum = 0.
      c3sum = 0.
      c4sum = 0.
      ncc = 0
      do kint = ktop,kbot
       c1sum = c1sum+CINI1(kint)
       c2sum = c2sum+CINI2(kint)
       c3sum = c3sum+CINI3(kint)
       c4sum = c4sum+CINI4(kint)
       ncc = ncc+1
      end do
      c1sum = c1sum/real(ncc)
      c2sum = c2sum/real(ncc)
      c3sum = c3sum/real(ncc)
      c4sum = c4sum/real(ncc)
      conc1(k,j,i) = c1sum
      cc1ref(k,j,i) = c1sum
      conc2(k,j,i) = c2sum
      cc2ref(k,j,i) = c2sum
      conc3(k,j,i) = c3sum
      cc3ref(k,j,i) = c3sum
      conc4(k,j,i) = c4sum
      cc4ref(k,j,i) = c4sum
      END DO
    END IF
  END DO
END DO


!   ---initialise Lagrangian floats
! step 1a: set positions to zero
DO n = 1,ntrac
 xpos(n) = 0.0
 ypos(n) = 0.0
 zpos(n) = 0.0
 XT(n) = 0.0
 YT(n) = 0.0
 ZT(n) = 0.0
 IT(n) = 0
 JT(n) = 0
 KT(n) = 0
 posH(n) = 0.0
 posU(n) = 0.0
 posC(n) = 0.0
END DO

! Allocate floats JK

! initialise random function
randm = random(1)

xlen = 1
ylen = 5

n1 = 0.0

do n = 1,100*1000
  i = 45+INT(random(0)*xlen)
  j = 5+INT(random(0)*ylen)
  depp = dep(j,i)
 ! IF((depp>200.0).AND.(depp<=400.0))THEN
 ! IF(depp>400.0)THEN
  IF(n1<2000)THEN
    n1 = n1 + 1
    IT(n1) = i
    randm = random(0)-0.5
    XT(n1) = randm*gx2(1,1)
    xpos(n1) = REAL(i)*gx2(1,1)+XT(n1)-0.5*gx2(1,1)
    JT(n1) = j
    randm = random(0)-0.5
    YT(n1) = randm*gy2(1)
    ypos(n1) = REAL(j)*gy2(1)+YT(n1)-0.5*gy2(1)
    zdist = 2.*random(0)
    kpos = max(1,int(zdist))
    KT(n1) = kpos 
    randm = random(0)-0.5
    ZT(n1) = randm*gzsc(kpos,j,i)
    zpos(n1) = REAL(kpos)*gzsc(kpos,j,i)+ZT(n1)-0.5*gzsc(kpos,j,i)
!  END IF
  END IF
END DO

write(6,*)' REGION 1. No of floats allocated =', n1

n1 = 0.0

do n = 1,100*1000
  i = 45+INT(random(0)*xlen)
  j = 10+INT(random(0)*ylen)
  depp = dep(j,i)
!  IF((depp>400.0).AND.(depp<=600.0))THEN
!  IF(depp>400.0)THEN
  IF((n1<2000))THEN
    n1 = n1 + 1
    IT(n1+2000) = i
    randm = random(0)-0.5
    XT(n1+2000) = randm*gx2(1,1)
    xpos(n1+2000) = REAL(i)*gx2(1,1)+XT(n1+2000)-0.5*gx2(1,1)
    JT(n1+2000) = j
    randm = random(0)-0.5
    YT(n1+2000) = randm*gy2(1)
    ypos(n1+2000) = REAL(j)*gy2(1)+YT(n1+2000)-0.5*gy2(1)
     zdist = 2.*random(0)
    kpos = max(1,int(zdist))
    KT(n1+2000) = KPOS
    randm = random(0)-0.5
    ZT(n1+2000) = randm*gzsc(kpos,j,i)
    zpos(n1+2000) = REAL(kpos)*gzsc(kpos,j,i)+ZT(n1+2000)-0.5*gzsc(kpos,j,i)
  END IF
!  END IF
END DO

write(6,*)' REGION 2. No of floats allocated =', n1

do n = 1,ntrac
 ITini(n) = IT(n)
 JTini(n) = JT(n)
 KTini(n) = KT(n)
 XTini(n) = XT(n)
 YTini(n) = YT(n)
 ZTini(n) = ZT(n)
 zpos(n) = REAL(KT(n)-0.5)*gzsc(KT(n),JT(n),IT(n))+ZT(n)
 xposini(n) = xpos(n)
 yposini(n) = ypos(n)
 zposini(n) = zpos(n)
 depp = dep(JT(n),IT(n))
 depos(n) = (1.0-zpos(n))*depp
 trac_start(n) = .false.
 trac_stop(n) = .false.
end do

ibatch = 5

RETURN

END SUBROUTINE preproc

!****************************

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:38:38

!     *defmod.f*  INITIALISATION OF THE MODEL

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 11 Aug 1998       @(COHERENS)defmod.f 8.4

!       DESCRIPTION - SERIES OF "EXAMPLE" SUBROUTINES FOR MODEL SETUP

!       REFERENCE - Section II-2.1 of the User Documentation

!       SUBROUTINES - DEFCON, DEFGRID, DEFOB2D, DEFOB3D, DEFICS,
!                     WRFCDAT, WROBDAT

!***********************************************************************

!=======================================================================

SUBROUTINE defcon
!***********************************************************************

!    *DEFCON*

!       AUTHOR - PATRICK LUYTEN AND KEVIN RUDDICK

!       LAST UPDATE - 11 Aug 1998       @(COHERENS)defmod.f 8.4

!       DESCRIPTION - DEFINE CONTROL PARAMETERS (SWITCHES,
!                     DATE/TIME VARIABLES, MODEL PARAMETERS)

!       REFERENCE - Section II-2.1 of the User Documentation

!       CALLING PROGRAM - PREPROC

!       EXTERNALS -

!***********************************************************************
!*    LOCAL VARIABLES

REAL :: cfl

!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CFL*      REAL      CFL-NUMBER (AT THE TOP OF THE CONE)

!-----------------------------------------------------------------------


!     1. SWITCHES
!     -----------

!     ---advection scheme for momentum (0/1/2/3/4)
iadvc = 1
!     ---advection scheme for scalars (0/1/2/3/4)
iadvs = 3
!     ---type of turbulence model (0/1/2)
ioptk = 1
!     ---disabling/enabling 2-D and 3-D current calculations (0/1)
iopt2 = 1
iopt3 = 1
!     ---horizontal pressure gradient and equation of state (0/1/2)
ioptd = 1
!     ---bottom stress formulation (0/1/2)
ibstr = 1

!     3. DATE/TIME PARAMETERS
!     -----------------------

!     ---start/end date (MMDDHHMM)
ibdate = 01010000
iedate = 02410000

!     4. TIME STEP
!     ------------

!     ---time step for 2-D mode [s]

delt = 5.0 !!!2.9

!     10. BACKGROUND OR UNIFORM MIXING COEFFICIENTS
!     ---------------------------------------------

!     ---background/uniform eddy viscosity [m2/s]
vismol = 1.e-4
!     ---background/uniform eddy diffusivity [m2/s]
difmol = 1.e-4

! horizontal eddy viscosity

hedvun = 1.0
heddun = 1.0

ic3d = 20

RETURN

END SUBROUTINE defcon

!=======================================================================

SUBROUTINE bcsin
!************************************************************************

!    *BCSIN*       READ OPEN BOUNDARY CONDITIONS/DAT AT INITIAL TIME,
!                  DETERMINE ORIENTATION OF OPEN BOUNDARY FACES

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 11 May 1998       @(COHERENS)bcsin.f 8.4

!       DESCRIPTION - READ OPEN BOUNDARY ('obc') FILE
!                   - CHECK POSSIBLE ERRORS IN OPEN BOUNDARY FILES
!                   - DETERMINE ORIENTATION OF OPEN BOUNDARY FACES

!       REFERENCE -

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS - BBCIN, BC2IN, CBCIN, ERRBCS, ERROR, HBCIN, PBCIN

!************************************************************************
!------------------------------------------------------------------------


!*    LOCAL VARIABLES

LOGICAL :: erflag
CHARACTER (LEN=60) :: bcsfil
INTEGER :: i, icon, ii, iunit, j, jj, nrline

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ERFLAG*   LOGICAL   ERROR FLAG (TRUE IF UNIT IUNIT IS CONNECTED
!                         TO A FILE)
!    *BCSFIL*   CHAR      NAME OF INPUT FILE
!    *IUNIT*    INTEGER   UNIT NUMBER OF INPUT FILE
!    *NRLINE*   INTEGER   LINE NUMBER IN INPUT FILE

!------------------------------------------------------------------------
!
!     5. DETERMINE ORIENTATION OF OPEN BOUNDARY FACES
!     -----------------------------------------------
!
!     5.1 U-BOUNDARIES (WEST OR EAST)
!     -------------------------------
!
      WESTOB(0) = .FALSE.
      DO II=1,NOBU
         I = IOBU(II)
         J = JOBU(II)
         IF (I.EQ.1) THEN
            WESTOB(II) = .TRUE.
         ELSEIF (NWD(J,I-1).EQ.0) THEN
            WESTOB(II) = .TRUE.
         ELSE
            WESTOB(II)= .FALSE.
         ENDIF
       END DO

!
!     5.2 V-BOUNDARIES (SOUTH OR NORTH)
!     ---------------------------------
!

      SOUTOB(0) = .FALSE.
      DO JJ=1,NOBV
         I = IOBV(JJ)
         J = JOBV(JJ)
         IF (J.EQ.1) THEN
            SOUTOB(JJ) = .TRUE.
         ELSEIF (NWD(J-1,I).EQ.0) THEN
            SOUTOB(JJ) = .TRUE.
         ELSE
            SOUTOB(JJ) = .FALSE.
         ENDIF
      END DO

RETURN

END SUBROUTINE bcsin

!=======================================================================

SUBROUTINE defgrid
!***********************************************************************

!    *DEFGRID*    TEST CASE CONES

!       AUTHOR - PATRICK LUYTEN AND KEVIN RUDDICK

!       LAST UPDATE - 11 Aug 1998       @(COHERENS)defmod.f 8.4

!       DESCRIPTION - DEFINE GRID PARAMETERS AND ARRAYS

!       REFERENCE - Section II-2.1 of the User Documentation

!       CALLING PROGRAM - PREPROC

!       EXTERNALS -

!***********************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, j, ii, jj
REAL :: dhun


!     1. GRID COORDINATES
!     -------------------

!     ---X-direction [m or degrees]
dhun = 1000.0
DO  i=1,nc+1
  gx0(i) = dhun*(i-1)
END DO

dhun = 1000.0
!     ---Y-direction [m or degrees]
! 
DO j = 1,nr+1
   gy0(j) = dhun*(j-1)
END DO

!     3. POINTER ARRAYS
!     -----------------

!     ---cell centres (0/1)
DO  i=1,nc
  DO  j=1,nr
     IF (dep(j,i) > 0.0) THEN
      nwd(j,i) = 1
    ELSE
      nwd(j,i) = 0
    END IF
  END DO
END DO

!    ---U-faces (0/1/2/3)
      DO J=1,NR
         DO I=1,NC
            NPIX(J,I+1) = 1
            IF (NWD(J,I) == 0) THEN
               NPIX(J,I) = 0
               NPIX(J,I+1) = 0
            ENDIF
          END DO
       END DO

       DO J = 1,39
         NPIX(J,1) = 2
       END DO

       DO J = 1,39
         NPIX(J,NC+1) = 2
       END DO

!     ---V-faces (0/1/2/3)
      DO I=1,NC
         DO J=1,NR
            NPIY(J+1,I) = 1
            IF (NWD(J,I) == 0) THEN
               NPIY(J,I) = 0
               NPIY(J+1,I) = 0
            ENDIF
      END DO
         NPIY(1,I) = 2
      !   NPIY(NR+1,I) = 2
      END DO

!     4. OPEN BOUNDARY INDEX ARRAYS
!     -----------------------------
!    ---u-open boundaries
!     --west boundary
ii = 0
 DO  j=1,39
   ii = ii+1
   iobu(ii) = 1
   jobu(ii) = j
 END DO

!     --east boundary
 DO  j=1,39
   ii = ii+1
   iobu(ii) = nc+1
   jobu(ii) = j
 END DO

!    ---v-open boundaries
!     --south boundary
 ii = 0
 DO  i=1,50
   ii = ii+1
   iobv(ii) = i
   jobv(ii) = 1
 END DO

!DO  i=1,50
!   ii = ii+1
!   iobv(ii) = i
!   jobv(ii) = nr+1
! END DO

!
!  TYPE OF OPEN BOUNDARY CONDITIONS 2D
!
!     ---U-open boundaries (0/1/2/3/4)
!     --western and eastern boundary
       DO II=1,NOBU
             ITYPOBU(II) = 2
       END DO
!
!     ---V-open boundaries (0/1/2/3/4)
!     --southern boundary
       DO JJ=1,NOBV
          ITYPOBV(JJ) = 2
       END DO
!
RETURN

END SUBROUTINE defgrid

SUBROUTINE inicon
!***********************************************************************

!    *INICON*      SET UNIVERSAL AND MODEL CONSTANTS FOR PHYSICS AND TURBULENCE

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 5 Nov 1998 @(COHERENS)inicon.f 8.4

!       DESCRIPTION - PARAMETERS, SET IN INICON, CANNOT BE RESET BY THE USER
!                     IN DEFCON

!       REFERENCE - Table 5.1.2 of the User Documentation

!       CALLING PROGRAM - INITC, PREPROC

!       EXTERNALS -

!***********************************************************************
!     1. PHYSICAL CONSTANTS
!     ---------------------


ckar = 0.4
cp = 3987.5
rearth = 6371000.0


!     2. TURBULENCE
!     -------------

!     ---Mellor-Yamada
IF (itcpar == 1) THEN
  clev2(1) = 24.5
  clev2(2) = -3.26
  clev2(3) = -12.8
  clev2(4) = 65.7
  clev25(1) = 0.556
  clev25(2) = 2.18
  clev25(3) = 20.4
  clev25(4) = 53.1
  clev25(5) = 0.699
  clev25(6) = 17.3
  cmu = clev25(1)**4
  eps0 = clev25(1)**3
  sigk = 1.96
  e1 = 1.8
  e2 = e1+2.0*ckar*ckar/(sigk*SQRT(cmu))-1.0
  gamin = -0.046
  IF (ilim == 0) THEN
    gamax = 1.0E+10
  ELSE
    gamax = 0.560
  END IF
  
!     ---k-epsilon
ELSE IF (itcpar == 2) THEN
  clev2(1) = 0.6480
  clev2(2) = -0.1085
  clev2(3) = -0.02290
  clev2(4) = 0.03953
  clev25(1) = 0.108
  clev25(2) = 0.0229
  clev25(3) = 0.471
  clev25(4) = 0.0275
  clev25(5) = 0.177
  clev25(6) = 0.403
  cmu = clev25(1)
  eps0 = cmu**0.75
  sigk = 1.0
  sige = 1.3
  c2e = 1.92
  c1e = c2e-ckar*ckar/(SQRT(cmu)*sige)
  c3e1 = 0.2
  c3e2 = 1.0
  gamin = -1.725
  IF (ilim == 0) THEN
    gamax = 1.0E+11
  ELSE
    gamax = 6.903
  END IF
END IF


RETURN

END SUBROUTINE inicon
!

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:40:01

SUBROUTINE initc
!************************************************************************

!    *INITC*       INITIALISE PARAMETERS AND ARRAYS,
!                  SET UP INITIAL CONDITIONS FOR THE MODEL RUN

!       AUTHOR - PATRICK LUYTEN AND KEVIN RUDDICK

!       LAST UPDATE - 5 May 1999 @(COHERENS)initc.f 8.4

!       DESCRIPTION -
!          1) INITIALISE MODEL CONSTANTS (UNIVERSAL, TURBULENCE, BIOLOGY)
!             NOT DEFINED IN DEFCON
!          2) READ MODEL BATHYMETRY AND GRID DATA
!          3) EVALUATE GRID SIZES, CORIOLIS, ...
!          4) CALCULATE CFL LIMIT ON TIME-STEP
!          5) CALCULATE TOTAL NUMBER OF TIME-STEPS (NSTEP)
!          6) READ INITIAL ARRAYS
!          7) INITIALISE PROGRAM ARRAYS (BLOCKED CELLS/INTERFACES, SPATIALLY
!             UNIFORM ARRAYS, ZERO WORKSPACE ARRAYS, DEFAULT BOUNDARY ARRAYS)

!       REFERENCE - Section V-1.9 of the User Documentation

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS - ERRICS, INIBIO, INICON, INPA, NDIFF,
!                   RDGRD,  STAUCH, TRANSH, TRANSV, ZERO

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, ivp, j, k, n, nsecs
REAL :: cormax, depmax, dxymin

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *NSECS*    INTEGER   NUMBER OF SECONDS BETWEEN START AND END DATE
!    *DEPMAX*   REAL      MAXIMUM WATER DEPTH
!    *DXYMIN*   REAL      MINIMUM GRID SPACING

!------------------------------------------------------------------------

!     1. INITIALISE UNIVERSAL AND MODEL CONSTANTS
!     -------------------------------------------

CALL inicon

!     2.2 EVALUATION OF CELL SIZES, GEOGRAPHICAL LOCATION AND CORIOLIS FREQ.
!     ----------------------------------------------------------------------


CALL transh
CALL transv

DO  i=1,nc
  DO  j=1,nr
    DO  k=1,nz
      gz1(k,j,i) = gz2(k,j,i)
    END DO
  END DO
END DO


!     3. TIME STEP
!     ------------

!     3.1 CFL LIMIT FOR 2-D TIME STEP
!     -------------------------------

WRITE (*,*) 'Calculating CFL-limit for 2-D time step'
!     ---minimum horizontal grid spacing
dxymin = gx2(1,1)
DO  i=1, nc
  DO  j=1, nr
    IF (gx2(j,i) < dxymin) dxymin = gx2(j,i)
    IF (gy2(j) < dxymin) dxymin = gy2(j)
  END DO
END DO
WRITE (*,*) '  minimum horizontal grid spacing : ',dxymin

!     ---maximum water depth
depmax = 0.0
cormax = 0.0

DO  j=1,nr
  IF (ABS(coriol(j)) > cormax) cormax = ABS(coriol(j))
  DO  i=1,nc
    IF (dep(j,i) > depmax) depmax = dep(j,i)
  END DO
END DO
312   CONTINUE
WRITE (*,*) '  maximum water depth : ',depmax

!     ---CFL-limit
IF (cormax > 1.0E-20) THEN
  dtmax = MIN(1.0/cormax,0.5*dxymin/SQRT(g*depmax))
ELSE
  dtmax = 0.5*dxymin/SQRT(g*depmax)
END IF
WRITE (*,*) '  CFL-limit : ',dtmax
WRITE (*,*) '  2-D time step : ',delt


!     3.2 CALCULATE TOTAL NUMBER OF TIME-STEPS (NSTEP)
!     ------------------------------------------------
!     5. INITIALISE ARRAYS ON DRY CELLS/FACES
!     ---------------------------------------

!     5.1 CELL CENTRES
!     ----------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 0) THEN
      
!           ---2-D arrays
      dep(j,i) = 0.0
      zeta2(j,i) = 0.0
      atcf1(j,i) = 0.0
      atcf2(j,i) = 0.0
      r1opt(j,i) = 0.0
      r2opt(j,i) = 0.0
      hexp(j,i) = 0.0
      cdz0(j,i) = 0.0
      fs(j,i) = 0.0
      gs(j,i) = 0.0
      sstot(j,i) = 0.0
!           ---3-D arrays
      DO  k=1,nz
        ro(k,j,i) = 0.0
        roref(k,j,i) = 0.0
        sbet(k,j,i) = 0.0
        tbet(k,j,i) = 0.0
        w2phys(k,j,i) = 0.0
        heddyvc(k,j,i) = 0.0
        heddydc(k,j,i) = 0.0
      END DO
    END IF
  END DO
END DO



!     5.2 U-NODES
!     -----------


DO  i=1,nc+1
  DO  j=1,nr
    IF (npix(j,i) == 0) THEN
!           ---2-D arrays
      ud2(j,i) = 0.0
!           ---3-D arrays
      DO  k=1,nz
        u2(k,j,i) = 0.0
        u2f(k,j,i) = 0.0
        heddyvu(k,j,i) = 0.0
        heddydu(k,j,i) = 0.0
      END DO
    END IF
  END DO
END DO


!     5.3 V-NODES
!     -----------


DO  i=1,nc
  DO  j=1,nr+1
    IF (npiy(j,i) == 0) THEN
!           ---2-D arrays
      vd2(j,i) = 0.0
!           ---3-D arrays
      DO  k=1,nz
        v2(k,j,i) = 0.0
        v2f(k,j,i) = 0.0
        heddyvv(k,j,i) = 0.0
        heddydv(k,j,i) = 0.0
      END DO
    END IF
  END DO
END DO


!     5.4 W-NODES
!     -----------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 0) THEN
!           ---3-D arrays
      DO  k=1,nz+1
        w2(k,j,i) = 0.0
        smu(k,j,i) = 0.0
        shb(k,j,i) = 0.0
        dissw(k,j,i) = 0.0
        tkew(k,j,i) = 0.0
        zlw(k,j,i) = 0.0
        veddyv(k,j,i) = 0.0
        veddyd(k,j,i) = 0.0
        veddyk(k,j,i) = 0.0
        veddye(k,j,i) = 0.0
      END DO
    END IF
  END DO
END DO


!     6. INITIALISE SPATIALLY UNIFORM VALUES
!     --------------------------------------

!     6.1 CELL CENTRES
!     ----------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      
!           ---2-D arrays
      atcf1(j,i) = atcf1un
      atcf2(j,i) = atcf2un
      r1opt(j,i) = r1optun
      r2opt(j,i) = r2optun
      hexp(j,i) = hexpun
      
      IF (ioptm == 1) THEN
        fs(j,i) = fsun
        gs(j,i) = gsun
        sstot(j,i) = SQRT(fsun**2+gsun**2)
      ELSE
        fs(j,i) = 0.0
        gs(j,i) = 0.0
        sstot(j,i) = 0.0
      END IF
!           ---3-D arrays
      DO  k=1,nz

        IF (ioptd == 0) THEN
          ro(k,j,i) = r0ref
          sbet(k,j,i) = 0.0
          tbet(k,j,i) = 0.0
        ELSE IF (ioptd == 1) THEN
          sbet(k,j,i) = sbetun
          tbet(k,j,i) = tbetun
        END IF
        IF (iodif == 1) THEN
          heddyvc(k,j,i) = hedvun
          heddydc(k,j,i) = heddun
        END IF
      END DO
    END IF
  END DO
END DO


!     6.2 U-NODES
!     -----------


DO  i=1,nc+1
  DO  j=1,nr
    IF (npix(j,i) > 0) THEN
!           ---3-D arrays
      DO  k=1,nz
        u2f(k,j,i) = u2(k,j,i)
        IF (iodif == 1) THEN
          heddyvu(k,j,i) = hedvun
          heddydu(k,j,i) = heddun
        END IF
      END DO
    END IF
  END DO
END DO


!     6.3 V-NODES
!     -----------


DO  i=1,nc
  DO  j=1,nr+1
    IF (npiy(j,i) > 0) THEN
!           ---3-D arrays
      DO  k=1,nz
        v2f(k,j,i) = v2(k,j,i)
        IF (iodif == 1) THEN
          heddyvv(k,j,i) = hedvun
          heddydv(k,j,i) = heddun
        END IF
      END DO
    END IF
  END DO
END DO



!     6.4 W-NODES
!     -----------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
!           ---3-D arrays
      veddyv(1,j,i) = vismol
      veddyd(1,j,i) = difmol
      veddyk(1,j,i) = vismol
      veddye(1,j,i) = vismol
      IF (ioptk == 0) THEN
        DO  k=1,nz+1
          veddyv(k,j,i) = vismol
          veddyd(k,j,i) = difmol
        END DO
      END IF
      veddyv(nz+1,j,i) = vismol
      veddyd(nz+1,j,i) = difmol
      veddyk(nz+1,j,i) = vismol
      veddye(nz+1,j,i) = vismol
    END IF
  END DO
END DO

!     7.0 ZERO INITIAL ARRAYS
!     -----------------------

!     7.1 CELL CENTRES
!     ----------------

!     ---2-D arrays
CALL zero33(p2,1,nr,nc)
CALL zero33(windu2,1,nr,nc)
CALL zero33(windv2,1,nr,nc)
CALL zero33(sst2,1,nr,nc)
CALL zero33(sat2,1,nr,nc)
CALL zero33(hum2,1,nr,nc)
CALL zero33(cloud2,1,nr,nc)
CALL zero33(rain2,1,nr,nc)
CALL zero33(bstot,1,nr,nc)
CALL zero33(cdb,1,nr,nc)
CALL zero33(cdb100,1,nr,nc)
CALL zero33(qnsol,1,nr,nc)
CALL zero33(qsol,1,nr,nc)
CALL zero33(evapr,1,nr,nc)
CALL zero33(ssalfl,1,nr,nc)
CALL zero33(dheddyvc,1,nr,nc)
!     ---3-D arrays
CALL zero33(atcfcor2,nz,nr,nc)
CALL zero33(par,nz,nr,nc)
CALL zero33(qheat,nz,nr,nc)
CALL zero33(buoy,nz,nr,nc)
CALL zero33(w2phys,nz,nr,nc)
IF (iodif /= 1) THEN
  CALL zero33(heddyvc,nz,nr,nc)
  CALL zero33(heddydc,nz,nr,nc)
END IF

!     7.2 U-NODES
!     -----------

!     ---2-D arrays
CALL zero33(uadhdev,1,nr,nc+1)
CALL zero33(uah2d,1,nr,nc+1)
CALL zero33(udh2d,1,nr,nc+1)
CALL zero33(udp,1,nr,nc+1)
CALL zero33(udqden,1,nr,nc+1)
CALL zero33(fb,1,nr,nc+1)
CALL zero33(fbk,1,nr,nc+1)
CALL zero33(dheddyvu,1,nr,nc+1)
!     ---3-D arrays
CALL zero33(uqden,nz,nr,nc+1)
IF (iodif /= 1) THEN
  CALL zero33(heddyvu,nz,nr,nc+1)
  CALL zero33(heddydu,nz,nr,nc+1)
END IF


!     7.3 V-NODES
!     -----------

!     ---2-D arrays
CALL zero33(vadhdev,1,nr+1,nc)
CALL zero33(vah2d,1,nr+1,nc)
CALL zero33(vdh2d,1,nr+1,nc)
CALL zero33(vdp,1,nr+1,nc)
CALL zero33(vdqden,1,nr+1,nc)
CALL zero33(gb,1,nr+1,nc)
CALL zero33(gbk,1,nr+1,nc)
CALL zero33(dheddyvv,1,nr+1,nc)

!     ---3-D arrays
CALL zero33(vqden,nz,nr+1,nc)
IF (iodif /= 1) THEN
  CALL zero33(heddyvv,nz,nr+1,nc)
  CALL zero33(heddydv,nz,nr+1,nc)
END IF


!     7.4 W-NODES
!     -----------

!     ---3-D arrays
CALL zero33(buprod,nz+1,nr,nc)
CALL zero33(shprod,nz+1,nr,nc)


!     7.5 WORKSPACE ARRAYS
!     --------------------


CALL zero33(zeros1,1,nr,nc)
CALL zero33(zeros2,1,nr,nc)
CALL zero33(zeros3,1,nr,nc)
CALL zero33(zeros4,1,nr,nc)
CALL zero33(zeros5,1,nr,nc)
CALL zero33(zeros6,1,nr,nc)



!     8. "QVP"-ARRAYS FOR TURBULENCE
!     ------------------------------------------

DO  ivp=0,nvprof
  turvp(0,ivp) = -99.0
  DO  k=1,nz
    turvp(k,ivp) = -100.0
  END DO
  turvp(nz+1,ivp) = -100.0
END DO



RETURN

END SUBROUTINE initc

!=======================================================================

SUBROUTINE default
!***********************************************************************

!    *DEFAULT*    DEFAULT VALUES FOR MODEL SETUP

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 5 May 1999       @(COHERENS)preproc.f    8.4

!       DESCRIPTION - SET ALL MODEL SWITCHES AND PARAMETERS BY DEFAULT
!                   - ALL OPEN BOUNDARY ARRAYS ARE SET TO A ZERO GRADIENT
!                     CONDITIONS
!                   - ALL OTHER ARRAYS (USED FOR MODEL INPUT) ARE SET TO ZERO

!       REFERENCE - Table 4.2.2 (switches) and 4.2.3 (model constants) of the
!                   User Documentation

!       CALLING PROGRAM - PREPROC

!       EXTERNALS - IZERO, ZERO

!***********************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, icon, ivp, j, k, n



!     1. CONTROL PARAMETERS
!     ---------------------

!     1.1 SWITCHES
!     ------------


!     ---advection scheme (upwind)
iadvc = 1
iadvs = 1
iadvwb = 1
!     ---horizontal diffusion (2==SMAGORINSKY)
iodif = 2
!     ---bottom stress (quadratic)
ibstr = 1
!     ---surface drag coefficient (constant)
idrag = 0
!     ---grid dimension (3-D)
igrdim = 3
!     ---grid (Cartesian)
igtrh = 0
!     ---equation of state (linear)
ioptd = 1
!     ---temperature (on)
iopthe = 1
!     ---turbulence scheme (on)
ioptk = 1
!     ---met forcing (on)
ioptm = 1
!     ---salinity (on)
ioptsa = 1
!     ---currents (on)
iopt2 = 1
iopt3 = 1
!     ---heat exchange coefficient (constant)
itdif = 0

!     ---turbulence switches
itform = 1
ntrans = 2
itcpar = 2
istpar = 1
ileng = 4
ilim = 1
iahdht = 1

!     1.2 DATE/TIME PARAMETERS
!     ------------------------

!     ---counters
ic3d  = 10
ic2bc = 0
ichbc = 0
icbbc = 0
icpbc = 0
iccbc = 0
icmet = 1
!     ---dates
ibyear = 1998
ieyear = 1998


maxfiles = 256


!     1.4 MODEL PARAMETERS
!     --------------------

!     1.4.1 PHYSICAL MODEL
!     --------------------

!     ---grid parameters
dlaref = -35.
dloref = 135.0
depun = 500.0
!     ---reference values
r0ref = 1026.0
sref = 34.0
tref = 20.0
!     ---expansion coefficients
sbetun = 7.6E-04
tbetun = 1.8E-04
!     ---optical parameters
atcf1un = 10.0
atcf2un = 2.06
r1optun = 0.54
r2optun = 0.4
epssal = 0.05714
hexpun = 7.0
!     ---linear bottom friction coefficient
cdlin = 0.e-4
!     ---bottom roughness length
cdz0un = 2.e-3
!     ---background viscosity, diffusivity
vismol = 1.e-5
difmol = 1.e-5
!     ---parameters for horizontal diffusion
hedvun = 1.0
heddun = 1.0
cm0 = 0.1
cs0 = 0.1
!     ---surface stress (uniform)
fsun = 0.0
gsun = 0.0

!     1.4.2 TURBULENCE PARAMETERS
!     ---------------------------

!     ---Pacanowski-Philander
ppa = 5.0
ppn = 2.0
ppved0 = 0.01
ppvisbg = 1.0E-04
ppdifbg = 1.0E-05
rmaxnut = 3.0
!     ---Munk-Anderson
amved0 = 0.06
ama = 10.0
amb = 3.33
amn1 = 0.5
amn2 = 1.5
rmaxlat = 4.0
!     ---flow-dependent formulations
adk1 = 0.0025
adk2 = 2.0E-05
adcnu = 2.0
adlam = 0.0
add1 = 0.0
add2 = 0.0
adr1 = 1.0
adr2 = 1.0
!     ---mixing length formulations
ablac = 0.2
bxing = -2.0
!     ---roughness lengths
z0bot = 0.0
z0sur = 0.0
!     ---limiting conditions
tkemin = 1.0E-06



!     2. GRID PARAMETERS
!     ------------------

!     2.1 REGULAR SIGMA-GRID IN THE VERTICAL
!     --------------------------------------
!     2.2 GRID ARRAYS
!     ---------------


CALL zero33 (gx0,1,1,nc+1)
CALL zero33 (gy0,1,nr+1,1)
CALL zero33 (dep,1,nr,nc)
CALL izero33 (nwd,1,nr,nc)
CALL izero33 (npix,1,nr,nc+1)
CALL izero33 (npiy,1,nr+1,nc)
CALL zero33 (cdz0,1,nr,nc)


!     3. OPEN BOUNDARY CONDITIONS
!     ---------------------------

!     3.1 U-OPEN BOUNDARIES
!     ---------------------

CALL izero33(iobu(0),nobu+1,1,1)
CALL izero33(jobu(0),nobu+1,1,1)
CALL izero33(itypobu(0),nobu+1,1,1)
CALL izero33(ivpobu(0),nobu+1,1,1)
CALL izero33(lstobu(0),nobu+1,1,1)

CALL zero33(r1obu(0),nobu+1,1,1)
CALL zero33(r2obu(0),nobu+1,1,1)
CALL zero33(ampobu(0,0),nobu+1,ncon+1,1)
CALL zero33(phaobu(0,0),nobu+1,ncon+1,1)
CALL zero33(qstobu(0),nobu+1,1,1)
CALL zero33(rmassu(0,1),nobu+1,nz,1)


!     3.2 V-OPEN BOUNDARIES
!     ---------------------


CALL izero33(iobv(0),nobv+1,1,1)
CALL izero33(jobv(0),nobv+1,1,1)
CALL izero33(itypobv(0),nobv+1,1,1)
CALL izero33(ivpobv(0),nobv+1,1,1)
CALL izero33(lstobv(0),nobv+1,1,1)
CALL zero33(r1obv(0),nobv+1,1,1)
CALL zero33(r2obv(0),nobv+1,1,1)
CALL zero33(ampobv(0,0),nobv+1,ncon+1,1)
CALL zero33(phaobv(0,0),nobv+1,ncon+1,1)
CALL zero33(qstobv(0),nobv+1,1,1)
CALL zero33(rmassv(0,1),nobv+1,nz,1)


!     3.3 "QVP"-ARRAYS
!     ----------------

DO  ivp=0,nvprof
  uvp(0,ivp) = -99.0
  svp(0,ivp) = -99.0
  tvp(0,ivp) = -99.0
  cvp(0,ivp) = -99.0
  p2bvp(0,ivp) = -99.0
  p2cvp(0,ivp) = -99.0
  p2mvp(0,ivp) = -99.0
  p2nvp(0,ivp) = -99.0
  p2nhsvp(0,ivp) = -99.0
  p2nosvp(0,ivp) = -99.0
  sedc1vp(0,ivp) = -99.0
  DO  k=1,nz
    uvp(k,ivp) = -100.0
    svp(k,ivp) = -100.0
    tvp(k,ivp) = -100.0
    cvp(k,ivp) = -100.0
    p2bvp(k,ivp) = -100.0
    p2cvp(k,ivp) = -100.0
    p2mvp(k,ivp) = -100.0
    p2nvp(k,ivp) = -100.0
    p2nhsvp(k,ivp) = -100.0
    p2nosvp(k,ivp) = -100.0
    sedc1vp(k,ivp) = -100.0
  END DO
END DO
331   CONTINUE

DO  n=0,nconc
  DO  ivp=0,nvprof
    convp(0,ivp,n) = -99.0
    DO  k=1,nz
      convp(k,ivp,n) = -100.0
    END DO
  END DO
END DO
332   CONTINUE

!     4. INITIAL CONDITIONS
!     ---------------------

!     ---currents
CALL zero33 (u2,nz,nr,nc+1)
CALL zero33(v2,nz,nr+1,nc)
CALL zero33(w2,nz+1,nr,nc)
CALL zero33(u2f,nz,nr,nc+1)
CALL zero33(v2f,nz,nr+1,nc)
CALL zero33(ud2,1,nr,nc+1)
CALL zero33(vd2,1,nr+1,nc)

!     ---bottom stress
CALL zero33(fb,1,nr,nc+1)
CALL zero33(gb,1,nr+1,nc)

!     ---surface elevation
CALL zero33(zeta2,1,nr,nc)

!     ---salinity, temperature
CALL zero33(s,nz,nr,nc)
CALL zero33(t,nz,nr,nc)
CALL zero33(conc1,nz,nr,nc)
CALL zero33(conc2,nz,nr,nc)

!     ---turbulence arrays
CALL zero33(veddyv,nz+1,nr,nc)
CALL zero33(veddyd,nz+1,nr,nc)
CALL zero33(veddyk,nz+1,nr,nc)
CALL zero33(veddye,nz+1,nr,nc)
CALL zero33(smu,nz+1,nr,nc)
CALL zero33(shb,nz+1,nr,nc)
CALL zero33(tkew,nz+1,nr,nc)
CALL zero33(zlw,nz+1,nr,nc)
CALL zero33(dissw,nz+1,nr,nc)

!     ---contaminants arrays (Eulerian module)
DO  n=0,nconc
  DO  i=1,nc
    DO  j=1,nr
      DO  k=1,nz
        concn1(k,j,i,n) = 0.0
      END DO
    END DO
  END DO
END DO

!     5. TIDAL PHASES AND FREQUENCIES
!     -------------------------------


DO  icon=0,ncon
  phase0(icon) = 0.0
  sigma(icon)  = 0.0
END DO


!     6. METEOROLOGICAL/WAVE DATA
!     ---------------------------

CALL zero33(windu2,1,nr,nc)
CALL zero33(windv2,1,nr,nc)
CALL zero33(p2,1,nr,nc)
CALL zero33(sat2,1,nr,nc)
CALL zero33(hum2,1,nr,nc)
CALL zero33(cloud2,1,nr,nc)
CALL zero33(evapr,1,nr,nc)
CALL zero33(ssalfl,1,nr,nc)
CALL zero33(rain2,1,nr,nc)
CALL zero33(qnsol,1,nr,nc)
CALL zero33(qsol,1,nr,nc)

RETURN

END SUBROUTINE default

!=======================================================================
!C=======================================================================

SUBROUTINE inctur
!***********************************************************************

!    *INCTUR*     INTIALISES TURBULENCE ARRAYS BY DEFAULT

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 21 Apr 1998       @(COHERENS)preproc.f    8.4

!       DESCRIPTION - NO INITIALISATION IS PERFORMED IF TURBULENCE IS
!                     ALREADY DEFINED IN DEFICS
!                   - INITIALISES TURBULENCE ARRAYS BY DEFAULT
!                     (CONSTANT VALUES IF IOPTK=0)
!                   - T.K.E. VARIES LINEARLY BETWEEN BOTTOM AND SURFACE
!                     VALUE
!                   - TURBULENCE LENGTH SCALE AND DISSIPATION RATE ARE
!                     CALCULATED IN SUBROUTINE TLENG USING SWITCH ILENG

!       REFERENCE - Section III-1.6.7 of the User Documentation

!       CALLING PROGRAM - PREPROC

!       EXTERNALS - BSTRES, BUOFR2, DISLEN, FNSHB, FNSMU, METIN, TLENG

!***********************************************************************
!*    LOCAL VARIABLES

LOGICAL :: flag
INTEGER :: i, iioptw, j, k
REAL :: ghan, tkewmin, z1
!REAL :: buofr2, fnshb, fnsmu
PARAMETER (tkewmin = 1.0E-06)

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *FLAG*     LOGICAL   FLAG TO DETERMINE WHETHER TURBULENCE ARRAYS ARE
!                         ALREADY INITIALISED
!    *TKEWMIN*  REAL      MINIMUM VALUE FOR T.K.E. [m2/s2]
!    *IILIM*    INTEGER   USED TO STORE VALUE OF ILIM SINCE TURBULENCE ARRAYS
!                         ARE INITIALISED WITHOUT LIMITING CONDITION
!    *GHAN*     REAL      = G_h (MELLOR-YAMADA), = ALPHA_N (K-EPS)

!------------------------------------------------------------------------

!     1. INITIALISE SURFACE AND BOTTOM STRESS
!     ---------------------------------------

!     ---check if turbulence is already initialised
flag = .false.
DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        IF (.NOT.flag) THEN
          IF (veddyv(k,j,i) > 0.0) flag = .true.
        END IF
      END DO
    END IF
  END DO
END DO
IF (flag) GO TO 1000
!     ---surface stress
IF (ioptm > 1) CALL metin
!     ---bottom stress
iioptw = ioptw
ioptw = 0
CALL bstres
ioptw = iioptw


!     2. UNIFORM EDDY COEFFICIENTS (IOPTK=0)
!     --------------------------------------

DO  i=1,nc
  DO  j=1,nr
    DO  k=1,nz
      IF (nwd(j,i) == 1) THEN
        veddyv(k,j,i) = vismol
        veddyd(k,j,i) = difmol
      ELSE
        veddyv(k,j,i) = 0.0
        veddyd(k,j,i) = 0.0
      END IF
    END DO
  END DO
END DO

IF (ioptk == 0) GO TO 1000


!     3. TURBULENCE ARRAYS
!     --------------------

!     ---turbulence kinetic energy
DO  i=1,nc
  DO  j=1,nr
    DO  k=1,nz+1
      IF (nwd(j,i) == 1) THEN
        tkew(k,j,i) = ((sstot(j,i)-bstot(j,i))*gz0(k,j,i) + bstot(j,i))/SQRT(cmu)
        tkew(k,j,i) = MAX(tkew(k,j,i),tkewmin)
      ELSE
        tkew(k,j,i) = 0.0
      END IF
    END DO
  END DO
END DO

!     ---mixing length, dissipation rate
CALL tleng
CALL dislen


!     4. EVALUATE TURBULENT DIFFUSION COEFFICIENTS
!     --------------------------------------------

!     4.1 K-L MODEL
!     -------------

IF (itcpar == 1) THEN
  
  DO  i=1,nc
    DO  j=1,nr
      DO  k=2,nz
        IF (nwd(j,i) == 1) THEN
          IF (tkew(k,j,i) > 0.0) THEN
            ghan = buofr2(k,j,i)*(zlw(k,j,i)/tkew(k,j,i))*zlw(k,j,i)
            IF (ghan > gamax) THEN
              ghan = gamax
              z1 = ghan*tkew(k,j,i)**1.5/(zlw(k,j,i)*buofr2(k,j,i))
            ELSE IF (ghan < gamin) THEN
              ghan = gamin
              z1 = ghan*tkew(k,j,i)**1.5/(zlw(k,j,i)*buofr2(k,j,i))
            ELSE
              z1 = SQRT(tkew(k,j,i))*zlw(k,j,i)
            END IF
            smu(k,j,i) = fnsmu(k,j,i,ghan)
            shb(k,j,i) = fnshb(k,j,i,ghan)
            z1 = SQRT(tkew(k,j,i))*zlw(k,j,i)
            veddyv(k,j,i) = smu(k,j,i)*z1
            veddyd(k,j,i) = shb(k,j,i)*z1
          ELSE
            smu(k,j,i) = 0.0
            shb(k,j,i) = 0.0
            veddyv(k,j,i) = 0.0
            veddyd(k,j,i) = 0.0
          END IF
          veddyk(k,j,i) = veddyv(k,j,i)/sigk + vismol
          veddyv(k,j,i) = veddyv(k,j,i) + vismol
          veddyd(k,j,i) = veddyd(k,j,i) + difmol
        END IF
      END DO
    END DO
  END DO
  
  
!     4.2 K-EPS MODEL
!     ---------------
  
ELSE IF (itcpar == 2) THEN
  
  DO  i=1,nc
    DO  j=1,nr
      DO  k=2,nz
        IF (nwd(j,i) == 1) THEN
          IF (tkew(k,j,i) > 0.0) THEN
            ghan = buofr2(k,j,i)*(tkew(k,j,i)/dissw(k,j,i))**2
            IF (ghan > gamax) THEN
              ghan = gamax
              z1 = ghan*dissw(k,j,i)/buofr2(k,j,i)

            ELSE IF (ghan < gamin) THEN
              ghan = gamin
              z1 = ghan*dissw(k,j,i)/buofr2(k,j,i)
            ELSE
              z1 = (tkew(k,j,i)/dissw(k,j,i))*tkew(k,j,i)
            END IF
            smu(k,j,i) = fnsmu(k,j,i,ghan)
            shb(k,j,i) = fnshb(k,j,i,ghan)
            veddyv(k,j,i) = smu(k,j,i)*z1
            veddyd(k,j,i) = shb(k,j,i)*z1
          ELSE
            smu(k,j,i) = 0.0
            shb(k,j,i) = 0.0
            veddyv(k,j,i) = 0.0
            veddyd(k,j,i) = 0.0
          END IF
          veddyk(k,j,i) = veddyv(k,j,i)/sigk + vismol
          veddye(k,j,i) = veddyv(k,j,i)/sige + vismol
          veddyv(k,j,i) = veddyv(k,j,i) + vismol
          veddyd(k,j,i) = veddyd(k,j,i) + difmol
        END IF
      END DO
    END DO
  END DO
  
END IF


1000 CONTINUE


RETURN

END SUBROUTINE inctur

!=======================================================================

SUBROUTINE wrcon
!***********************************************************************

!    *WRCON*      WRITE THE 'con'-FILE WITH MODEL SWITCHES AND PARAMETERS

!       AUTHOR - PATRICK LUYTEN AND KEVIN RUDDICK

!       LAST UPDATE - 5 May 1999       @(COHERENS)preproc.f    8.4

!       DESCRIPTION - WRITE THE '.con'-FILE (PARAMETERS DEFINED IN DEFCON,
!                     DEFAULT OR 'defmod.par')
!                   - WRITE ESTIMATED STORAGE REQUIREMENTS TO STDOUT

!       REFERENCE -

!       CALLING PROGRAM - PREPROC

!       EXTERNALS - OPENF

!***********************************************************************
!*    LOCAL VARIABLES

CHARACTER (LEN=11) :: iname
INTEGER :: nnodes(4)
INTEGER :: iarr, icon, ires, ivar, n, nbytes2, nbytes3
INTEGER :: nbytper, nbyttot, ncontmax, ninfi, nofiles, noufi

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *INAME*    CHAR      TEMPORARY NAME FOR OUTPUT FILES
!    *NNODES*   INTEGER   NO. OF VALUES IN X-, Y-, Z- DIRECTION AND TIME
!                         WRITTEN TO EACH OUTPUT FILE
!    *NBYTES2*  INTEGER   ESTIMATED NO. OF BYTES STORED IN 2-D RESULTS FILE
!    *NBYTES3*  INTEGER   ESTIMATED NO. OF BYTES STORED IN 3-D RESULTS FILE
!    *NBYTPER*  INTEGER   NO. OF BYTES STORED PER VALUE IN OUTPUT FILES
!    *NBYTTOT*  INTEGER   ESTIMATED TOTAL NO. OF BYTES STORED IN OUPUT FILES
!    *NINFI*    INTEGER   NO. OF OUTPUT FILES FOR INITIALISING THE MODEL
!    *NOFILES*  INTEGER   NO. OF MODEL OUTPUT FILES (TIME SERIES,
!                         HARMONIC ANALYSIS, TIME-AVERAGED OUTPUT)
!    *NOUFI*    INTEGER   NO. OF FINAL OUPUT FILES

!-----------------------------------------------------------------------

!     1. OPEN .con FILE
!     -----------------

OPEN(20,FILE='header.txt',FORM='formatted', STATUS='unknown')


!     2. SIMULATION TITLE
!     -------------------

WRITE (20,9001) 'Title of the simulation'
WRITE (20,9002) 'TITLE'


!     3. CONTROL PARAMETERS
!     ---------------------

WRITE (20,9001) 'Control parameters'


!     3.1 SWITCHES
!     ------------

WRITE (20,9003) 'Switches :'
WRITE (20,9004) 'IADVC', iadvc
WRITE (20,9004) 'IADVS', iadvs
WRITE (20,9004) 'IADVWB',iadvwb
WRITE (20,9004) 'IBSTR', ibstr
WRITE (20,9004) 'IDRAG', idrag
WRITE (20,9004) 'IFLUFF',ifluff
WRITE (20,9004) 'IGRDIM',igrdim
WRITE (20,9004) 'IGTRH', igtrh
WRITE (20,9004) 'IODIF', iodif
WRITE (20,9004) 'IOPTB', ioptb
WRITE (20,9004) 'IOPTC', ioptc
WRITE (20,9004) 'IOPTD', ioptd
WRITE (20,9004) 'IOPTHE',iopthe
WRITE (20,9004) 'IOPTK', ioptk
WRITE (20,9004) 'IOPTM', ioptm
WRITE (20,9004) 'IOPTP', ioptp
WRITE (20,9004) 'IOPTS', iopts
WRITE (20,9004) 'IOPTSA',ioptsa
WRITE (20,9004) 'IOPTW', ioptw
WRITE (20,9004) 'IOPT2', iopt2
WRITE (20,9004) 'IOPT3', iopt3
WRITE (20,9004) 'IOUTF', ioutf
WRITE (20,9004) 'IOUTP', ioutp
WRITE (20,9004) 'IOUTS', iouts
WRITE (20,9004) 'IOUTT', ioutt
WRITE (20,9004) 'ITDIF', itdif

WRITE (20,9003) 'Turbulence switches :'
WRITE (20,9004) 'ITFORM', itform
WRITE (20,9004) 'NTRANS', ntrans
WRITE (20,9004) 'ITCPAR', itcpar
WRITE (20,9004) 'ISTPAR', istpar
WRITE (20,9004) 'ILENG',  ileng
WRITE (20,9004) 'ILIM',   ilim
WRITE (20,9004) 'IAHDHT', iahdht


!     3.2 DATE/TIME PARAMETERS
!     ------------------------

WRITE (20,9003) '2D time step :'
WRITE (20,9005) 'DELT', delt

WRITE (20,9003) 'Counters :'
WRITE (20,9004) 'IC3D',   ic3d
WRITE (20,9004) 'IC2BC',  ic2bc
WRITE (20,9004) 'ICHBC',  ichbc
WRITE (20,9004) 'ICBBC',  icbbc
WRITE (20,9004) 'ICPBC',  icpbc
WRITE (20,9004) 'ICCBC',  iccbc
WRITE (20,9004) 'ICMET',  icmet
WRITE (20,9004) 'ICWAV',  icwav

WRITE (20,9003) 'Dates :'
WRITE (20,9006) 'IBDATE', ibdate
WRITE (20,9006) 'IEDATE', iedate
WRITE (20,9004) 'IBYEAR', ibyear
WRITE (20,9004) 'IEYEAR', ieyear


!     5. MODEL PARAMETERS
!     --------------------


WRITE (20,9001) 'Model parameters'


!     5.1 PHYSICAL MODEL
!     ------------------

WRITE (20,9003) 'Physical model :'
WRITE (20,9005) 'DLAREF',  dlaref
WRITE (20,9005) 'DLOREF',  dloref
WRITE (20,9005) 'DEPUN',   depun
WRITE (20,9005) 'R0REF',   r0ref
WRITE (20,9005) 'SREF',    sref
WRITE (20,9005) 'TREF',    tref
WRITE (20,9005) 'SBETUN',  sbetun
WRITE (20,9005) 'TBETUN',  tbetun
WRITE (20,9005) 'ATCF1UN', atcf1un
WRITE (20,9005) 'ATCF2UN', atcf2un
WRITE (20,9005) 'R1OPTUN', r1optun
WRITE (20,9005) 'R2OPTUN', r2optun
WRITE (20,9005) 'EPSSAL',  epssal
WRITE (20,9005) 'HEXPUN',  hexpun
WRITE (20,9005) 'CDLIN',   cdlin
WRITE (20,9005) 'VISMOL',  vismol
WRITE (20,9005) 'DIFMOL',  difmol
WRITE (20,9005) 'HEDVUN',  hedvun
WRITE (20,9005) 'HEDDUN',  heddun
WRITE (20,9005) 'CM0',     cm0
WRITE (20,9005) 'CS0',     cs0
WRITE (20,9005) 'FSUN',    fsun
WRITE (20,9005) 'GSUN',    gsun
WRITE (20,9005) 'HSUN',    hsun
WRITE (20,9005) 'TWUN',    twun


!     5.2 TURBULENCE PARAMETERS
!     -------------------------

WRITE (20,9003) 'Turbulence model :'
WRITE (20,9005) 'PPA',     ppa
WRITE (20,9005) 'PPN',     ppn
WRITE (20,9005) 'PPVED0',  ppved0
WRITE (20,9005) 'PPVISBG', ppvisbg
WRITE (20,9005) 'PPDIFBG', ppdifbg
WRITE (20,9005) 'RMAXNUT', rmaxnut
WRITE (20,9005) 'AMVED0',  amved0
WRITE (20,9005) 'AMA',     ama
WRITE (20,9005) 'AMB',     amb
WRITE (20,9005) 'AMN1',    amn1
WRITE (20,9005) 'AMN2',    amn2
WRITE (20,9005) 'RMAXLAT', rmaxlat
WRITE (20,9005) 'ADK1',    adk1
WRITE (20,9005) 'ADK2',    adk2
WRITE (20,9005) 'ADCNU',   adcnu
WRITE (20,9005) 'ADLAM',   adlam
WRITE (20,9005) 'ADD1',    add1
WRITE (20,9005) 'ADD2',    add2
WRITE (20,9005) 'ADR1',    adr1
WRITE (20,9005) 'ADR2',    adr2
WRITE (20,9005) 'ABLAC',   ablac
WRITE (20,9005) 'BXING',   bxing
WRITE (20,9005) 'Z0BOT',   z0bot
WRITE (20,9005) 'Z0SUR',   z0sur
WRITE (20,9005) 'TKEMIN',  tkemin
CLOSE(20)

RETURN


!     8. FORMATS
!     ----------


9001 FORMAT ('///',1X,a,1X,'///')
9002 FORMAT (2X,a,t18,' : ',a)
9003 FORMAT (a)
9004 FORMAT (2X,a,t18,' : ',i8)
9005 FORMAT (2X,a,t18,' : ',1PE14.7)
9006 FORMAT (2X,a,t18,' : ',i8.8)
9007 FORMAT (2X,i3,2X,a1,2X,a11)
9008 FORMAT (2X,a,t18,' : ',l1)
9010 FORMAT (2X,a)
9011 FORMAT (3(3(2X,i6),/),3(2X,i8))
9012 FORMAT (3(2X,i8))


9101 FORMAT (2X,a6,'_',i1,'.',a3,i1,a1,' : ',i10,' bytes')
9102 FORMAT (2X,'Total amount of bytes : ',i10)

END SUBROUTINE wrcon


SUBROUTINE transh
!************************************************************************

!    *TRANSH*     EVALUATE GRID PARAMETERS AND RELATED ARRAYS

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 26 Mar 1999 @(COHERENS)transh.f 8.4

!       DESCRIPTION - EVALUATE HORIZONTAL GRID SPACINGS, GEOGRAPHICAL
!                     COORDINATES, CORIOLIS FREQUENCY AND OTHER WORK
!                     SPACE ARRAYS
!                   - GRID SPACINGS ARE IN m
!                   - CORIOLIS FREQUENCY IS CONSTANT IN CARTESIAN CASE
!                   - SPHCUR AND SPHCURV EQUAL TAN(PHI)/R (SPHERICAL) OR 0
!                     (CARTESIAN)
!                   - COSPHI AND COSPHIV EQUAL COS(PHI) (SPHERICAL) OR 1
!                     (CARTESIAN)

!       REFERENCE -

!       CALLING PROGRAM - INITC, PREPROC

!       EXTERNALS -

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, j
REAL :: coriun, raxis
!REAL :: gx0c, gy0c

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CORIUN*   REAL      UNIFORM CORIOLIS FREQUENCY [1/s]
!    *RAXIS*    REAL      DISTANCE OF GRID POINT TO EARTH'S ROTATION AXIS [m]

!------------------------------------------------------------------------

!     1. EVALUATE GRID SPACINGS
!     -------------------------

!     ---Cartesian grid

DO  i=1,nc
  DO  j=1,nr
    gx2(j,i) = gx0(i+1) - gx0(i)
  END DO
END DO
DO  j=1,nr
  gy2(j) = gy0(j+1) - gy0(j)
END DO


!     2. GEOGRAPHICAL COORDINATES, CORIOLIS AND SPHERICAL CURVATURE FACTORS
!     ---------------------------------------------------------------------

coriun = 4.0*pi*SIN(dlaref*conv)/86164.0
coriun = -1.0e-4

DO  i=1,nc
  dlon(i) = dloref
END DO
DO  j=1,nr
  dlat(j) = dlaref
  coriol(j)  = coriun
  coriolv(j) = coriun
  cosphi(j)  = 1.0
  cosphiv(j) = 1.0
  sphcur(j)  = 0.0
  sphcurv(j) = 0.0
END DO
cosphiv(nr+1) = 1.0
coriolv(nr+1) = coriun
sphcurv(nr+1) = 0.0


RETURN

END SUBROUTINE transh


! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:45:08

SUBROUTINE zero33(var,nk,nj,ni)
!************************************************************************

!    *ZERO*       INITIALISE A REAL ARRAY

!       AUTHOR - KEVIN RUDDICK

!       LAST UPDATE - 7 Dec 1994 @(COHERENS)Util.f 6.1

!       DESCRIPTION - THE ELEMENTS OF THE ARRAY VAR ARE SET TO ZERO

!       REFERENCE -

!       CALLING PROGRAM - VARIOUS

!       EXTERNALS - NONE

!************************************************************************

!*    ARGUMENTS

INTEGER, INTENT(IN)                      :: nk
INTEGER, INTENT(IN)                      :: nj
INTEGER, INTENT(IN)                      :: ni
REAL, INTENT(IN OUT)                        :: var(nk,nj,ni)



!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *VAR*      REAL      ARRAY TO BE INITIALISED                 (IN/OUT)
!    *NI*       INT       FIRST  DIMENSION OF VAR                     (IN)
!    *NJ*       INT       SECOND DIMENSION OF VAR                     (IN)
!    *NK*       INT       THIRD  DIMENSION OF VAR                     (IN)

!------------------------------------------------------------------------

!------------------------------------------------------------------------


!*    LOCAL VARIABLES

INTEGER :: i, j, k



!     1. SET ALL ARRAY ELEMENTS TO ZERO
!     ---------------------------------

DO  i=1,ni
  DO  j=1,nj
    DO  k=1,nk
      var(k,j,i) = 0.0
    END DO
  END DO
END DO


RETURN

END SUBROUTINE zero33

!========================================================================

SUBROUTINE izero33(ivar,nk,nj,ni)
!************************************************************************

!    *IZERO*      INITIALISE AN INTEGER ARRAY

!       AUTHOR - KEVIN RUDDICK

!       LAST UPDATE - 7 Dec 1994 @(COHERENS)Util.f 6.1

!       DESCRIPTION - THE ELEMENTS OF THE ARRAY IZERO ARE SET TO ZERO

!       REFERENCE -

!       CALLING PROGRAM - VARIOUS

!       EXTERNALS - NONE

!************************************************************************

!*    ARGUMENTS


INTEGER, INTENT(IN)                      :: nk
INTEGER, INTENT(IN)                      :: nj
INTEGER, INTENT(IN)                      :: ni
INTEGER, INTENT(IN OUT)                     :: ivar(nk,nj,ni)



!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *IVAR*     INT       ARRAY TO BE INITIALISED                 (IN/OUT)
!    *NI*       INT       FIRST  DIMENSION OF IVAR                    (IN)
!    *NJ*       INT       SECOND DIMENSION OF IVAR                    (IN)
!    *NK*       INT       THIRD  DIMENSION OF IVAR                    (IN)

!------------------------------------------------------------------------

!------------------------------------------------------------------------

!*    LOCAL VARIABLES

INTEGER :: i, j, k



!     1. SET ALL ARRAY ELEMENTS TO ZERO
!     ---------------------------------

DO  i=1,ni
  DO  j=1,nj
    DO  k=1,nk
      ivar(k,j,i) = 0.0
    END DO
  END DO
END DO


RETURN

END SUBROUTINE izero33

!=======================================================================

SUBROUTINE searho
!***********************************************************************

!    *SEARHO*      UPDATE DENSITY, BUOYANCY AND EXPANSION COEFFICIENTS

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 5 May 1999 @(COHERENS)densty.f 8.4

!       DESCRIPTION - UPDATES DENSITY AND BUOYANCY FOR LINEAR AND GENERAL
!                     EQUATION OF STATE
!                   - CALCULATES EXPANSION COEFFICIENTS FOR TEMPERATURE
!                     AND SALINITY AT ZERO PRESSURE USING THE GENERAL
!                     EQUATION OF STATE (IOPTD=2)

!       REFERENCE - Section III-1.5 of the User Documentation

!       CALLING PROGRAM - MAINPROG, PREPROC

!***********************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, j, k
REAL :: rdens

!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *A*        REAL      COEFFICIENTS IN GENERAL EQUATION OF STATE
!    *DRS0*     REAL      DERIVATIVE OF RHO WITH RESPECT TO SALINITY
!    *DRT0*     REAL      DERIVATIVE OF RHO WITH RESPECT TO TEMPERATURE

!-----------------------------------------------------------------------

!     1. LINEAR EQUATION OF STATE
!     ---------------------------
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
!          rdens = SBET(K,J,I)*(S(k,j,i)-SREF)-TBET(K,J,I)*(T(K,J,I)-TREF)        
          rdens = SBET(K,J,I)*S(k,j,i)-TBET(K,J,I)*T(K,J,I)      
          ro(k,j,i) = r0ref*(rdens+1.0)
          rdens = SBET(K,J,I)*s(k,j,i)-TBET(K,J,I)*t(k,j,i)   
          buoy(k,j,i) = -g*rdens
        END DO
      END IF
    END DO
  END DO

RETURN

END SUBROUTINE searho

SUBROUTINE rhoini
!***********************************************************************

!    *SEARHO*      UPDATE DENSITY, BUOYANCY AND EXPANSION COEFFICIENTS

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 5 May 1999 @(COHERENS)densty.f 8.4

!       DESCRIPTION - UPDATES DENSITY AND BUOYANCY FOR LINEAR AND GENERAL
!                     EQUATION OF STATE
!                   - CALCULATES EXPANSION COEFFICIENTS FOR TEMPERATURE
!                     AND SALINITY AT ZERO PRESSURE USING THE GENERAL
!                     EQUATION OF STATE (IOPTD=2)

!       REFERENCE - Section III-1.5 of the User Documentation

!       CALLING PROGRAM - MAINPROG, PREPROC

!***********************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, j, k
REAL :: rdens

!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *A*        REAL      COEFFICIENTS IN GENERAL EQUATION OF STATE
!    *DRS0*     REAL      DERIVATIVE OF RHO WITH RESPECT TO SALINITY
!    *DRT0*     REAL      DERIVATIVE OF RHO WITH RESPECT TO TEMPERATURE

!-----------------------------------------------------------------------

!     1. LINEAR EQUATION OF STATE
!     ---------------------------
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
!          rdens = SBET(K,J,I)*(S(k,j,i)-SREF)-TBET(K,J,I)*(T(K,J,I)-TREF)     
          rdens = SBET(K,J,I)*S(k,j,i)-TBET(K,J,I)*T(K,J,I)        
!          rdens = 0.
          roref(k,j,i) = r0ref*(rdens+1.0)
        END DO
      END IF
    END DO
  END DO

RETURN

END SUBROUTINE rhoini


! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:36:56

SUBROUTINE bstres
!************************************************************************

!    *BSTRES*      BOTTOM STRESS CALCULATIONS

!       AUTHOR - KEVIN RUDDICK AND PATRICK LUYTEN

!       LAST UPDATE - 6 May 1998 @(COHERENS)bstres.f 8.4

!       DESCRIPTION - CALCULATES BOTTOM STRESS, BOTTOM DRAG AND FRICTION
!                     COEFFICIENT FOR 2-D AND 3-D (SLIP) MOMENTUM EQUATIONS
!                   - WAVE EFFECTS ARE INCLUDED BY CALLING WAVCUR IF IOPTW>0

!       REFERENCE - Section III-1.6.2 of the User Documentation

!       CALLING PROGRAM - CRRNT3C, CRRNT3P, INCTUR, MAINPROG

!       EXTERNALS - WAVCUR

!************************************************************************
!*    LOCAL VARIABLES

LOGICAL :: call1
INTEGER :: i, j
!REAL :: cu2atv, cv2atu, u2atc, v2atc

SAVE call1
DATA call1 /.true./

!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL

!-----------------------------------------------------------------------

!     1. INITIALISE ROUGHNESS LENGTH ON FIRST CALL
!     --------------------------------------------

IF (call1) THEN
  IF (ibstr == 0) THEN
    WRITE (*,*) 'using zero bottom stress condition'
  ELSE IF (ibstr == 1) THEN
    DO  i=1,nc
      DO  j=1,nr
        IF (nwd(j,i) == 1) THEN
          cdb(j,i) = ((LOG(0.5*gz0(2,j,i)*dep(j,i)/cdz0(j,i))) /ckar)**(-2)
          cdb100(j,i) = (LOG(cdz0(j,i))/ckar)**(-2)
        ELSE
          cdb(j,i) = 0.0
          cdb100(j,i) = 0.0
        END IF
      END DO
    END DO
    WRITE (*,*) 'using quadratic bottom stress condition'
  ELSE IF (ibstr == 2) THEN
    WRITE (*,*) 'Using constant linear drag coefficient ', cdlin,' m/s'
   
  END IF
  !call1 = .false.
  
END IF




!     3. U-COMPONENT OF BOTTOM STRESS
!     -------------------------------

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      IF (ibstr <= 1) THEN
        fbk(j,i) = 0.5*(cdb(j,i)+cdb(j,i-1))*  &
            SQRT(u2(1,j,i)**2 + cv2atu(1,j,i)**2)
      ELSE IF (ibstr == 2) THEN
        fbk(j,i) = cdlin
      END IF
    END IF
    fb(j,i) = -fbk(j,i)*u2(1,j,i)
  END DO
END DO


!     4. V-COMPONENT OF BOTTOM STRESS
!     -------------------------------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      IF (ibstr <= 1) THEN
        gbk(j,i) = 0.5*(cdb(j,i)+cdb(j-1,i))*  &
            SQRT(v2(1,j,i)**2 + cu2atv(1,j,i)**2)
      ELSE IF (ibstr == 2) THEN
        gbk(j,i) = cdlin
      END IF
    END IF
    gb(j,i) = -gbk(j,i)*v2(1,j,i)
  END DO
END DO

!     5. TOTAL BOTTOM STRESS
!     ----------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      IF (ibstr <= 1) THEN
        bstot(j,i) = cdb(j,i)*(u2atc(1,j,i)**2+v2atc(1,j,i)**2)
      ELSE IF (ibstr == 2) THEN
        bstot(j,i) = cdlin*SQRT(u2atc(1,j,i)**2+v2atc(1,j,i)**2)
      END IF
    ELSE
      bstot(j,i) = 0.0
    END IF
  END DO
END DO

RETURN

END SUBROUTINE bstres

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:43:44

SUBROUTINE densty
!************************************************************************

!    *DENSTY*      CALCULATE BAROCLINIC PRESSURE GRADIENT (3-D APPLICATION)

!       AUTHOR - KEVIN RUDDICK AND PATRICK LUYTEN

!       LAST UPDATE - 5 May 1999 @(COHERENS)densty.f 8.4

!       DESCRIPTION - CALCULATES BAROCLINIC PRESSURE GRADIENT (UQDEN/VQDEN)
!                     AND ITS DEPTH-INTEGRAL (UDQDEN/VDQDEN)

!       REFERENCE - Section III-4.3.9 of the User Documentation

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS -

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, j, k
REAL :: uqden1(nz), uqden2(nz+1), vqden1(nz), vqden2(nz+1)
!REAL :: gx2u, gy2v, gz2u, gz2v, h2atc

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *RPRESS*   REAL      REDUCED PRESSURE (W NODES) (m2/s2)
!    *UQDEN1*   REAL      (MINUS) THE FIRST TERM ON THE R.H.S OF (3.4.140)
!    *UQDEN2*   REAL      THE FACTOR G1 AS GIVEN BY (3.4.143)
!    *VQDEN1*   REAL      (MINUS) THE FIRST TERM ON THE R.H.S OF (3.4.141)
!    *VQDEN2*   REAL      THE FACTOR G2 AS GIVEN BY (3.4.144)

!------------------------------------------------------------------------

!     1. REDUCED PRESSURE AT W NODES FROM HYDROSTATIC EQUATION
!     --------------------------------------------------------

DO  i=1,nc
  DO  j=1,nr
    rpress(nz+1,j,i) = 0.0
    IF (nwd(j,i) == 1) THEN
      DO  k=nz,1,-1
        rpress(k,j,i) = rpress(k+1,j,i) - buoy(k,j,i)*gz2(k,j,i)
      END DO
    ELSE
      DO  k=1,nz
        rpress(k,j,i) = 0.0
      END DO
    END IF
  END DO
END DO

!     2. CALCULATE HORIZONTAL GRADIENTS OF REDUCED PRESSURE
!     -----------------------------------------------------

!     2.1 INTERNAL U NODES
!     --------------------

DO  i=2,nc
  DO  j=1,nr
    
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        uqden1(k) = 0.5*((rpress(k,j,i)   + rpress(k+1,j,i))*gz2(k,j,i)-  &
            (rpress(k,j,i-1) + rpress(k+1,j,i-1))*gz2(k,j,i-1))/ gx2u(j,i)
      END DO
      
      uqden2(nz+1) = 0.0
      DO  k=1,nz
        uqden2(k) = -0.5*(rpress(k,j,i) + rpress(k,j,i-1))*  &
            (gz0(k,j,i)*(h2atc(j,i)-h2atc(j,i-1)) - (dep(j,i)-dep(j,i-1)))/gx2u(j,i)
      END DO
      
      DO  k=1,nz
        uqden(k,j,i) =-( uqden1(k) + (uqden2(k+1)-uqden2(k))) / gz2u(k,j,i) - uqden0(k,j,i)
      END DO
      
    END IF
  END DO
END DO

!     2.2 INTERNAL V NODES
!     --------------------

DO  i=1,nc
  DO  j=2,nr
    
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        vqden1(k) =  &
            0.5*((rpress(k,j  ,i) + rpress(k+1,j  ,i))*gz2(k,j  ,i) -  &
            (rpress(k,j-1,i) + rpress(k+1,j-1,i))*gz2(k,j-1,i))/ gy2v(j)
      END DO
      
      vqden2(nz+1) = 0.0
      DO  k=1,nz
        vqden2(k) = -0.5*(rpress(k,j,i) + rpress(k,j-1,i))*  &
            (gz0(k,j,i)*(h2atc(j,i)-h2atc(j-1,i)) -(dep(j,i)-dep(j-1,i)))/gy2v(j)
      END DO
      
      DO  k=1,nz
        vqden(k,j,i) =-( vqden1(k) + (vqden2(k+1)-vqden2(k))) / gz2v(k,j,i) - vqden0(k,j,i)
      END DO

    END IF
  END DO

END DO

!     3. DEPTH-INTEGRATE HORIZONTAL GRADIENTS OF REDUCED PRESSURE
!     -----------------------------------------------------------

!     3.1 INTERNAL U NODES
!     --------------------

DO  i=2,nc
  DO  j=1,nr
    udqden(j,i) = 0.0
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        udqden(j,i) = udqden(j,i) + gz2u(k,j,i)*uqden(k,j,i)
      END DO
    END IF
  END DO
END DO


!     3.2 INTERNAL V NODES
!     --------------------

DO  i=1,nc

  DO  j=2,nr
    vdqden(j,i) = 0.0
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        vdqden(j,i) = vdqden(j,i) + gz2v(k,j,i)*vqden(k,j,i)
      END DO
    END IF
  END DO
END DO


RETURN

END SUBROUTINE densty



SUBROUTINE transv
!************************************************************************

!    *TRANSV*     UPDATE VERTICAL GRID SPACINGS AT NEW 3-D TIME STEP

!       AUTHOR - KEVIN RUDDICK

!       LAST UPDATE - 27 Sep 1995 @(COHERENS)transv.f 7.0

!       DESCRIPTION - THE VERTICAL GRID SPACING EQUALS THE JACOBIAN J

!       REFERENCE -

!       CALLING PROGRAM - INITC, MAINPROG, PREPROC

!       EXTERNALS -

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, j, k
!REAL :: h2atc



!     1. PERFORM SIGMA TRANSFROMATION (AT EACH TIME STEP)
!     ---------------------------------------------------

DO  i=1,nc
  DO  j=1,nr
    DO  k=1,nz
      gz2(k,j,i) = (gz0(k+1,j,i)-gz0(k,j,i))*h2atc(j,i)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE transv

!
SUBROUTINE metin
!************************************************************************

!    *METIN*       READ MET DATA, UPDATE SURFACE FLUXES AND SOLAR IRRADIANCE

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 25 May 1999 @(COHERENS)metin.f 8.4

!       DESCRIPTION - READS IN ATMOSPHERIC INPUT DATA FROM 'met'-FILE AT
!                     SELECTED INTERVALS
!                   - TYPE AND FORM OF INPUT DATA SELECTED BY IOPTM
!                   - UPDATES SURFACE FLUXES
!                   - EVALUATES SOLAR IRRADIANCE

!       REFERENCE -

!       CALLING PROGRAM - INCTUR, MAINPROG

!       EXTERNALS - ERROR, RDARR, SOLRAD, SURFLX

!************************************************************************
!*    LOCAL VARIABLES

LOGICAL :: call1, erflag
CHARACTER (LEN=60) :: metfil
INTEGER :: i, iioptm, j
REAL :: cloud2un, hum2un, prevapun, qnsolun, qsolun
REAL :: sat2un, windu2un, windv2un,fff
SAVE call1
DATA call1 /.true./

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL
!    *ERFLAG*   LOGICAL   ERROR FLAG (TRUE IF UNITS IOMET OR IOWAV ARE
!                         CONNECTED TO A FILE)
!    *METFIL*   CHAR      NAME OF MET FIL
!    *CLOUD2UN* REAL      HORIZ. UNIFORM CLOUD COVERAGE
!    *HUM2UN*   REAL      HORIZ. UNIFORM RELATIVE HUMIDITY
!    *PREVAPUN* REAL
!                IOPTSA = 1 -> HORIZ. UNIFORM PRECIP - EVAP RATE [kg/(m2*s)]
!                IOPTSA = 2 -> HORIZ. UNIFORM PRECIPITATION RATE [kg/(m2*s)]
!    *QNSOLUN*  REAL      HORIZ. UNIFORM NON-SOLAR HEAT FLUX [W/m2]
!    *QSOLUN*   REAL      HORIZ. UNIFORM SOLAR HEAT FLUX [W/m2]
!    *SAT2UN*   REAL      HORIZ. UNIFORM SEA SURFACE TEMPERATURE [deg C]
!    *WINDU2UN* REAL      X-COMPONENT OF HORIZ. UNIFORM 10m WIND [m/s]
!    *WINDV2UN* REAL      Y-COMPONENT OF HORIZ. UNIFORM 10m WIND [m/s]

!------------------------------------------------------------------------

!     1. READ FILE HEADER ON FIRST CALL,INITIALISE FLUXES AND WAVE DATA
!     -----------------------------------------------------------------

!     2. READ MET DATA FOR NEW TIME STEP
!     ----------------------------------


!        --read data
!         IF (METFORM.EQ.'A') THEN
!            READ (IOMET,9002) WINDU2UN, WINDV2UN, SAT2UN,
!     1                        HUM2UN,   CLOUD2UN, PREVAPUN

fff = REAL(nt)/(24.*3600.)
fff = AMIN1(fff,1.)

sat2un = 20.
hum2un = 0.5
cloud2un = 0.
prevapun = 0.

DO  i=1,nc
  DO  j=1,nr
    windu2(j,i) = 0.
    windv2(j,i) = 0.
    sat2(j,i) = sat2un
    sst2(j,i) = t(nz,j,i)+tref
    hum2(j,i) = hum2un
    cloud2(j,i) = cloud2un
    rain2(j,i) = prevapun
  END DO
END DO
!     3. ZERO VALUES ON DRY CELLS
!     ---------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 0) THEN
      windu2(j,i) = 0.0
      windv2(j,i) = 0.0
      p2(j,i) = 0.0
      sat2(j,i) = 0.0
      hum2(j,i) = 0.0
      cloud2(j,i) = 0.0
      evapr(j,i) = 0.0
      rain2(j,i) = 0.0
    END IF
  END DO
END DO


!     4. UPDATE SURFACE STRESS AND HEAT FLUXES
!     ----------------------------------------


!     ---non-solar heat and salinity fluxes
CALL surflx
!     ---solar heat fluxes
CALL solrad

RETURN
END SUBROUTINE metin

!=======================================================================

SUBROUTINE surflx
!***********************************************************************

!    *SURFLX*      EVALUATE SURFACE STRESS, NON-SOLAR HEAT FLUX AND
!                  SALINITY FLUX AT SEA SURFACE

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE -  25 May 1999 @(COHERENS)metin.f 8.4

!       DESCRIPTION - EVALUATE SURFACE STESS
!                   - EVALUATE LATENT, SENSIBLE AND LONG-WAVE RADIATION FLUXES
!                   - EVALUATE SALINITY FLUXES

!       REFERENCE - Section III-1.6.1 of the User Documentation

!       CALLING PROGRAM - METIN

!       EXTERNALS - CD, CE, FLUXCO, HUMID

!***********************************************************************
!*    LOCAL VARIABLES

REAL, PARAMETER :: rhoa = 1.28
LOGICAL :: call1
REAL :: drag(13,11), evap(nr,nc), exch(13,11)
INTEGER :: i, j
REAL :: avap, cds, ces, cpair, hvap, qa, qlat, qnlw, qs, qsen
REAL :: svap, tdif, tkelv, tv, wind
!REAL :: cd, ce

SAVE call1, drag, exch
DATA call1 /.true./

!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL
!    *AVAP*     REAL      SATURATED VAPOUR PRESSURE AT AIR TEMPERATURE [mb]
!    *CDS*      REAL      SURFACE DRAG COEFFICIENT
!    *CES*      REAL      HEAT EXCHANGE COEFFICIENT
!    *CPAIR*    REAL      SPECIFIC HEAT OF AIR                   [J/(kg K)]
!    *DRAG*     REAL      ARRAY WHICH CONTAINS VALUES OF DRAG COEFFICIENT
!                         AT SPECIFIED INTERVALS OF WIND SPEED AND AIR/SEA
!                         TEMPERATURE DIFFERENCE (ITDIF=1 ONLY)
!    *EVAP*     REAL      EVAPORATION RATE                        [kg/m2/s]
!    *EXCH*     REAL      ARRAY WHICH CONTAINS VALUES OF THERMAL EXCHANGE
!                         COEFFICIENT AT SPECIFIED INTERVALS OF WIND SPEED
!                         AND AIR/SEA TEMPERATURE DIFFERENCE (ITDIF=1 ONLY)
!    *HVAP*     REAL      LATENT HEAT OF VAPORIZATION                [J/kg]
!    *QA*       REAL      SPECIFIC HUMIDITY OF AIR
!    *QLAT*     REAL      LATENT HEAT FLUX                           [W/m2]
!    *QNLW*     REAL      LONG WAVE RADIATION FLUX                   [W/m2]
!    *QS*       REAL      SPECIFIC HUMIDITY AT SEA SURFACE TEMPERATURE
!    *QSEN*     REAL      SENSIBLE HEAT FLUX                         [W/m2]
!    *RHOA*     REAL      AIR DENSITY                               [kg/m3]
!    *SVAP*     REAL      VAPOUR PRESSURE AT SEA SURFACE TEMPERATURE [mbar]
!    *TDIF*     REAL      AIR-SEA TEMPERATURE DIFFERENCE            [deg C]
!    *TKELV*    REAL      SEA SURFACE TEMPERATURE                   [deg K]
!    *TV*       REAL      TEMPORARY WORK SPACE
!    *WIND*     REAL      WIND SPEED AT 10 m HEIGHT                   [m/s]

!-----------------------------------------------------------------------

!     1. INITIALISE ARRAYS ON FIRST CALL
!     ----------------------------------


IF (call1) THEN
  IF (iopthe == 0.AND.itdif == 1) THEN
    itdif = 0
    WRITE (0,'(A)') 'WARNING : switch ITDIF is set to '//  &
        'zero since heat equation is not solved'
  END IF
  IF (itdif == 1) THEN
 !   CALL fluxco(drag,exch)
  END IF
  call1 = .false.
END IF


!     2. SURFACE STRESS
!     -----------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      wind = SQRT(windu2(j,i)**2+windv2(j,i)**2)
      tdif = sat2(j,i)-sst2(j,i)
!      cds = cd(wind,tdif,drag,itdif)
      cds = 1.25e-3    
      fs(j,i) = rhoa*cds*wind*windu2(j,i)/r0ref
      gs(j,i) = rhoa*cds*wind*windv2(j,i)/r0ref
      sstot(j,i) = rhoa*cds*wind*wind/r0ref
    END IF
  END DO
END DO

!     3. NON-SOLAR HEAT FLUXES
!     ------------------------

DO  i=1,nc
  DO  j=1,nr
      qnsol(j,i) = 0.
  END DO
END DO

!     4. SALINITY FLUX
!     ----------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
!               SSALFL(J,I) = (S(NZ,J,I)+SREF)*(EVAP(J,I)-RAIN2(J,I))
!     1                      /((1.0-0.001*S(NZ,J,I))*R0REF)
      
      ssalfl(j,i) = 0.
      
    END IF
  END DO
END DO

RETURN

END SUBROUTINE surflx

!=======================================================================

SUBROUTINE solrad

!*    LOCAL VARIABLES

INTEGER :: i, j, k

DO  i=1,nc
  DO  j=1,nr
    
       qsol(j,i) = 0.0
   
  END DO
END DO


RETURN

END SUBROUTINE solrad

SUBROUTINE dislen
!************************************************************************

!    *DISLEN*   EVALUATE DISSIPATION RATE (K-L THEORY)

!     AUTHOR - PATRICK LUYTEN

!     LAST UPDATE  - 7 May 1998         @(COHERENS)turben.f 8.4

!     DESCRIPTION - EVALUATES DISSIPATION RATE [W/kg] AS A FUNCTION OF
!                   MIXING LENGTH

!     REFERENCE - Section III-1.2.2 of the User Documentation

!     CALLING PROGRAM - INCTUR, VEDDY2

!     EXTERNALS -

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, j, k


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        IF (zlw(k,j,i) > 0.0) THEN
          dissw(k,j,i) = eps0*tkew(k,j,i)**1.5/zlw(k,j,i)
        ELSE
          dissw(k,j,i) = 0.0
        END IF
      END DO
    END IF
  END DO
END DO


RETURN

END SUBROUTINE dislen

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 11:32:41

SUBROUTINE tleng
!************************************************************************

!    *TLENG*      ALGEBRAIC MIXING LENGTH FORMULATIONS

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 7 May 1998        @(COHERENS)tleng.f 8.4

!       DESCRIPTION - EVALUATES TURBULENCE MIXING LENGTH USING AN ALGEBRAIC
!                     PRESCRIPTION WITHOUT SOLVING A SEPARATE TRANSPORT EQN.
!                   - LENGTH SCALE PRESCRIPTION SELECTED BY ILENG
!                   1) ILENG = 1 : PARABOLIC LAW
!                   2) ILENG = 2 : "MODIFIED" PARABOLIC LAW
!                   3) ILENG = 3 : "XING" FORMULATION
!                   4) ILENG = 4 : "BLACKADAR" FORMULATION

!       REFERENCE - Section III-1.2.2b of the User Documentation
!        1) Xing J. and Davies A.M., 1996. Application of a range of
!           turbulence energy models to the determination of M4 tidal current
!           profiles. Cont. Shelf Res., 16, 517-547.
!        2) Blackadar A.K., 1962. The vertical distribution of wind and
!           turbulent exchange in a neutral atmosphere. J. Geophys. Res., 67,
!           3095-3102.
!        3) Mellor G.L. and Yamada T., 1974. A hierarchy of turbulence closure
!           models for planetary boundary layers. J. Atmos. Sc., 31, 1791-1806.

!       CALLING PROGRAM - INCTUR, VEDDY2

!       EXTERNALS -

!************************************************************************
!*    LOCAL VARIABLES

LOGICAL :: call1
INTEGER :: i, j, k
REAL :: botint, topint, zexp, zlasym, zl1, zl2, z1
!REAL :: gzsc, gz2w, h2atc

SAVE call1
DATA  call1 /.true./

!-----------------------------------------------------------------------
!     NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL
!    *ZLASYM*   REAL      BLACKADAR'S ASYMPTOTIC MIXING LENGTH [m]
!    *ZEXP*     REAL      XING'S EXPONENTAL FACTOR

!-----------------------------------------------------------------------

!     1. OUTPUT MESSAGES ON FIRST CALL
!     --------------------------------


IF (call1) THEN
  IF (ileng == 1) THEN
    WRITE (*,*) 'Parabolic mixing length'
  ELSE IF (ileng == 2) THEN
    WRITE (*,*) '"Modified" parabolic mixing length'
  ELSE IF (ileng == 3) THEN
    WRITE (*,*) '"Xing" mixing length'
  ELSE IF (ileng == 4) THEN
    WRITE (*,*) '"Blackadar" mixing length'
  END IF
  call1 = .false.
END IF


!     2. EVALUATE TURBULENCE LENGTH SCALE
!     -----------------------------------

!     2.1 PARABOLIC LAW
!     ---------------------------

IF (ileng == 1) THEN
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz+1
          zl1 = gz0(k,j,i)*h2atc(j,i) + z0bot
          zl2 = h2atc(j,i)*(1.0-gz0(k,j,i)) + z0sur
          zlw(k,j,i) = (ckar*zl1*zl2)/(zl1+zl2)
        END DO
      END IF
    END DO
  END DO
  211   CONTINUE
  
  
!     2.2 "MODIFIED" PARABOLIC LAW
!     ---------------------------
  
ELSE IF (ileng == 2) THEN
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz+1
          zl1 = gz0(k,j,i)*h2atc(j,i) + z0bot
          zl2 = h2atc(j,i)*(1.0-gz0(k,j,i)) + z0sur
          zlw(k,j,i) = (ckar*zl1)*SQRT(zl2/(zl1+zl2))
        END DO
      END IF
    END DO
  END DO
  221   CONTINUE
  
  
!     2.3 "XING" MIXING LENGTH
!     ------------------------
  
  
ELSE IF (ileng == 3) THEN
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        zexp = EXP(bxing*gz0(1,j,i))
        zl1 = zexp*gz0(1,j,i)*h2atc(j,i) + z0bot
        zl2 = h2atc(j,i)*(1.0-gz0(1,j,i)) + z0sur
        zlw(1,j,i) = (ckar*zl1*zl2)/(zl1+zl2)
        DO  k=2,nz+1
          zexp  = zexp*EXP(bxing*gzsc(k-1,j,i))
          zl1 = zexp*gz0(k,j,i)*h2atc(j,i) + z0bot
          zl2 = h2atc(j,i)*(1.0-gz0(k,j,i)) + z0sur
          zlw(k,j,i) = (ckar*zl1*zl2)/(zl1+zl2)
        END DO
      END IF
    END DO
  END DO
  231   CONTINUE
  
  
!     2.4 "BLACKADAR" FORMULA
!     -----------------------
  
  
ELSE IF (ileng == 4) THEN
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        
!        ---asymptotic length scale
        z1 = SQRT(tkew(1,j,i))*0.5*gz2(1,j,i)
        topint = (1.0-0.25*gzsc(1,j,i))*z1
        botint = z1
        DO  k=2,nz
          z1 = SQRT(tkew(k,j,i))*gz2w(k,j,i)
          topint = topint +(1.0-gz0(k,j,i))*z1
          botint = botint + z1
        END DO
        z1 = SQRT(tkew(nz+1,j,i))*0.5*gz2(nz,j,i)
        topint = topint + 0.25*gzsc(nz,j,i)*z1
        botint = botint + z1
        IF (botint > 0.0) THEN
          zlasym = ablac*h2atc(j,i)*topint/botint
        ELSE
          zlasym = 0.5*ablac*h2atc(j,i)
        END IF
        
!        ---length scale
        DO  k=1,nz+1
          zl1 = gz0(k,j,i)*h2atc(j,i) + z0bot
          zl2 = h2atc(j,i)*(1.0-gz0(k,j,i)) + z0sur
          zlw(k,j,i) = (ckar*zl1*zl2)/(zl1+zl2+(ckar*zl1*zl2)/zlasym)
        END DO
        
      END IF
    END DO
  END DO
  241   CONTINUE
  
END IF


RETURN

END SUBROUTINE tleng


END MODULE cohini