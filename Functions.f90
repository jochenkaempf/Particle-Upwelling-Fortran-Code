MODULE functions

USE param

CONTAINS

!************************************************************************

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:40:50

!    *Intp_in.f*  INTERPOLATION FUNCTIONS

!       AUTHOR - KEVIN RUDDICK

!       LAST UPDATE - 6 Nov 1997 @(COHERENS)Interp.f 8.1

!       DESCRIPTION - PERFORMS SIMPLE INTERPOLATIONS OF GRID AND
!                     PHYSICAL VARIABLES TO OBTAIN THOSE VALUES NOT
!                     STORED DIRECTLY ON AN ARAKAWA 'C' GRID.
!                   - FOR EFFICIENT RUNNING THESE FUNCTIONS SHOULD
!                     BE COMPILED INLINE. IF THIS IS NOT POSSIBLE
!                     OWING TO COMPILER LIMITATIONS, IT MAY BE PREFERABLE
!                     TO DEFINE THESE QUANTITIES AS COMMON ARRAYS
!                     RATHER THAN AS THE FOLLOWING FUNCTIONS.
!                   - THE TYPE OF GRID LOCATION WHERE THE QUANTITY IS
!                     INTERPOLATED, IS INDICATED BY THE LAST LETTER(S)
!                     OF THE FUNCTION NAME:
!                      C:  CELL CENTRE
!                      U:  U-NODE
!                      V:  V-NODE
!                      W:  W-NODE
!                      UX: X-NODE
!                      VX: V-NODE

!       REFERENCE - Section V-1.12 of the User Documentation

!       CALLING PROGRAM - VARIOUS

!       EXTERNALS -

!************************************************************************

!=======================================================================

REAL FUNCTION cud2atv(j,i)

!     'wet points only' one-sided Coriolis interpolation adjacent to
!     solid boundaries
!*    ARGUMENTS

INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                  :: i

!*    LOCAL VARIABLES

INTEGER :: ncor
REAL :: sum


sum = 0.0
ncor= 0
IF (j <= nr) THEN
  IF (npix(j,i) == 1) THEN
    sum = sum + ud2(j,i)
    ncor= ncor+ 1
  END IF
  IF (npix(j,i+1) == 1) THEN
    sum = sum + ud2(j,i+1)
    ncor= ncor+ 1
  END IF
END IF
IF (j >= 2) THEN
  IF (npix(j-1,i) == 1) THEN
    sum = sum + ud2(j-1,i)
    ncor= ncor+ 1
  END IF
  IF (npix(j-1,i+1) == 1) THEN
    sum = sum + ud2(j-1,i+1)
    ncor= ncor+ 1
  END IF
END IF

IF (ncor > 0) THEN
  cud2atv = sum/ncor
ELSE
  cud2atv = 0.0
END IF

RETURN

END FUNCTION cud2atv

!=======================================================================

REAL FUNCTION cu2atv(k,j,i)

!     'wet points only' one-sided Coriolis interpolation adjacent to
!     solid boundaries
!*    ARGUMENTS


INTEGER, INTENT(IN)                      :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                      :: i


!*    LOCAL VARIABLES

INTEGER :: ncor
REAL :: sum

sum = 0.0
ncor= 0
IF (j <= nr) THEN
  IF (npix(j,i) == 1) THEN
    sum = sum + u2(k,j,i)
    ncor= ncor+ 1
  END IF
  IF (npix(j,i+1) == 1) THEN
    sum = sum + u2(k,j,i+1)
    ncor= ncor+ 1
  END IF
END IF
IF (j >= 2) THEN
  IF (npix(j-1,i) == 1) THEN
    sum = sum + u2(k,j-1,i)
    ncor= ncor+ 1
  END IF
  IF (npix(j-1,i+1) == 1) THEN
    sum = sum + u2(k,j-1,i+1)
    ncor= ncor+ 1
  END IF
END IF

IF (ncor > 0) THEN
  cu2atv = sum/ncor
ELSE
  cu2atv = 0.0
END IF

RETURN

END FUNCTION cu2atv

!=======================================================================

REAL FUNCTION cvd2atu(j,i)

!     'wet points only' one-sided Coriolis interpolation adjacent to
!     solid boundaries

!*    ARGUMENTS


INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                      :: i


!*    LOCAL VARIABLES

INTEGER :: ncor
REAL :: sum


sum = 0.0
ncor= 0
IF (i <= nc) THEN
  IF (npiy(j,i) == 1) THEN
    sum = sum + vd2(j,i)
    ncor= ncor+ 1
  END IF
  IF (npiy(j+1,i) == 1) THEN
    sum = sum + vd2(j+1,i)
    ncor= ncor+ 1
  END IF
END IF
IF (i >= 2) THEN
  IF (npiy(j,i-1) == 1) THEN
    sum = sum + vd2(j,i-1)
    ncor= ncor+ 1
  END IF
  IF (npiy(j+1,i-1) == 1) THEN
    sum = sum + vd2(j+1,i-1)
    ncor= ncor+ 1
  END IF
END IF

!IF((j==2).and.(i==2))THEN
! write(6,*)"I'M INSIDE", vd2(j,i),vd2(j+1,i),vd2(j,i-1),vd2(j+1,i-1)
!end IF

IF (ncor > 0) THEN
  cvd2atu = sum/ncor
ELSE
  cvd2atu = 0.0
END IF

RETURN

END FUNCTION cvd2atu

!=======================================================================

REAL FUNCTION cv2atu(k,j,i)

!     'wet points only' one-sided Coriolis interpolation adjacent to
!     solid boundaries

!*    ARGUMENTS


INTEGER, INTENT(IN)                         :: k
INTEGER, INTENT(IN)                         :: j
INTEGER, INTENT(IN)                         :: i


!*    LOCAL VARIABLES

INTEGER :: ncor
REAL :: sum


sum = 0.0
ncor= 0
IF (i <= nc) THEN
  IF (npiy(j,i) == 1) THEN
    sum = sum + v2(k,j,i)
    ncor= ncor+ 1
  END IF
  IF (npiy(j+1,i) == 1) THEN
    sum = sum + v2(k,j+1,i)
    ncor= ncor+ 1
  END IF
END IF
IF (i >= 2) THEN
  IF (npiy(j,i-1) == 1) THEN
    sum = sum + v2(k,j,i-1)
    ncor= ncor+ 1
  END IF
  IF (npiy(j+1,i-1) == 1) THEN
    sum = sum + v2(k,j+1,i-1)
    ncor= ncor+ 1
  END IF
END IF

IF (ncor > 0) THEN
  cv2atu = sum/ncor
ELSE
  cv2atu = 0.0
END IF

RETURN

END FUNCTION cv2atu

!=======================================================================

REAL FUNCTION depatu(j,i)

INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                      :: i


depatu = 0.5*(dep(j,i)+dep(j,i-1))

RETURN

END FUNCTION depatu

!=======================================================================

REAL FUNCTION depatv(j,i)


INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                  :: i


depatv = 0.5*(dep(j,i)+dep(j-1,i))

RETURN

END FUNCTION depatv

!=======================================================================

REAL FUNCTION gx0c(i)

INTEGER, INTENT(IN)                  :: i


gx0c = 0.5*(gx0(i)+gx0(i+1))

RETURN

END FUNCTION gx0c

!=======================================================================

REAL FUNCTION gx2u(j,i)

INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                      :: i


gx2u = 0.5*( gx2(j,i) + gx2(j,i-1) )

RETURN

END FUNCTION gx2u

!=======================================================================

REAL FUNCTION gx2v(j,i)


INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                  :: i


gx2v = 0.5*( gx2(j,i) + gx2(j-1,i) )

RETURN

END FUNCTION gx2v

!=======================================================================

REAL FUNCTION gy0c(j)

INTEGER, INTENT(IN)                  :: j


gy0c = 0.5*(gy0(j)+gy0(j+1))

RETURN

END FUNCTION gy0c

!=======================================================================

REAL FUNCTION gy2u(j)

INTEGER, INTENT(IN)                  :: j


gy2u = gy2(j)

RETURN

END FUNCTION gy2u

!=======================================================================

REAL FUNCTION gy2v(j)

INTEGER, INTENT(IN)                  :: j


gy2v = 0.5*( gy2(j) + gy2(j-1) )

RETURN

END FUNCTION gy2v

!=======================================================================

REAL FUNCTION gzsc(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                  :: i


gzsc = gz0(k+1,j,i)-gz0(k,j,i)

RETURN

END FUNCTION gzsc

!=======================================================================

REAL FUNCTION gz0c(k,j,i)


INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                  :: i

gz0c = 0.5*(gz0(k,j,i) + gz0(k+1,j,i))

RETURN

END FUNCTION gz0c

!=======================================================================

REAL FUNCTION gz1u(k,j,i)


INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                      :: i


gz1u = 0.5*( gz1(k,j,i) + gz1(k,j,i-1) )

RETURN

END FUNCTION gz1u

!=======================================================================

REAL FUNCTION gz1v(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i


gz1v = 0.5*( gz1(k,j,i) + gz1(k,j-1,i) )

RETURN

END FUNCTION gz1v

!=======================================================================

REAL FUNCTION gz1w(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i


gz1w = 0.5*( gz1(k,j,i) + gz1(k-1,j,i) )

RETURN

END FUNCTION gz1w

!=======================================================================

REAL FUNCTION gz2u(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                      :: i


gz2u = 0.5*( gz2(k,j,i) + gz2(k,j,i-1) )

RETURN

END FUNCTION gz2u

!=======================================================================

REAL FUNCTION gz2ux(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                      :: i


gz2ux = 0.25*( gz2(k,j,i) + gz2(k,j,i-1) + gz2(k-1,j,i) + gz2(k-1,j,i-1) )

RETURN

END FUNCTION gz2ux

!=======================================================================

REAL FUNCTION gz2v(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i


gz2v = 0.5*( gz2(k,j,i) + gz2(k,j-1,i) )

RETURN

END FUNCTION gz2v

!=======================================================================

REAL FUNCTION gz2vx(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i


gz2vx = 0.25*( gz2(k,j,i) + gz2(k,j-1,i) + gz2(k-1,j,i) + gz2(k-1,j-1,i) )

RETURN

END FUNCTION gz2vx

!=======================================================================

REAL FUNCTION gz2w(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i


gz2w = 0.5*( gz2(k,j,i) + gz2(k-1,j,i) )

RETURN

END FUNCTION gz2w

!=======================================================================

REAL FUNCTION h1atc(j,i)

INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                  :: i


h1atc = dep(j,i) + zeta1(j,i)

RETURN

END FUNCTION h1atc

!=======================================================================

REAL FUNCTION h1atu(j,i)

INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                  :: i


h1atu = 0.5*( dep(j,i) + zeta1(j,i) + dep(j,i-1) + zeta1(j,i-1) )

RETURN

END FUNCTION h1atu

!=======================================================================

REAL FUNCTION h1atv(j,i)

INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                  :: i


h1atv = 0.5*( dep(j,i) + zeta1(j,i) + dep(j-1,i) + zeta1(j-1,i) )

RETURN

END FUNCTION h1atv

!=======================================================================

REAL FUNCTION h2atc(j,i)

INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                  :: i


h2atc = dep(j,i) + zeta2(j,i)

RETURN

END FUNCTION h2atc

!=======================================================================

REAL FUNCTION h2atu(j,i)

INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                      :: i


h2atu = 0.5*( dep(j,i) + zeta2(j,i) + dep(j,i-1) + zeta2(j,i-1) )

RETURN

END FUNCTION h2atu

!=======================================================================

REAL FUNCTION h2atv(j,i)

INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                  :: i


h2atv = 0.5*( dep(j,i) + zeta2(j,i) + dep(j-1,i) + zeta2(j-1,i) )

RETURN

END FUNCTION h2atv

!=======================================================================

REAL FUNCTION ud2atc(j,i)

INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                      :: i


ud2atc = 0.5*( ud2(j,i) + ud2(j,i+1) )

RETURN

END FUNCTION ud2atc

!=======================================================================

REAL FUNCTION ud2atv(j,i)

INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                      :: i


ud2atv = 0.25*( ud2(j,i) + ud2(j,i+1) + ud2(j-1,i) + ud2(j-1,i+1) )

RETURN

END FUNCTION ud2atv

!=======================================================================

REAL FUNCTION u1atc(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                      :: i


u1atc = 0.5*( u1(k,j,i) + u1(k,j,i+1) )

RETURN

END FUNCTION u1atc

!=======================================================================

REAL FUNCTION u1atv(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                      :: i


u1atv = 0.25*( u1(k,j,i) + u1(k,j,i+1) + u1(k,j-1,i) + u1(k,j-1,i+1) )

RETURN

END FUNCTION u1atv

!=======================================================================

REAL FUNCTION u2atc(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                      :: i


u2atc = 0.5*( u2(k,j,i) + u2(k,j,i+1) )

RETURN

END FUNCTION u2atc

!=======================================================================

REAL FUNCTION u2atv(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                      :: i


u2atv = 0.25*( u2(k,j,i) + u2(k,j,i+1) + u2(k,j-1,i) + u2(k,j-1,i+1) )

RETURN

END FUNCTION u2atv

!=======================================================================

REAL FUNCTION vd2atc(j,i)

INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                  :: i


vd2atc = 0.5*( vd2(j,i) + vd2(j+1,i) )

RETURN

END FUNCTION vd2atc

!=======================================================================

REAL FUNCTION vd2atu(j,i)

INTEGER, INTENT(IN)                  :: j
INTEGER, INTENT(IN)                  :: i


vd2atu = 0.25*( vd2(j,i) + vd2(j+1,i) + vd2(j,i-1) + vd2(j+1,i-1) )

RETURN

END FUNCTION vd2atu

!========================================================================

REAL FUNCTION v1atc(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i


v1atc = 0.5*( v1(k,j,i) + v1(k,j+1,i) )

RETURN

END FUNCTION v1atc

!=======================================================================

REAL FUNCTION v1atu(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i


v1atu = 0.25*( v1(k,j,i) + v1(k,j+1,i) + v1(k,j,i-1) + v1(k,j+1,i-1) )

RETURN

END FUNCTION v1atu

!========================================================================

REAL FUNCTION v2atc(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i


v2atc = 0.5*( v2(k,j,i) + v2(k,j+1,i) )

RETURN

END FUNCTION v2atc

!=======================================================================

REAL FUNCTION v2atu(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i


v2atu = 0.25*( v2(k,j,i) + v2(k,j+1,i) + v2(k,j,i-1) + v2(k,j+1,i-1) )

RETURN

END FUNCTION v2atu

!=======================================================================

REAL FUNCTION w2atc(k,j,i)

INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i


w2atc = 0.5*( w2(k,j,i) + w2(k+1,j,i) )

RETURN

END FUNCTION w2atc

REAL FUNCTION cd(wind,tdif,drag,itype)

!***********************************************************************

!    *CD*           CALCULATE SURFACE DRAG COEFFICIENT

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE -  7 Dec 1994 @(COHERENS)metin.f 6.1

!       DESCRIPTION - OPTIONS SELECTED BY IDRAG (WIND DEPENDENCE) AND
!                     ITYPE (DEPENDENCE ON AIR-SEA TEMPERATURE DIFFERENCE)
!        IDRAG = 0 -> CONSTANT VALUE
!              = 1 -> LARGE AND POND (1981)
!              = 2 -> SMITH AND BANKE (1975)
!              = 3 -> GEERNAERT ET AL. (1986)
!              = 4 -> CHARNOCK'S RELATION
!        ITYPE = 0 -> NO DEPENDENCE ON TDIF
!              = 1 -> AS FUNCTION OF TDIF
!                 - USE BILINEAR INTERPOLATION IF ITYPE = 1

!       REFERENCE  - Section III-1.6.4 of the User Documentation
!                  - Geernaert G.L., 1990. Bulk parameterizations
!                    for the wind stress and heat fluxes.
!                    In: G.L. Geernaert and W.J. Plant (Editors),
!                    Surface Waves and Fluxes, Vol. 1 - Current theory.
!                    Kluwer Academic Publishers, Dordrecht, pp. 91-172.

!       CALLING PROGRAM - SURFLX

!       EXTERNALS - CHARNO

!***********************************************************************
!*    ARGUMENTS


REAL, INTENT(IN)                         :: wind
REAL, INTENT(IN)                         :: tdif
REAL, INTENT(IN)                         :: drag(13,11)
INTEGER, INTENT(IN)                  :: itype



!REAL :: charno

!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ITYPE*    INTEGER   SWITCH TO SELECT DEPENDENCE ON TDIF         (IN)
!    *DRAG*     REAL      ARRAY USED FOR BILINEAR EVALUATION OF DRAG
!                         COEFFICIENT                                 (IN)
!    *TDIF*     REAL      AIR-SEA TEMPERATURE DIFFERENCE              (IN)
!                         (USED ONLY IF ITYPE=1)
!    *WIND*     REAL      WIND SPEED AT 10 m HEIGHT   [m/s]            (IN)

!-----------------------------------------------------------------------

!*    LOCAL VARIABLES

INTEGER :: irow, jcol
REAL :: f11, f12, f21, f22, z1, z2


!     1. NEUTRAL VALUE
!     ----------------


IF (itype == 0) THEN
  IF (idrag == 0) THEN
    cd = 0.0013
  ELSE IF (idrag == 1) THEN
    IF (wind < 11.0) THEN
      cd = 0.0012
    ELSE
      cd = 0.001*(0.49+0.065*wind)
    END IF
  ELSE IF (idrag == 2) THEN
    cd = 0.001*(0.63+0.066*wind)
  ELSE IF (idrag == 3) THEN
    cd = 0.001*(0.43+0.097*wind)
  ELSE IF (idrag == 4) THEN
    cd = charno(wind)
  END IF
  
  
!     2. WITH DEPENDENCE ON TDIF
!     --------------------------
  
  
ELSE IF (itype == 1) THEN
  irow = MIN(12.0,0.4*wind+1)
  jcol = MIN(10.0,tdif+6)
  jcol = MAX0(1,jcol)
  f11 = drag(irow,jcol)
  f12 = drag(irow,jcol+1)
  f21 = drag(irow+1,jcol)
  f22 = drag(irow+1,jcol+1)
  z1 = 0.4*wind-irow+1
  z2 = tdif-jcol+6
  cd = (1.0-z1)*(1.0-z2)*f11+z2*(1.0-z1)*f12+z1*(1.0-z2)*f21 +z1*z2*f22
END IF

RETURN

END FUNCTION cd

!=======================================================================

REAL FUNCTION charno(wind)
!***********************************************************************

!    *CHARNO*      CALCULATE THE SURFACE DRAG COEFFICIENT USING
!                  CHARNOCK'S RELATION

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE -  6 May 1998 @(COHERENS)metin.f 8.4

!       DESCRIPTION - USE A NEWTON-RAPHSON ITERATION PROCEDURE

!       REFERENCE - Charnock H., 1955. Wind stress on a water surface.
!                   Quarterly Journal Royal Meteorological Society,
!                   81, 639-640.

!       CALLING PROGAM - CD

!***********************************************************************
!*    ARGUMENTS


REAL, INTENT(IN)                     :: wind


!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *WIND*     REAL      WIND SPEED AT 10m HEIGHT [m/s]              (IN)

!-----------------------------------------------------------------------

!*    LOCAL VARIABLES

REAL ::  c, df,  f, x, w10
REAL, PARAMETER :: alchar = 0.014
REAL, PARAMETER :: eps = 1.0E-06

!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ALCHAR*   REAL      CONSTANT IN CHARNOCK'S RELATION
!    *EPS*      REAL      PRECISION OF THE ALGORITHM

!-----------------------------------------------------------------------

c = LOG(10.0*g/alchar)
x = 0.001
w10 = MAX(0.1,wind)

10   f = LOG(x)+ckar/SQRT(x)+2.0*LOG(w10)-c
df = 1.0/x-0.5*ckar*x**(-1.5)
x = ABS(x-f/df)
IF (ABS(f/df) > eps) GO TO 10
charno = x

RETURN

END FUNCTION charno

!=======================================================================

REAL FUNCTION ce(wind,tdif,exch,itype)
!***********************************************************************

!    *CE*           CALCULATE THERMAL EXCHANGE COEFFICIENT

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE -  7 Dec 1994 @(COHERENS)metin.f 6.1

!       DESCRIPTION - OPTIONS SELECTED BY ITYPE
!        ITYPE = 0 -> NO DEPENDENCE ON TDIF
!              = 1 -> AS FUNCTION OF TDIF
!                 - USE BILINEAR INTERPOLATION IF ITYPE = 1

!       REFERENCE  - Section III-1.6.4 of the User Documentation
!                  - Geernaert G.L., 1990. Bulk parameterizations
!                    for the wind stress and heat fluxes.
!                    In: G.L. Geernaert and W.J. Plant (Editors),
!                    Surface Waves and Fluxes, Vol. 1 - Current theory.
!                    Kluwer Academic Publishers, Dordrecht, pp. 91-172.

!       CALLING PROGRAM - SURFLX

!       EXTERNALS -

!***********************************************************************

!*    ARGUMENTS


REAL, INTENT(IN)                         :: wind
REAL, INTENT(IN)                         :: tdif
REAL, INTENT(IN)                         :: exch(13,11)
INTEGER, INTENT(IN)                  :: itype




!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ITYPE*    INTEGER   SWITCH TO SELECT DEPENDENCE ON TDIF         (IN)
!    *EXCH*     REAL      ARRAY USED FOR BILINEAR EVALUATION OF
!                         THERMAL EXCHANGE COEFFICIENT                (IN)
!    *TDIF*     REAL      AIR-SEA TEMPERATURE DIFFERENCE              (IN)
!                         (USED ONLY IF ITYPE=1)
!    *WIND*     REAL      WIND SPEED AT 10m HEIGHT   [m/s]            (IN)

!-----------------------------------------------------------------------

!*    LOCAL VARIABLES

INTEGER :: irow, jcol
REAL :: f11, f12, f21, f22, z1, z2


!     1. NEUTRAL VALUE
!     ----------------


IF (itype == 0) THEN
  ce = 0.00113
  
  
!     2. WITH DEPENDENCE ON TDIF
!     --------------------------
  
  
ELSE IF (itype == 1) THEN
  irow = MIN(12.0,0.4*wind+1)
  jcol = MIN(10.0,tdif+6)
  jcol = MAX0(1,jcol)
  f11 = exch(irow,jcol)
  f12 = exch(irow,jcol+1)
  f21 = exch(irow+1,jcol)
  f22 = exch(irow+1,jcol+1)
  z1 = 0.4*wind-irow+1
  z2 = tdif-jcol+6
  ce = (1.0-z1)*(1.0-z2)*f11+z2*(1.0-z1)*f12+z1*(1.0-z2)*f21 +z1*z2*f22
END IF

RETURN

END FUNCTION ce

REAL FUNCTION fnlim(adv2,adv1,adv2up,adv1up,ifnlim)

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:45:01

!************************************************************************

!     *FNLIM*     EVALUATE WEIGHT FACTOR FOR TVD SCHEME

!       AUTHOR - KEVIN RUDDICK AND PATRICK LUYTEN

!       LAST UPDATE - 3 Jul 1995 @(COHERENS)Upwind.f 6.2

!       DESCRIPTION - CALCULATES EITHER SUPERBEE OR MONOTONIC LIMITER
!                     DEPENDING ON THE VALUE OF IFNLIM

!       REFERENCE - Section III-4.3.3a of the User Documentation

!       CALLING PROGRAM - HAD2DU, HAD2DV, XADVDIS, XHAD3DU, XHAD3DV,
!                         YADVDIS, YHAD3DU, YHAD3DV, ZADVDIS

!       EXTERNALS -

!************************************************************************


!*    ARGUMENTS


REAL, INTENT(IN)                         :: adv2
REAL, INTENT(IN)                         :: adv1
REAL, INTENT(IN)                         :: adv2up
REAL, INTENT(IN)                         :: adv1up
INTEGER, INTENT(IN)                  :: ifnlim




!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ADV2*     REAL      ADVECTIVE FLUX (2ND ORDER SCHEME)           (IN)
!    *ADV1*     REAL      ADVECTIVE FLUX (1ST ORDER SCHEME)           (IN)
!    *ADV2UP*   REAL      ADVECTIVE FLUX (2ND ORDER SCHEME) AT UPWIND POINT
!                                                                     (IN)
!    *ADV1UP*   REAL      ADVECTIVE FLUX (1ST ORDER SCHEME) AT UPWIND POINT
!                                                                     (IN)
!    *IFNLIM*   INTEGER   TYPE OF LIMITING FUNCTION                   (IN)
!                         (= 3 => SUPERBEE; = 4 => MONOTONIC)

!------------------------------------------------------------------------

!*    LOCAL VARIABLES

REAL :: adv21, adv21up,  r
REAL, PARAMETER :: epsmin=1.0E-12

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *R*        REAL      RATIO OF FLUXES
!    *ADV21*    REAL      DIFF OF FLUXES (2ND-1ST ORDER)
!    *ADV21UP*  REAL      DIFF OF FLUXES (2ND-1ST ORDER) AT UPWIND POINT
!------------------------------------------------------------------------

adv21   = adv2   - adv1
adv21up = adv2up - adv1up

IF (ABS(adv21) > epsmin) THEN
  r=adv21up/adv21
  IF (ifnlim == 3) THEN
    fnlim=MAX(0.0,MIN(2.0*r,1.0),MIN(r,2.0))
  ELSE IF (ifnlim == 4) THEN
    fnlim=(r+ABS(r))/(1.0+ABS(r))
  END IF
ELSE
  fnlim=0.0
END IF

RETURN

END FUNCTION fnlim

REAL FUNCTION random(iflag)
!***********************************************************************

!    *RANDOM*     RANDOM GENERATOR RETURNING A RANDOM NUMBER BETWEEN
!                 0 AND 1

!       AUTHOR - JAN ROSS

!       LAST UPDATE - 11 Jun 1998 @(COHERENS)Util.f 8.4

!       DESCRIPTION - RETURNS RANDOM NUMBER BETWEEN 0 AND 1
!                   - USES CHAOS EQUATION X=R*(X-X*2) WITH ITERATIONS
!                     OF X BETWEEN 0 AND 1
!                   - RANDOM GENERATOR BECOMES INITIALISED WHEN IFLAG=0
!                   - RANDOM IS USED IN THE PARTICLE MODULE

!       REFERENCE -

!       CALLING PROGRAM - MCOSB, MCROB, MC3D, SEDLAG

!       EXTERNALS -

!***********************************************************************

!*    ARGUMENTS


INTEGER, INTENT(IN)                  :: iflag


!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *IFLAG*    INTEGER   INITIALISES THE RANDOM GENERATOR IF = 0     (IN)

!------------------------------------------------------------------------

!*    LOCAL VARIABLES

DOUBLE PRECISION ::    x
DOUBLE PRECISION, PARAMETER :: dpi=3.141592741D0
DOUBLE PRECISION, PARAMETER :: r=3.999999D0
DOUBLE PRECISION, PARAMETER :: rr=r*5.d0
INTEGER :: i, iend

SAVE x


IF (iflag /= 0) THEN
  x = 0.49D0
END IF

iend = INT(rr-x*5.d0)
DO  i=1,iend
  x = r*(x-x*x)
END DO
random = REAL(ACOS((x-0.5D0)*2.d0)/dpi)


RETURN

END FUNCTION random

!=======================================================================

REAL FUNCTION ran3(iseed)
!***********************************************************************

!    *RAN3*       RANDOM GENERATOR RETURNING A RANDOM NUMBER BETWEEN
!                 0 AND 1

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 16 Mar 1998 @(COHERENS)Util.f 8.3

!       DESCRIPTION - RETURNS RANDOM NUMBER BETWEEN 0 AND 1
!                   - SET ISEED TO A NEGATIVE VALUE TO (RE)INITIALISE
!                     THE SEQUENCE

!        REFERENCE - Press W.H., Flannery B.P., Teukolsky S.A. and Vetterling
!                    W.T., 1989. Numerical Recipes. The art of scientific
!                    computing. Cambridge University Press, Cambridge, 702 pp.

!       CALLING PROGRAM - RAN3 IS USED TO SET UP THE WIND DIRECTION IN THE
!                         FORCING OF THE CSBIO TEST CASE

!       EXTERNALS -

!***********************************************************************

!*    ARGUMENTS


INTEGER, INTENT(IN OUT)                     :: iseed


!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ISEED*    INTEGER   (RE)INITIALISES THE RANDOM GENERATOR IF < 0 (IN)

!------------------------------------------------------------------------

INTEGER :: i, iff, ii, inext, inextp, k
INTEGER ::  mj, mk
INTEGER :: ma(55)

INTEGER, PARAMETER :: mbig=1000000000
INTEGER, PARAMETER :: mseed=161803398
INTEGER, PARAMETER :: mz=0
REAL, PARAMETER :: fac=1.0/mbig

SAVE inext, inextp, ma
DATA iff /0/


IF (iseed < 0.OR.iff == 0) THEN
  iff = 1
  mj = mseed - IABS(iseed)
  mj = MOD(mj,mbig)
  ma(55) = mj
  mk = 1
  DO  i=1,54
    ii = MOD(21*i,55)
    ma(ii) = mk
    mk = mj - mk
    IF (mk < mz) mk = mk + mbig
    mj = ma(ii)
  END DO
  DO  k=1,4
    DO  i=1,55
      ma(i) = ma(i) - ma(1+MOD(i+30,55))
      IF (ma(i) < mz) ma(i) = ma(i) + mbig
    END DO
  END DO
  inext = 0
  inextp = 31
  iseed = 1
END IF
inext = inext + 1
IF (inext == 56) inext = 1
inextp = inextp + 1
IF (inextp == 56) inextp = 1
mj = ma(inext)-ma(inextp)
IF (mj < mz) mj = mj+mbig
ma(inext) = mj
ran3 = mj*fac


RETURN

END FUNCTION ran3

!=======================================================================

REAL FUNCTION buofr2(k,j,i)
!***********************************************************************

!    *BUOFR2*   CALCULATE SQUARED BUOYANCY FREQUENCY [1/s**2] AT W-NODES

!     AUTHOR - PATRICK LUYTEN AND KEVIN RUDDICK

!     LAST UPDATE - 3 Jul 1995          @(COHERENS)veddy_in.f 6.2

!     DESCRIPTION - EVALUATE THE SQUARED BUOYANCY (BRUNT-VAISALA)
!                   FREQUENCY AT W NODES
!                 - ROUND-OFF ERROR IS REDUCED BY AVOIDING SREF AND
!                   TREF IN CALCULATION
!                 - BUOFR2 IS HORIZONTALLY SMOOTHED TO AVOID SPURIOUS
!                   HORIZONTAL OSCILLATIONS BY THE TURBULENCE SCHEME

!     REFERENCE - Section III-4.5.1 of the User Documentation

!     CALLING PROGRAM - DISMIN, FNTURF, INCTUR, PRODUC, RICH, VEDDY2, ZLMAX

!     EXTERNALS -

!***********************************************************************
!*    ARGUMENTS


INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                      :: i


!*    LOCAL VARIABLES

REAL :: buofrc, buofre, buofrn, buofrs, buofrw
REAL :: gradsc, gradse, gradsn, gradss, gradsw
REAL :: gradtc, gradte, gradtn, gradts, gradtw
REAL :: tbetc,  tbete,  tbetn,  tbets,  tbetw
REAL :: sbetc,  sbete,  sbetn,  sbets,  sbetw
!REAL :: gz2w


sbetc = 0.5*(sbet(k,j,i)+sbet(k-1,j,i))
tbetc = 0.5*(tbet(k,j,i)+tbet(k-1,j,i))
gradsc = (s(k,j,i)-s(k-1,j,i))/gz2w(k,j,i)
gradtc = (tr(k,j,i)-tr(k-1,j,i))/gz2w(k,j,i)
buofrc = g*(-sbetc*gradsc+tbetc*gradtc)

IF ((j > 1).AND.(npiy(j,i) == 1)) THEN
  sbets = 0.5*(sbet(k,j-1,i)+sbet(k-1,j-1,i))
  tbets = 0.5*(tbet(k,j-1,i)+tbet(k-1,j-1,i))
  gradss = (s(k,j-1,i)-s(k-1,j-1,i))/gz2w(k,j-1,i)
  gradts = (tr(k,j-1,i)-tr(k-1,j-1,i))/gz2w(k,j-1,i)
  buofrs = g*(-sbets*gradss+tbets*gradts)
ELSE
  buofrs = buofrc
END IF

IF ((j < nr).AND.(npiy(j+1,i) == 1)) THEN
  sbetn = 0.5*(sbet(k,j+1,i)+sbet(k-1,j+1,i))
  tbetn = 0.5*(tbet(k,j+1,i)+tbet(k-1,j+1,i))
  gradsn = (s(k,j+1,i)-s(k-1,j+1,i))/gz2w(k,j+1,i)
  gradtn = (tr(k,j+1,i)-tr(k-1,j+1,i))/gz2w(k,j+1,i)
  buofrn = g*(-sbetn*gradsn+tbetn*gradtn)
ELSE
  buofrn = buofrc
END IF

IF ((i > 1).AND.(npix(j,i) == 1)) THEN
  sbetw = 0.5*(sbet(k,j,i-1)+sbet(k-1,j,i-1))
  tbetw = 0.5*(tbet(k,j,i-1)+tbet(k-1,j,i-1))
  gradsw = (s(k,j,i-1)-s(k-1,j,i-1))/gz2w(k,j,i-1)
  gradtw = (tr(k,j,i-1)-tr(k-1,j,i-1))/gz2w(k,j,i-1)
  buofrw = g*(-sbetw*gradsw+tbetw*gradtw)
ELSE
  buofrw = buofrc
END IF

IF ((i < nc).AND.(npix(j,i+1) == 1)) THEN
  sbete = 0.5*(sbet(k,j,i+1)+sbet(k-1,j,i+1))
  tbete = 0.5*(tbet(k,j,i+1)+tbet(k-1,j,i+1))
  gradse = (s(k,j,i+1)-s(k-1,j,i+1))/gz2w(k,j,i+1)
  gradte = (tr(k,j,i+1)-tr(k-1,j,i+1))/gz2w(k,j,i+1)
  buofre = g*(-sbete*gradse+tbete*gradte)
ELSE
  buofre = buofrc
END IF

buofr2 = (2.0*buofrc+buofrs+buofrn+buofrw+buofre)/6.0


RETURN

END FUNCTION buofr2

!====================================================================

REAL FUNCTION shfr2(k,j,i)
!***********************************************************************

!    *SHFR2*    CALCULATE THE SQUARED SHEAR FREQUENCY [1/s**2] AT W-NODES

!     AUTHOR - PATRICK LUYTEN AND KEVIN RUDDICK

!     LAST UPDATE - 7 Dec 1994          @(COHERENS)veddy_in.f 6.1

!     DESCRIPTION - SHFR2 WILL BE REDUCED BY HALF NEAR SOLID (BUT NOT
!                   OPEN) BOUNDARIES BECAUSE OF THE TWO-SIDED INTERPOLATION

!     REFERENCE - Section III-4.5.1 of the User Documentation

!     CALLING PROGRAM - FNTURF, PRODUC, RICH

!     EXTERNALS -

!***********************************************************************
!*    ARGUMENTS


INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                      :: i


!*    LOCAL VARIABLES

REAL :: dzu2, dzv2
!REAL :: gz2w

!     ---bracketing to reduce round-off error
dzu2 = (u2(k,j,i)-u2(k-1,j,i))+(u2(k,j,i+1)-u2(k-1,j,i+1))
dzv2 = (v2(k,j,i)-v2(k-1,j,i))+(v2(k,j+1,i)-v2(k-1,j+1,i))
shfr2 = (dzu2**2+dzv2**2)*((0.5/gz2w(k,j,i))**2)

RETURN

END FUNCTION shfr2
!===========================================================================

REAL FUNCTION fnturf(k,j,i)
!***********************************************************************

!    *FNTURF*   EVALUATE K/L**2 [1/s**2] IN THE ZERO-EQUATION TURBULENCE SCHEME

!     AUTHOR - PATRICK LUYTEN

!     LAST UPDATE - 31 Mar 1998          @(COHERENS)veddy_in.f 8.4

!     DESCRIPTION - =TKEW(I,J,K)/ZLW(I,J,K)**2

!     REFERENCE  - Section III-1.2.2b of the User Documentation
!                - Press W.H., Flannery B.P., Teukolsky S.A. and Vetterling
!                  W.T., 1989. Numerical Recipes. The art of scientific
!                  computing. Cambridge University Press, Cambridge, 702 pp.

!     CALLING PROGRAM - VEDDY2

!     EXTERNALS - BUOFR2, FVEDMA, GVEDMA, SHFR2

!************************************************************************
!*    ARGUMENTS


INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i



!*    LOCAL VARIABLES

REAL :: b0, c0, disc, root1, root2, turf
!REAL :: buofr2, fvedma, gvedma, shfr2


IF (istpar == 1) THEN
  b0 = clev2(1)*buofr2(k,j,i)+clev2(2)*shfr2(k,j,i)
  c0 = buofr2(k,j,i)*(clev2(3)*shfr2(k,j,i) + clev2(4)*buofr2(k,j,i))
  disc = b0**2 - 4*c0
  IF (disc <= 0.0) THEN
    fnturf = 0.0
  ELSE
    root1 = -0.5*(b0+SIGN(1.0,b0)*SQRT(disc))
    IF (root1 /= 0.0) THEN
      root2 = c0/root1
    ELSE
      root2 = 0.0
    END IF
    turf = MAX(root1,root2)
    fnturf = MAX(turf,0.0)
  END IF
  
ELSE IF (istpar == 2) THEN
  fnturf = (shfr2(k,j,i)*clev25(1)*fvedma(k,j,i)  &
      - buofr2(k,j,i)*clev25(5)*gvedma(k,j,i))
END IF

IF ((itcpar == 1).AND.(istpar == 2)) fnturf = fnturf/eps0
IF (itcpar == 2) fnturf = fnturf/eps0**2


RETURN

END FUNCTION fnturf

!=========================================================================

REAL FUNCTION fnsmu(k,j,i,x)
!*************************************************************************

!    *FNSMU*    EVALUATE STABILITY COEFFICIENT FOR MOMENTUM

!     AUTHOR - PATRICK LUYTEN

!     LAST UPDATE - 31 Mar 1998          @(COHERENS)veddy_in.f 8.4

!     DESCRIPTION - EVALUATES STABILITY COEFFICIENT SMU (S_m IN K-L,
!                   S_u IN K-EPS THEORY) AS FUNCTION OF X
!                   (G_h IN K-L, ALPHA_N IN K-EPS THEORY)
!                   X>0 STABLE ; <0 UNSTABLE STRATIFICATION
!                 - FORM SELECTED BY ISTPAR
!        = 1 => AS FUNCTION OF X (G_h IN K-L, ALPHA_N IN K-EPS THEORY)
!        = 2 => AS FUNCTION OF RICHARDSON NUMBER

!     REFERENCE -  Section III-1.22a of the User Documentation

!     CALLING PROGRAM - INCTUR, VEDDY2

!     EXTERNALS - FVEDMA

!**************************************************************************
!*    ARGUMENTS


INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i
REAL, INTENT(IN)                         :: x



!*    LOCAL VARIABLES

REAL :: smubot, smutop
!REAL :: fvedma


IF (istpar == 1) THEN
  smutop = clev25(1) + clev25(2)*x
  smubot = 1.0 + clev25(3)*x + clev25(4)*x*x
  fnsmu = smutop/smubot
ELSE IF (istpar == 2) THEN
  fnsmu = clev25(1)*fvedma(k,j,i)
END IF


RETURN

END FUNCTION fnsmu

!==========================================================================

REAL FUNCTION fnshb(k,j,i,x)
!**************************************************************************

!    *FNSHB*    EVALUATE STABILITY COEFFICIENT FOR SCALARS

!     AUTHOR - PATRICK LUYTEN

!     LAST UPDATE - 31 Mar 1998           @{COHERENS}veddy_in.f  8.4

!     DESCRIPTION - EVALUATES STABILITY COEFFICIENT SHB (SQRT(2)*S_h IN K-L,
!                   S_b IN K-EPS THEORY) AS FUNCTION OF X
!                   (G_h IN K-L, ALPHA_N IN K-EPS THEORY)
!                   X>0 STABLE ; <0 UNSTABLE STRATIFICATION
!                 - FORM SELECTED BY ISTPAR
!        = 1 => AS FUNCTION OF X (=G_h IN K-L, ALPHA_N IN K-EPS THEORY)
!        = 2 => AS FUNCTION OF RICHARDSON NUMBER

!     REFERENCE - Section III-1.2.2a of the User Documentation

!     CALLING PROGRAM - VEDDY2

!     EXTERNALS - GVEDMA

!**************************************************************************
!*    ARGUMENTS


INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i
REAL, INTENT(IN)                     :: x


!REAL :: gvedma


IF (istpar == 1) THEN
  fnshb = clev25(5)/(1.0 + clev25(6)*x)
ELSE IF (istpar == 2) THEN
  fnshb = clev25(5)*gvedma(k,j,i)
END IF


RETURN

END FUNCTION fnshb

!=======================================================================

REAL FUNCTION rich(k,j,i)
!***********************************************************************

!    *RICH*     EVALUATE RICHARDSON NUMBER

!     AUTHOR - PATRICK LUYTEN AND KEVIN RUDDICK

!     LAST UPDATE - 31 Mar 1998          @(COHERENS)veddy_in.f 8.4

!     DESCRIPTION -

!     REFERENCE -

!     CALLING PROGRAM - FVEDMA, GVEDMA, VEDDY1

!     EXTERNALS - BUOFR2, SHFR2

!***********************************************************************
!*    ARGUMENTS


INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i



!*    LOCAL VARIABLES

REAL ::  rbot,  rtop
!REAL :: buofr2, shfr2
REAL, PARAMETER :: epsmin = 1.0E-20
REAL, PARAMETER :: ricmax = 1.0E+03


IF (k == 1.OR.k == nz+1) THEN
  rich = 0.0
ELSE
  rtop = buofr2(k,j,i)
  rbot = shfr2(k,j,i)
  IF ((ABS(rtop) > rbot*ricmax).OR.(rbot < epsmin)) THEN
    rich = ricmax*SIGN(1.0,rbot)
  ELSE
    rich = rtop/rbot
  END IF
END IF


RETURN

END FUNCTION rich

!=======================================================================

REAL FUNCTION fvedma(k,j,i)
!***********************************************************************

!    *FVEDMA*   EVALUATE THE FACTOR f_m(Ri) IN THE MUNK-ANDERSON FORMULATION

!     AUTHOR - PATRICK LUYTEN

!     LAST UPDATE - 31 Mar 1998          @(COHERENS)veddy_in.f 8.4

!     DESCRIPTION -

!     REFERENCE - Section III-1.2.1a of the User Documentation

!     CALLING PROGRAM - FNSMU, FNTURF, VEDDY1

!     EXTERNALS - RICH

!***********************************************************************


!*    ARGUMENTS


INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i



!*    LOCAL VARIABLES

REAL :: denom, epsmin
!REAL :: rich


epsmin = rmaxnut**(-1.0/amn1)
denom = MAX(1.0+ama*rich(k,j,i),epsmin)
fvedma = 1.0/(denom**amn1)


RETURN

END FUNCTION fvedma

!=======================================================================

REAL FUNCTION gvedma(k,j,i)
!***********************************************************************

!    *GVEDMA*   EVALUATE THE FACTOR g_m(Ri) IN THE MUNK-ANDERSON FORMULATION

!     AUTHOR - PATRICK LUYTEN

!     LAST UPDATE - 31 Mar 1998          @(COHERENS)veddy_in.f 8.4

!     DESCRIPTION -

!     REFERENCE - Section III-1.2.1a of the User Documentation

!     CALLING PROGRAM - FNSHB, FNTURF, VEDDY1

!     EXTERNALS - RICH

!***********************************************************************

!*    ARGUMENTS


INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i



!*    LOCAL VARIABLES

REAL :: denom, epsmin
!REAL :: rich


epsmin = rmaxlat**(-1.0/amn2)
denom = MAX(1.0+amb*rich(k,j,i),epsmin)
gvedma = 1.0/(denom**amn2)


RETURN

END FUNCTION gvedma

!=======================================================================

REAL FUNCTION zlmax(k,j,i)
!***********************************************************************

!    *ZLMAX*    EVALUATE MAXIMUM MIXING LENGTH (K-EPS THEORY)

!     AUTHOR - PATRICK LUYTEN

!     LAST UPDATE - 7 Dec 1994          @(COHERENS)veddy_in.f 6.1

!     DESCRIPTION - CALCULATES UPPER LIMIT FOR TURBULENCE LENGTH SCALE [m]
!                   IN THE CASE OF STABLE STRATIFICATION

!     REFERENCE - Section III-1.2.2c of the User Documentation

!     CALLING PROGRAM - VEDDY2

!     EXTERNALS - BUOFR2

!***********************************************************************
!*    ARGUMENTS


INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i


!*    LOCAL VARIABLES

REAL :: buofr
!REAL :: buofr2


buofr = buofr2(k,j,i)
IF (buofr > 0.0) THEN
  zlmax = SQRT(gamax*tkew(k,j,i)/buofr)
ELSE
  zlmax = 1.0E+10
END IF

RETURN

END FUNCTION zlmax


!========================================================================

REAL FUNCTION dismin(k,j,i)
!************************************************************************

!    *DISMIN*   EVALUATE MINIMUM DISSIPATION RATE (K-EPS THEORY)

!     AUTHOR - PATRICK LUYTEN

!     LAST UPDATE - 7 Dec 1994          @(COHERENS)veddy_in.f 6.1

!     DESCRIPTION - CALCULATES LOWER LIMIT FOR DISSIPATION RATE [W/kg]
!                   IN THE CASE OF STABLE STRATIFICATION

!     REFERENCE - Section III-1.2.2c of the User Documentation

!     CALLING PROGRAM - VEDDY2

!     EXTERNALS - BUOFR2

!*************************************************************************
!*    ARGUMENTS


INTEGER, INTENT(IN)                  :: k
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                  :: i


!*    LOCAL VARIABLES

REAL :: buofr
!REAL :: buofr2


buofr = buofr2(k,j,i)
IF (buofr > 0.0) THEN
  dismin = tkew(k,j,i)*SQRT(buofr/gamax)
ELSE
  dismin = 0.0
END IF

RETURN

END FUNCTION dismin

END MODULE functions
