MODULE cohrun
USE param
USE functions
USE cohini

CONTAINS

SUBROUTINE RUNCOH(NNT)

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
INTEGER, INTENT(IN) :: NNT
INTEGER :: NT,K,J,I,N,L,nend,NMAX
REAL :: ad2,etac,zzz,ztop,zbot,eee,tcut

nt = NNT
  
 ttime = real(nt)*delt/(24.*3600.)

 WRITE(6,*)'HENRY time/days = ',ttime

   ad = MIN(ttime/2.0,1.0)
 !  IF(ttime>1.0)ad = 0.0

ugeo = 0.0 !!!!-0.1*ad

! wind forcing
DO  i=1,nc
DO  j=1,nr
     windu2(j,i) = 0.0 !-ad*10.0
END DO
END DO

CALL surflx

etam =  0.0001

  FLOAT_ON = .true.
!
! CALL metin
!
!     2.3 SOLVE TURBULENCE AND 3D-MOMEMTUM EQUATIONS (PREDICTOR STEP)
!     ---------------------------------------------------------------

  IF (iopt3 == 1) THEN
    IF ((nt-1)/ic3d*ic3d == (nt-1)) THEN
      IF (iodif > 0) CALL heddy
      IF (igrdim == 1) CALL densty1
      IF ((ioptd > 0).AND.(igrdim == 3)) CALL densty
      IF (ioptk == 1) CALL veddy1
      IF (ioptk == 2) CALL veddy2
       CALL crrnt3p
    END IF
  END IF

!     2.4 SOLVE 2D-MOMENTUM EQUATIONS
!     -------------------------------
  IF (iopt2 == 1) THEN
    CALL contny
 !
CALL surflx

do i = 1,nc
  zeta2(1,i) = 0.0 
!  zeta2(nr-1,i) = -ad*0.05
!  zeta2(nr,i) = zeta2(nr-1,i)
end do

!do j = 1,nr
!  zeta2(j,1) = -ad*REAL(j-1)/REAL(nr-1)*0.05
!end do

   CALL crrnt2

  END IF

!     2.5 UPDATE CURRENTS (CORRECTOR STEP)
!     ------------------------------------
  
  IF (iopt3 == 1) THEN
    IF (nt/ic3d*ic3d == nt) THEN
      CALL transv
      IF (igrdim == 1) THEN
        CALL crrnt1c
      ELSE
        CALL crrnt3c
      END IF

      IF ((iadvc /= 0).OR.(iadvs /= 0)) CALL wcalc
    END IF
  END IF
!
!       2.6 SOLVE SALT AND TEMPERATURE EQUATIONS
!       ----------------------------------------
  
  IF (ioptd > 0) THEN
    IF (nt/ic3d*ic3d == nt) THEN
  
      IF (ioptsa > 0)  CALL salt  
      IF (iopthe > 0) CALL heat
      CALL euler

IF(ttime>4.0)THEN
fcount = fcount+1
DO  i=1,nc
DO  j=1,nr
  zetafinal(j,i) =  zetafinal(j,i)+zeta2(j,i)
  DO k=1,nz
    ufinal(k,j,i) = ufinal(k,j,i)+u2(k,j,i)
    vfinal(k,j,i) = vfinal(k,j,i)+v2(k,j,i)
    wfinal(k,j,i) = wfinal(k,j,i)+w2(k,j,i)
    tfinal(k,j,i) = tfinal(k,j,i)+t(k,j,i)
    sfinal(k,j,i) = sfinal(k,j,i)+s(k,j,i)
    rfinal(k,j,i) = rfinal(k,j,i)+ro(k,j,i)-r0ref
    pfinal(k,j,i) = pfinal(k,j,i)+rpress(k,j,i)
    hedfinal(k,j,i) = hedfinal(k,j,i)+heddydc(k,j,i)
    vedfinal(k,j,i) = vedfinal(k,j,i)+veddyd(k,j,i)
  END DO
  wfinal(nz+1,j,i) = wfinal(nz+1,j,i)+w2(nz+1,j,i)
END DO
END DO
END IF

!======================
! float prediction

IF(ttime>=2.0)THEN
 tcut = MIN((ttime-2.0)/5.0,1.0)*2000
 NMAX = MAX(tcut,1)
 ! write(6,*)NMAX,"<============NMAX"
 ! pause
 DO N = 1,NMAX
  trac_start(N) = .TRUE.
  IF(trac_stop(N)) trac_start(N) = .FALSE.
  trac_start(N+2000) = .TRUE.
  IF(trac_stop(N+2000)) trac_start(N+2000) = .FALSE.
 END DO
     CALL MC3D
END IF

    END IF
  END IF

!     2.7 RECOMPUTE DENSITY RELATED PARAMETERS
!     ----------------------------------------
  
  IF (ioptd > 0) THEN
    IF (nt/ic3d*ic3d == nt) CALL searho
  END IF


END SUBROUTINE RUNCOH


!=======================================================================

SUBROUTINE densty1
!************************************************************************

!    *DENSTY1*     CALCULATE PRESSURE GRADIENT (1-D APPLICATION)

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 5 May 1999 @(COHERENS)densty.f 8.4

!       DESCRIPTION - EVALUATES PRESSURE GRADIENT (UQDEN/VQDEN) AS
!                     -G* SURFACE SLOPE + BAROCLINIC COMPONENT

!       REFERENCE - Section III-1.1.4 of the User Documentation

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS -

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, icon, ind, j, k
REAL :: forceu, forcev, phase, zet

!**** SOME VARIABLES DEFINED IN bounds.inc AND crrnts.inc NOW HAVE A ****
!****              DIFFERENT MEANING                                 ****

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *IND*      INTEGER   INDEX WHICH EQUALS 1 FOR 1-D APPLICATIONS AND IS
!                         ZERO OTHERWISE TO PREVENT (EVENTUAL) ERRORS
!                         DURING COMPILATION
!    *FORCEU*   REAL      PRESSURE GRADIENT IN X-DIRECTION [m/s2]
!    *FORCEV*   REAL      PRESSURE GRADIENT IN Y-DIRECTION [m/s2]
!    *AMPOBU*   REAL      (0,*): AMPLITUDES OF SURFACE ELEVATION      [m]
!                         (1,0): X-COMP. OF BAROCLINIC PRESSURE GRADIENT
!                         (1,*): AMPLITUDES OF PRESSURE GRAD IN X-DIR [m/s2]
!    *AMPOBV*   REAL      (1,0): Y-COMP. OF BAROCLINIC PRESSURE GRADIENT
!                         (1,*): AMPLITUDES OF PRESSURE GRAD IN Y-DIR [m/s2]
!    *PHAOBU*   REAL      (0,*): PHASES OF SURFACE ELEVATION          [rad]
!                         (1,*): PHASES OF PRESSURE GRADIENT IN X-DIR [rad]
!    *PHAOBV*   REAL      (1,*): PHASES OF PRESSURE GRADIENT IN Y-DIR [rad]
!    *UQDEN*    REAL      PRESSURE GRADIENT IN X-DIRECTION [m/s2]
!    *VQDEN*    REAL      PRESSURE GRADIENT IN Y-DIRECTION [m/s2]

!------------------------------------------------------------------------

!     1. UPDATE PHASE
!     ---------------

DO  icon=1,ncon
  phase0(icon) = MOD(phase0(icon)+ del3*sigma(icon),2.0*pi)
END DO


!     2. UPDATE SURFACE ELEVATION AND GRID SPACING
!     --------------------------------------------

zet = 0.0
DO  icon=1,ncon
  phase = phase0(icon) - phaobu(0,icon)
  zet = zet + ampobu(0,icon)*COS(phase)
END DO
DO  i=1,nc
  DO  j=1,nr
    zeta2(j,i) = zet
  END DO
END DO

CALL transv


!     3. UPDATE PRESSURE GRADIENT
!     ---------------------------

IF (igrdim == 1) THEN
  ind = 1
ELSE
  ind = 0
END IF


!     3.1 U-COMPONENT
!     ---------------

forceu = -ampobu(ind,0)
DO  icon=1,ncon
  phase = phase0(icon) - phaobu(ind,icon)
  forceu = forceu-ampobu(ind,icon)*COS(phase)
END DO
DO  i=2,nc
  DO  j=1,nr
    DO  k=1,nz
      uqden(k,j,i) = forceu
    END DO
  END DO
END DO

!     3.2 V-COMPONENT
!     ---------------

forcev = -ampobv(ind,0)
DO  icon=1,ncon
  phase = phase0(icon) - phaobv(ind,icon)
  forcev = forcev-ampobv(ind,icon)*COS(phase)
END DO
DO  i=1,nc
  DO  j=2,nr
    DO  k=1,nz
      vqden(k,j,i) = forcev
    END DO
  END DO
END DO


RETURN

END SUBROUTINE densty1


!

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:35:57


SUBROUTINE xadvdis (phi,hedphi,uadv,gz,nzphi,ivpu,  &
    phivp,jadvs,jdifs,simp,sexp)
!************************************************************************

!    *XADVDIS*    EVALUATE THE ADVECTIVE AND DIFFUSION TERMS ALONG THE
!                 X-DIRECTION IN THE TRANSPORT EQUATION FOR A SCALAR
!                 QUANTITY (CELL-CENTERED OR AT W-NODES)

!       AUTHOR - KEVIN RUDDICK AND PATRICK LUYTEN

!       LAST UPDATE - 26 Mar 1999 @(COHERENS)advdis.f 8.4

!       DESCRIPTION - COMPUTES X-DIRECTION HORIZONTAL ADVECTION/DIFFUSION
!                     OF A 3-D SCALAR QUANTITY HELD AT THE CELL CENTRE
!                     (E.G. CONCENTRATION OF CONTAMINANT, SUSPENDED PARTICLES,
!                     SALINITY) OR AT W-NODES (TURBULENCE VARIABLE)
!                   - INCORPORATES OPEN BOUNDARY CONDITIONS
!                   - ADVECTIVE/DIFFUSIVE TERMS ARE EVALUATED EXPLICITLY
!                   - ADVECTION SCHEME SELECTED BY JADVS
!                      = 1 => FIRST ORDER UPWIND
!                      = 2 => LAX-WENDROFF
!                      = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                      = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!                   - DIFFUSION SELECTED BY JDIFS
!                      = 0 => OFF
!                      = 1 => UNIFORM
!                      = 2 => SMAGORINSKY
!                   - FLUXES ARE CALCULATED AT U- OR X-NODES

!       REFERENCE - Sections III-4.4.2-4.4.3 of the User Documentation

!       CALLING PROGRAM - TRANSPC, TRANSPW

!       EXTERNALS - FNLIM, ZERO

!************************************************************************
!*    ARGUMENTS

INTEGER, INTENT(IN)                      :: nzphi
REAL, INTENT(IN)                         :: phi(nzphi,nr,nc)
REAL, INTENT(IN)                         :: hedphi(nzphi,nr,nc+1)
REAL, INTENT(IN)                         :: uadv(nzphi,nr,nc+1)
REAL, INTENT(IN)                         :: gz(nzphi,nr,nc)
INTEGER, INTENT(IN)                      :: ivpu(0:nobu)
REAL, INTENT(IN)                         :: phivp(0:nzphi,0:nvprof)
INTEGER, INTENT(IN)                      :: jadvs
INTEGER, INTENT(IN)                      :: jdifs
REAL, INTENT(IN OUT)                        :: simp(nzphi,nr,nc)
REAL, INTENT(IN OUT)                        :: sexp(nzphi,nr,nc)

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *PHI*      REAL      SCALAR TO DE ADVECTED/DIFFUSED [PHI]        (IN)
!    *HEDPHI*   REAL      HORIZONTAL DIFFUSION COEFFICIENT [m2/s]     (IN)
!    *UADV*     REAL      ADVECTING VELOCITY [m/s]                    (IN)
!    *GZ*       REAL      VERTICAL GRID SPACING ABOUT PHI             (IN)
!    *NZPHI*    INT       VERTICAL DIMENSION OF PHI (=NZ OR NZ+1)     (IN)
!    *IVPU*     INT       NO. OF PROFILE AT U-O.B.                    (IN)
!    *PHIVP*    REAL      OPEN BOUNDARY CONDITION [PHI]               (IN)
!                       PHIVP(K,*)>PHIVP(0,*)  : IMPOSED INFLOW
!                       PHIVP(K,*)<= PHIVP(0,*): NORMAL GRADIENT CONDITION
!    *JADVS*    INTEGER   SWITCH TO SELECT TYPE OF ADVECTION SCHEME (0,1,2,3)
!                                                                     (IN)
!    *JDIFS*    INTEGER   SWITCH TO SELECT TYPE OF HORIZONTAL DIFFUSION
!                         (0,1,2)                                     (IN)
!    *SIMP*     REAL      LEFT HAND SIDE OF PHI-EQUATION /PHI(I,J,K)
!                                                                 (IN/OUT)
!    *SEXP*     REAL      RIGHT HAND SIDE OF PHI-EQUATION  [PHI]
!                                                                 (IN/OUT)

!------------------------------------------------------------------------

!*    LOCAL VARIABLES

INTEGER :: i, ii, ivp, j, k
REAL :: advfl(nz,nr,nc+1), advflw(nz,nr,nc+1)
REAL :: advfup(nz,nr,nc+1), diffl(nz,nr,nc+1)
REAL :: utemp(nz,nr,nc+1),htemp(nz,nr,nc+1),uadv2(nz,nr,nc+1)
REAL :: cfl, difc, psi
REAL :: uadvu, uwt
!REAL :: fnlim, gx2u

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ADVFUP*   REAL      UPWIND ADVECTIVE FLUX OF PHI IN X-DIRECTION AT U-NODE
!                                                                [PHI*m/s]
!    *ADVFLW*   REAL      LW ADVECTIVE FLUX OF PHI IN X-DIRECTION AT U-NODE
!                                                                [PHI*m/s]

!    *ADVFL*    REAL      TVD ADVECTIVE FLUX OF PHI IN X-DIRECTION AT U-NODE
!                                                                [PHI*m/s]
!    *DIFFL*    REAL      DIFFUSIVE FLUX OF PHI IN X-DIRECTION AT U-NODE
!                                                                [PHI*m/s]
!    *CFL*      REAL      COURANT-FRIEDRICH-LEVY NUMBER
!    *PSI*      REAL      FLUX LIMITER
!    *UWT*      REAL      SIGN OF ADVECTIVE VELOCITY

!------------------------------------------------------------------------

!     1. INITIALISATION OF TEMPORARY ARRAYS
!     -------------------------------------


CALL zero33(advfup,nz,nr,nc+1)
CALL zero33(advflw,nz,nr,nc+1)
CALL zero33(advfl,nz,nr,nc+1)
CALL zero33(diffl,nz,nr,nc+1)

DO  i=1,nc+1
DO  j=1,nr
DO  k=1,nz
  uadv2(k,j,i) = uadv(k,j,i)+ugeo
END DO
END DO
END DO


DO  i=1,nc+1
DO  j=1,nr
DO  k=1,nz
  utemp(k,j,i) = uadv2(k,j,i)
  htemp(k,j,i) = hedphi(k,j,i)
END DO
END DO
END DO

IF(trigger==1.0)CALL zero33(utemp,nz,nr,nc+1)
IF(trigger==2.0)call zero33(htemp,nz,nr,nc+1)

!     2. ADVECTIVE FLUXES AT U-NODES
!     ------------------------------

IF (jadvs > 0) THEN
  
  
!     2.1 UPWIND AND LW-FLUXES AT INTERIOR POINTS
!     -------------------------------------------
  
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) THEN
        DO  k=1,nz
          uwt = SIGN(1.0,utemp(k,j,i))
          cfl = utemp(k,j,i)*del3/gx2u(j,i)
          advfup(k,j,i) = 0.5*utemp(k,j,i)  &
              *((1.0+uwt)*phi(k,j,i-1)+(1.0-uwt)*phi(k,j,i))
          advflw(k,j,i) = 0.5*utemp(k,j,i)  &
              *((1.0+cfl)*phi(k,j,i-1)+(1.0-cfl)*phi(k,j,i))
        END DO
      END IF
    END DO
  END DO
  211   CONTINUE
  
  
!     2.2 ADVECTIVE FLUXES AT INTERIOR POINTS
!     ---------------------------------------
  
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) THEN
        DO  k=1,nz
          psi = 0.0
          IF (jadvs == 1) THEN
!          ---fully upwind
            psi = 0.0
          ELSE IF (jadvs == 2) THEN
!          ---Lax-Wendroff
            psi = 1.0
          ELSE IF ((jadvs == 3).OR.(jadvs == 4)) THEN
!          ---TVD scheme
            uadvu = utemp(k,j,i)
            IF (uadvu > 0.0) THEN
              psi = fnlim(advflw(k,j,i),  advfup(k,j,i),  &
                  advflw(k,j,i-1),advfup(k,j,i-1),jadvs)
            ELSE IF (uadvu < 0.0) THEN
              psi = fnlim(advflw(k,j,i),  advfup(k,j,i),  &
                  advflw(k,j,i+1),advfup(k,j,i+1),jadvs)
            END IF
          END IF
          advfl(k,j,i) = (advfup(k,j,i) +  &
              psi*(advflw(k,j,i)-advfup(k,j,i))) *0.5*(gz(k,j,i-1)+gz(k,j,i))
        END DO
      END IF
    END DO
  END DO
  221   CONTINUE
  
  
!     2.3 ADVECTIVE FLUXES AT U-OPEN BOUNDARIES
!     -----------------------------------------
  
  DO  ii=1,nobu
    i = iobu(ii)
    j = jobu(ii)
    ivp = ivpu(ii)
    DO  k=1,nz
      uwt = SIGN(1.0,utemp(k,j,i))
      IF (westob(ii)) THEN
        cfl = utemp(k,j,i)*del3/gx2(j,i)
        simp(k,j,i) = simp(k,j,i) - 0.5*(1.0-uwt)*cfl
        IF (phivp(k,ivp) <= phivp(0,ivp)) THEN
          sexp(k,j,i) = sexp(k,j,i) + 0.5*(1.0+uwt)*cfl*phi(k,j,i)
        ELSE
          sexp(k,j,i) = sexp(k,j,i) + 0.5*(1.0+uwt)*cfl*phivp(k,ivp)
        END IF
      ELSE
        cfl = utemp(k,j,i)*del3/gx2(j,i-1)
        simp(k,j,i-1) = simp(k,j,i-1) + 0.5*(1.0+uwt)*cfl
        IF (phivp(k,ivp) <= phivp(0,ivp)) THEN
          sexp(k,j,i-1) = sexp(k,j,i-1) - 0.5*(1.0-uwt)*cfl *phi(k,j,i-1)
        ELSE
          sexp(k,j,i-1) = sexp(k,j,i-1) - 0.5*(1.0-uwt)*cfl *phivp(k,ivp)
        END IF
      END IF
    END DO
  END DO
  231   CONTINUE
  
END IF


!     3. DIFFUSIVE FLUXES
!     -------------------

IF (jdifs > 0) THEN
  
  
!     3.1 INTERIOR POINTS
!     -------------------
  
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) THEN
        DO  k=1,nz
          difc = htemp(k,j,i)*(phi(k,j,i)-phi(k,j,i-1))/gx2u(j,i)
          diffl(k,j,i) = difc*0.5*(gz(k,j,i-1)+gz(k,j,i))
        END DO
      END IF
    END DO
  END DO
  311   CONTINUE
  
  
!     3.2 U-OPEN BOUNDARIES
!     ---------------------
  
  DO  ii=1,nobu
    i = iobu(ii)
    j = jobu(ii)
    ivp = ivpu(ii)
    DO  k=1,nz
      IF (phivp(k,ivp) > phivp(0,ivp)) THEN
        IF (westob(ii)) THEN
          simp(k,j,i) = simp(k,j,i) + del3*htemp(k,j,i)/(gx2(j,i)**2)
          sexp(k,j,i) = sexp(k,j,i) + del3*htemp(k,j,i)  &
              *phivp(k,ivp)/(gx2(j,i)**2)
        ELSE
          simp(k,j,i-1) = simp(k,j,i-1) + del3*htemp(k,j,i) /(gx2(j,i-1)**2)
          sexp(k,j,i-1) = sexp(k,j,i-1) + del3*htemp(k,j,i)  &
              *phivp(k,ivp)/(gx2(j,i-1)**2)
        END IF
      END IF
    END DO
  END DO
  321   CONTINUE
  
END IF


!     4. ADVECTIVE AND DIFFUSIVE TERMS
!     --------------------------------

!     4.1 ADVECTION
!     -------------
IF (jadvs > 0) THEN
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
   sexp(k,j,i)=sexp(k,j,i)-del3*(advfl(k,j,i+1)-advfl(k,j,i))/(gx2(j,i)*gz(k,j,i))
        END DO
      END IF
    END DO
  END DO
  411   CONTINUE
END IF

!     4.2 DIFFUSION
!     -------------

IF (jdifs > 0) THEN
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
          sexp(k,j,i) = sexp(k,j,i) + del3* (diffl(k,j,i+1)-diffl(k,j,i))  &
              /(gx2(j,i)*gz(k,j,i))
        END DO
      END IF
    END DO
  END DO
  421   CONTINUE
END IF

RETURN

END SUBROUTINE xadvdis

!========================================================================

SUBROUTINE yadvdis (phi,hedphi,vadv,gz,nzphi,ivpv,  &
    phivp,jadvs,jdifs,simp,sexp)
!************************************************************************

!    *YADVDIS*    EVALUATE THE ADVECTIVE AND DIFFUSION TERMS ALONG THE
!                 Y-DIRECTION IN THE TRANSPORT EQUATION FOR A SCALAR
!                 QUANTITY (CELL-CENTERED OR AT W-NODES)

!       AUTHOR - KEVIN RUDDICK AND PATRICK LUYTEN

!       LAST UPDATE - 26 Mar 1999 @(COHERENS)advdis.f 8.4

!       DESCRIPTION - COMPUTES Y-DIRECTION HORIZONTAL ADVECTION/DIFFUSION
!                     OF A 3-D SCALAR QUANTITY HELD AT THE CELL CENTRE
!                     (E.G. CONCENTRATION OF CONTAMINANT, SUSPENDED PARTICLES,
!                     SALINITY) OR AT W-NODES (TURBULENCE VARIABLE)
!                   - INCORPORATES OPEN BOUNDARY CONDITIONS
!                   - ADVECTIVE/DIFFUSIVE TERMS ARE EVALUATED EXPLICITLY
!                   - ADVECTION SCHEME SELECTED BY JADVS
!                      = 1 => FIRST ORDER UPWIND
!                      = 2 => LAX-WENDROFF
!                      = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                      = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!                   - DIFFUSION SELECTED BY JDIFS
!                      = 0 => OFF
!                      = 1 => UNIFORM
!                      = 2 => SMAGORINSKY
!                   - FLUXES ARE CALCULATED AT V- OR Y-NODES

!       REFERENCE - Sections III-4.4.2-4.4.3 of the User Documentation

!       CALLING PROGRAM - TRANSPC, TRANSPW


!       EXTERNALS - FNLIM, ZERO

!************************************************************************
!*    ARGUMENTS

INTEGER, INTENT(IN)                      :: nzphi
REAL, INTENT(IN)                         :: phi(nzphi,nr,nc)
REAL, INTENT(IN)                         :: hedphi(nzphi,nr+1,nc)
REAL, INTENT(IN)                         :: vadv(nzphi,nr+1,nc)
REAL, INTENT(IN)                         :: gz(nzphi,nr,nc)
INTEGER, INTENT(IN)                      :: ivpv(0:nobv)
REAL, INTENT(IN)                         :: phivp(0:nzphi,0:nvprof)
INTEGER, INTENT(IN)                      :: jadvs
INTEGER, INTENT(IN)                      :: jdifs
REAL, INTENT(IN OUT)                        :: simp(nzphi,nr,nc)
REAL, INTENT(IN OUT)                        :: sexp(nzphi,nr,nc)

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *PHI*      REAL      SCALAR TO DE ADVECTED/DIFFUSED [PHI]        (IN)
!    *HEDPHI*   REAL      HORIZONTAL DIFFUSION COEFFICIENT [m2/s]     (IN)

!    *VADV*     REAL      ADVECTING VELOCITY [m/s]                    (IN)
!    *GZ*       REAL      VERTICAL GRID SPACING ABOUT PHI             (IN)
!    *NZPHI*    INT       VERTICAL DIMENSION OF PHI (=NZ OR NZ+1)     (IN)
!    *IVPV*     INT       NO. OF PROFILE AT V-O.B.                    (IN)
!    *PHIVP*    REAL      OPEN BOUNDARY CONDITION [PHI]               (IN)
!                       PHIVP(K,*)>PHIVP(0,*)  : IMPOSED INFLOW
!                       PHIVP(K,*)<= PHIVP(0,*): NORMAL GRADIENT CONDITION
!    *JADVS*    INTEGER   SWITCH TO SELECT TYPE OF ADVECTION SCHEME (0,1,2,3)
!                                                                     (IN)
!    *JDIFS*    INTEGER   SWITCH TO SELECT TYPE OF HORIZONTAL DIFFUSION
!                         (0,1,2)                                     (IN)
!    *SIMP*     REAL      LEFT HAND SIDE OF PHI-EQUATION /PHI(I,J,K)
!                                                                 (IN/OUT)
!    *SEXP*     REAL      RIGHT HAND SIDE OF PHI-EQUATION  [PHI]
!                                                                 (IN/OUT)

!------------------------------------------------------------------------

!*    LOCAL VARIABLES

REAL :: advfup(nz,nr+1,nc), advflw(nz,nr+1,nc)
REAL :: advfl(nz,nr+1,nc), diffl(nz,nr+1,nc)
REAL :: vtemp(nz,nr+1,nc), htemp(nz,nr+1,nc)
INTEGER :: i, ivp, j, jj, k
REAL :: cfl, difc, psi
REAL :: vadvv, vwt
!REAL :: fnlim, gy2v

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ADVFUP*   REAL      UPWIND ADVECTIVE FLUX OF PHI IN Y-DIRECTION AT V-NODE
!                                                                [PHI*m/s]
!    *ADVFLW*   REAL      LW ADVECTIVE FLUX OF PHI IN Y-DIRECTION AT V-NODE
!                                                                [PHI*m/s]
!    *ADVFL*    REAL      TVD ADVECTIVE FLUX OF PHI IN Y-DIRECTION AT V-NODE
!                                                                [PHI*m/s]
!    *DIFFL*    REAL      DIFFUSIVE FLUX OF PHI IN Y-DIRECTION AT V-NODE
!                                                                [PHI*m/s]
!    *CFL*      REAL      COURANT-FRIEDRICH-LEVY NUMBER
!    *PSI*      REAL      FLUX LIMITER
!    *VWT*      REAL      SIGN OF ADVECTIVE VELOCITY

!------------------------------------------------------------------------

!     1. INITIALISATION OF TEMPORARY ARRAYS
!     -------------------------------------


CALL zero33(advfup,nz,nr+1,nc)
CALL zero33(advflw,nz,nr+1,nc)
CALL zero33(advfl,nz,nr+1,nc)
CALL zero33(diffl,nz,nr+1,nc)

DO  i=1,nc
DO  j=1,nr+1
DO  k=1,nz
  vtemp(k,j,i) = vadv(k,j,i)
  htemp(k,j,i) = hedphi(k,j,i)
END DO
END DO
END DO

IF(trigger==1.0)call zero33(vtemp,nz,nr+1,nc)
IF(trigger==2.0)call zero33(htemp,nz,nr+1,nc)

!     2. ADVECTIVE FLUXES AT V-NODES
!     ------------------------------

IF (jadvs > 0) THEN
  
  
!     2.1 UPWIND AND LW-FLUXES AT INTERIOR POINTS
!     -------------------------------------------
  
  
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        DO  k=1,nz
          vwt = SIGN(1.0,vtemp(k,j,i))
          cfl = vtemp(k,j,i)*del3/gy2v(j)
          advfup(k,j,i) = 0.5*vtemp(k,j,i)  &
              *((1.0+vwt)*phi(k,j-1,i)+(1.0-vwt)*phi(k,j,i))
          advflw(k,j,i) = 0.5*vtemp(k,j,i)  &
              *((1.0+cfl)*phi(k,j-1,i)+(1.0-cfl)*phi(k,j,i))
        END DO
      END IF
    END DO
  END DO
  211   CONTINUE
  
  
!     2.2 ADVECTIVE FLUXES AT INTERIOR POINTS
!     ---------------------------------------
  
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        DO  k=1,nz
          psi = 0.0
          IF (jadvs == 1) THEN
!         ---fully upwind
            psi = 0.0
          ELSE IF (jadvs == 2) THEN
!          ---Lax-Wendroff
            psi = 1.0
          ELSE IF ((jadvs == 3).OR.(jadvs == 4)) THEN
!          ---TVD scheme
            vadvv = vtemp(k,j,i)
            IF (vadvv > 0.0) THEN
              psi = fnlim(advflw(k,j,i),  advfup(k,j,i),  &
                  advflw(k,j-1,i),advfup(k,j-1,i),jadvs)
            ELSE IF (vadvv < 0.0) THEN
              psi = fnlim(advflw(k,j,i),  advfup(k,j,i),  &
                  advflw(k,j+1,i),advfup(k,j+1,i),jadvs)
            END IF
          END IF
          advfl(k,j,i) = (advfup(k,j,i) +  &
              psi*(advflw(k,j,i)-advfup(k,j,i)))  &
              *0.5*(gz(k,j-1,i)+gz(k,j,i))*cosphiv(j)
        END DO
      END IF
    END DO
  END DO
  221   CONTINUE
  
  
!     2.3 ADVECTIVE FLUXES AT V-OPEN BOUNDARIES
!     -----------------------------------------
  
  DO  jj=1,nobv
    i = iobv(jj)
    j = jobv(jj)
    ivp = ivpv(jj)
    DO  k=1,nz
      vwt = SIGN(1.0,vtemp(k,j,i))
      IF (soutob(jj)) THEN
        cfl = cosphiv(j)*vtemp(k,j,i)*del3/(gy2(j)*cosphi(j))
        simp(k,j,i) = simp(k,j,i) - 0.5*(1.0-vwt)*cfl
        IF (phivp(k,ivp) <= phivp(0,ivp)) THEN
          sexp(k,j,i) = sexp(k,j,i) + 0.5*(1.0+vwt)*cfl*phi(k,j,i)
        ELSE
          sexp(k,j,i) = sexp(k,j,i) + 0.5*(1.0+vwt)*cfl*phivp(k,ivp)
        END IF
      ELSE
        cfl = cosphiv(j)*vtemp(k,j,i)*del3/(gy2(j-1)*cosphi(j-1))
        simp(k,j-1,i) = simp(k,j-1,i) + 0.5*(1.0+vwt)*cfl
        IF (phivp(k,ivp) <= phivp(0,ivp)) THEN
          sexp(k,j-1,i) = sexp(k,j-1,i) - 0.5*(1.0-vwt)*cfl *phi(k,j-1,i)
        ELSE
          sexp(k,j-1,i) = sexp(k,j-1,i) - 0.5*(1.0-vwt)*cfl *phivp(k,ivp)
        END IF
      END IF
    END DO
  END DO
  231   CONTINUE
  
END IF


!     3. DIFFUSIVE FLUXES
!     -------------------

IF (jdifs > 0) THEN
  
  
!     3.1 INTERIOR POINTS
!     -------------------
  
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        DO  k=1,nz
          difc = htemp(k,j,i)*cosphiv(j) *(phi(k,j,i)-phi(k,j-1,i))/gy2v(j)
          diffl(k,j,i) = difc*0.5*(gz(k,j-1,i)+gz(k,j,i))
        END DO
      END IF
    END DO
  END DO
  311   CONTINUE
  
  
!     3.2 V-OPEN BOUNDARIES
!     ---------------------
  
  DO  jj=1,nobv
    i = iobv(jj)
    j = jobv(jj)
    ivp = ivpv(jj)
    DO  k=1,nz
      IF (phivp(k,ivp) > phivp(0,ivp)) THEN
        IF (soutob(jj)) THEN
          simp(k,j,i) = simp(k,j,i) + del3*htemp(k,j,i)  &
              *cosphiv(j)/(cosphi(j)*gy2(j)**2)
          sexp(k,j,i) = sexp(k,j,i) + del3*htemp(k,j,i)  &
              *cosphiv(j)*phivp(k,ivp)/(cosphi(j)*gy2(j)**2)
        ELSE
          simp(k,j-1,i) = simp(k,j-1,i) + del3*htemp(k,j,i)  &
              *cosphiv(j)/(cosphi(j-1)*gy2(j-1)**2)
          sexp(k,j-1,i) = sexp(k,j-1,i) + del3*htemp(k,j,i)  &
              *cosphiv(j)*phivp(k,ivp)/(cosphi(j-1)*gy2(j-1)**2)
        END IF
      END IF
    END DO
  END DO
  321   CONTINUE
  
END IF


!     4. ADVECTIVE AND DIFFUSIVE TERMS
!     --------------------------------

!     4.1 ADVECTION
!     -------------

IF (jadvs > 0) THEN
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
          sexp(k,j,i) = sexp(k,j,i) - del3* (advfl(k,j+1,i)-advfl(k,j,i))  &
              /(gy2(j)*gz(k,j,i)*cosphi(j))
        END DO
      END IF
    END DO
  END DO
  411   CONTINUE
END IF

!     4.2 DIFFUSION
!     -------------

IF (jdifs > 0) THEN
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
          sexp(k,j,i) = sexp(k,j,i) + del3* (diffl(k,j+1,i)-diffl(k,j,i))  &
              /(gy2(j)*gz(k,j,i)*cosphi(j))
        END DO
      END IF
    END DO
  END DO
  421   CONTINUE
END IF

RETURN

END SUBROUTINE yadvdis

!=======================================================================

SUBROUTINE zadvdis(phi,vedphi,wadv,gz, iblock,isur,ibot,nzphi,nrphi,ncphi,  &
    sflux,stran,sphi,bflux,btran,bphi, jadvw,simpa,simpb,simpc,sexp)
!************************************************************************

!    *ZADVDIS*     EVALUATE Z-DIRECTION ADVECTION-DISPERSION TERMS

!       AUTHOR - KEVIN RUDDICK AND PATRICK LUYTEN

!       LAST UPDATE - 17 Aug 1998 @(COHERENS)advdis.f 8.4

!       DESCRIPTION - COMPUTES Z-DIRECTION ADVECTION/DIFFUSION TERMS
!                     FOR A 3-D QUANTITY
!                   - INCORPORATES SURFACE AND BOTTOM BOUNDARY CONDITIONS
!                   - IS USED FOR QUANTITIES STORED AT EITHER U, V, W OR
!                     CENTRAL NODES BY PASSING THE VERTICAL DIFFUSION
!                     COEFFICIENTS AND VELOCITIES BETWEEN THE UNKNOWNS PHI
!                   - VERTICAL ADVECTION IS EVALUATED SEMI-IMPLICITLY
!                   - VERTICAL DIFFUSION IS EVALUATED FULLY IMPLCITLY
!                   - ADVECTION IS WEIGHTED BETWEEN CENTRAL
!                     (SECOND ORDER SPACE) AND UPWIND (FIRST ORDER, MONOTONIC)
!                   - IF DIRICHLET BOUNDARY CONDITIONS ARE APPLIED
!                     ZADVDIS MUST BE CALLED AFTER ALL OTHER ROUTINES
!                     AFFECTING COEFFICIENTS SIMPA/B/C, SEXP SINCE THESE ARE
!                     OVERWRITTEN
!                   - VARIOUS ADVECTION SCHEMES CAN BE SELECTED :
!                       JADVW = 1 => FIRST ORDER UPWIND
!                             = 2 => SECOND ORDER CENTRAL
!                             = 3 => TVD SCHEME WITH SUPERBEE LIMITING
!                             = 4 => TVD SCHEME WITH MONOTONIC LIMITING

!       REFERENCE - Sections III-4.3.7,4.3.8,4.4.4,4.4.5 of the User
!                    Documentation

!       CALLING PROGRAM - TRANSPC, TRANSPW, UCALC, VCALC

!       EXTERNALS - FNLIM, ZERO

!************************************************************************
!*    ARGUMENTS

INTEGER, INTENT(IN)                      :: nzphi
INTEGER, INTENT(IN)                      :: nrphi
INTEGER, INTENT(IN)                      :: ncphi

REAL, INTENT(IN)                         :: phi(nzphi,nrphi,ncphi)
REAL, INTENT(IN)                         :: vedphi(nz+1,nrphi,ncphi)
REAL, INTENT(IN)                         :: wadv(nz+1,nrphi,ncphi)
REAL, INTENT(IN)                         :: gz(nzphi,nrphi,ncphi)
INTEGER, INTENT(IN)                  :: iblock(nrphi,ncphi)
INTEGER, INTENT(IN)                  :: isur
INTEGER, INTENT(IN)                  :: ibot
REAL, INTENT(IN)                     :: sflux(nrphi,ncphi)
REAL, INTENT(IN)                     :: stran(nrphi,ncphi)
REAL, INTENT(IN)                         :: sphi(nrphi,ncphi)
REAL, INTENT(IN)                     :: bflux(nrphi,ncphi)
REAL, INTENT(IN)                     :: btran(nrphi,ncphi)
REAL, INTENT(IN)                         :: bphi(nrphi,ncphi)
INTEGER, INTENT(IN)                  :: jadvw
REAL, INTENT(IN OUT)                        :: simpa(nzphi,nrphi,ncphi)
REAL, INTENT(IN OUT)                        :: simpb(nzphi,nrphi,ncphi)
REAL, INTENT(IN OUT)                        :: simpc(nzphi,nrphi,ncphi)
REAL, INTENT(IN OUT)                        :: sexp(nzphi,nrphi,ncphi)



!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *PHI*      REAL      SCALAR TO DE ADVECTED/DIFFUSED [PHI]        (IN)
!    *VEDPHI*   REAL      VERT DIFFUSION COEFFICIENT BELOW PHI [m2/s] (IN)
!    *WADV*     REAL      TRANSFORMED VERT VELOCITY BELOW PHI [m/s]   (IN)
!    *GZ*       REAL      VERTICAL CELL SIZE (AROUND PHI) [m] (IN)
!    *IBLOCK*   INT       CELL/FACE FLAG (=0 IF BLOCKED)              (IN)
!    *ISUR*     INTEGER   TYPE OF SURFACE BOUNDARY CONDITION          (IN)
!    *IBOT*     INTEGER   TYPE OF BOTTOM  BOUNDARY CONDITION          (IN)

!     ISUR=0 (NEUMANN)   => SURFACE FLUX = SFLUX + STRAN*( SPHI - PHI)
!     ISUR=1 (DIRICHLET) => SURFACE PHI  = SPHI
!     ISUR=2 (DIRICHLET) => NEAR SURFACE PHI  = SPHI, SURFACE PHI  = 0.0

!     IBOT=0 (NEUMANN)   => BOTTOM  FLUX = BFLUX + BTRAN*( BPHI - PHI)
!     IBOT=1 (DIRICHLET) => BOTTOM  PHI  = SPHI
!     IBOT=2 (DIRICHLET) => NEAR BOTTOM  PHI  = SPHI, BOTTOM  PHI  = 0.0

!    *NCPHI*    INT       X-DIMENSION OF PHI (=NC OR NC+1)            (IN)
!    *NRPHI*    INT       Y-DIMENSION OF PHI (=NR OF NR+1)            (IN)
!    *NZPHI*    INT       VERTICAL DIMENSION OF PHI (=NZ OR NZ+1)     (IN)
!    *SFLUX*    REAL      SEE ABOVE [PHI*m/s]                         (IN)
!    *STRAN*    REAL      SURFACE TRANSFER COEFFICIENT [m/s]          (IN)
!    *SPHI*     REAL      SURFACE VALUE OF PHI [PHI]                  (IN)
!    *BFLUX*    REAL      SEE ABOVE [PHI*m/s]                         (IN)
!    *BTRAN*    REAL      BOTTOM TRANSFER COEFFICIENT [m/s]           (IN)
!    *BPHI*     REAL      BOTTOM VALUE OF PHI [PHI]                   (IN)
!    *JADVW*    INTEGER   ADVECTION SCHEME                            (IN)
!                         (1 => UPWIND ;
!                          2 => CENTRAL ;
!                          3 => TVD WITH SUPERBEE LIMITING FUNCTION ;
!                          4 => TVD WITH MONOTONIC LIMITING FUNCTION)
!    *SIMPA*    REAL      LEFT HAND SIDE OF PHI-EQUATION /PHI(K-1)(IN/OUT)
!    *SIMPB*    REAL      LEFT HAND SIDE OF PHI-EQUATION /PHI(K)  (IN/OUT)
!    *SIMPC*    REAL      LEFT HAND SIDE OF PHI-EQUATION /PHI(K+1)(IN/OUT)
!    *SEXP*     REAL      RIGHT HAND SIDE OF PHI-EQUATION  [PHI]  (IN/OUT)

!------------------------------------------------------------------------


!*    LOCAL VARIABLES

INTEGER :: i, j, k
REAL :: advfup(nz+2), advfce(nz+2)
REAL :: cfl,  dzc, dzl, dzu, fac, psi
REAL :: ved, wadvw, wwt
!REAL :: fnlim
REAL, PARAMETER :: dthimp = 1.0
REAL, PARAMETER :: thimp = 0.501

REAL :: wtemp(nz+1,nr+1,nc), vtemp(nz+1,nr,nc)

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CFL*      REAL      COURANT-FRIEDRICH-LEVY NUMBER
!    *DTHIMP*   REAL      IMPLICIT FACTOR FOR DIFFUSION
!                         (=0 => EXPLICIT (UNSTABLE) ; =+1 => FULLY IMPLICIT
!                          =0.5 => CRANK-NICOLSON)
!    *THIMP*    REAL      IMPLICIT FACTOR FOR ADVECTION
!                         (=0 => EXPLICIT ; =+1 => FULLY IMPLICIT
!                          =0.5 => CRANK-NICOLSON)
!    *ADVFUP*   REAL      UPWIND ADVECTIVE FLUX IN Z-DIRECTION BELOW PHI
!                                                              [PHI*m/s]
!    *ADVFCE*   REAL      CENTRAL ADVECTIVE FLUX IN Z-DIRECTION BELOW PHI
!                                                              [PHI*m/s]
!    *PSI*      REAL      FLUX LIMITER
!    *WWT*      REAL      SIGN OF ADVECTIVE VELOCITY

!------------------------------------------------------------------------


DO  i=1,nc
DO  j=1,nr
DO  k=1,nz+1
  wtemp(k,j,i) = wadv(k,j,i)
  vtemp(k,j,i) = vedphi(k,j,i)
END DO
END DO
END DO

IF(trigger==1.0)call zero33(wtemp,nz+1,nr,nc)
IF(trigger==2.0)call zero33(vtemp,nz+1,nr,nc)

!     1. CALCULATE ADVECTIVE FLUXES
!     -----------------------------

!     1.1 FACE BY FACE
!     ----------------

DO  i=1,ncphi
  DO  j=1,nrphi
    
    IF ((iblock(j,i) == 1).AND.(jadvw > 0)) THEN
      
!       ---first sweep calculate upwind and CE fluxes for a vertical
      CALL zero33(advfup,nz+2,1,1)
      CALL zero33(advfce,nz+2,1,1)
      DO  k=2,nzphi
        wwt = SIGN(1.0,wtemp(k,j,i))
        advfup(k) = 0.5*wtemp(k,j,i)  &
            *((1.0+wwt)*phi(k-1,j,i)+(1.0-wwt)*phi(k,j,i))
        advfce(k) = wtemp(k,j,i)*0.5*(phi(k-1,j,i)+phi(k,j,i))
      END DO
!       ---second sweep use limiters to determine required flux
      DO  k=2,nzphi
        wadvw = wtemp(k,j,i)
        psi = 0.0
        IF (jadvw == 1) THEN
!        --fully upwind
          psi = 0.0
        ELSE IF (jadvw == 2) THEN
!        --central
          psi = 1.0
        ELSE IF ((jadvw == 3).OR.(jadvw == 4)) THEN
!        --TVD scheme
          IF (wadvw > 0.0) THEN
            psi = fnlim(advfce(k),  advfup(k), advfce(k-1),advfup(k-1),jadvw)
          ELSE IF (wadvw < 0.0) THEN
            psi = fnlim(advfce(k),  advfup(k), advfce(k+1),advfup(k+1),jadvw)
          END IF
        END IF
        fac = (1.0-psi)*SIGN(1.0,wadvw)
        
        
!        ---flux below PHI
        
        cfl = wadvw*del3/gz(k,j,i)
        simpa(k,j,i) = simpa(k,j,i) - thimp*cfl*0.5*(1.0+fac)
        simpb(k,j,i) = simpb(k,j,i) - thimp*cfl*0.5*(1.0-fac)  
        sexp(k,j,i)  = sexp(k,j,i) + (1.0-thimp)*cfl*0.5*  &
            (phi(k,j,i)*(1.0-fac)+phi(k-1,j,i)*(1.0+fac))
     
!        ---flux above PHI
        
        cfl = wadvw*del3/gz(k-1,j,i)
        simpb(k-1,j,i) = simpb(k-1,j,i) + thimp*cfl*0.5*(1.0+fac)
        simpc(k-1,j,i) = simpc(k-1,j,i) + thimp*cfl*0.5*(1.0-fac)
        sexp(k-1,j,i)  = sexp(k-1,j,i) - (1.0-thimp)*cfl*0.5*  &
            (phi(k,j,i)*(1.0-fac)+phi(k-1,j,i)*(1.0+fac))
      END DO
      
    END IF
    
  END DO
END DO


!     2. CALCULATE DIFFUSIVE FLUXES
!     -----------------------------

!     2.1 INTERNAL (HORIZONTAL) FACES
!     -------------------------------

DO  i=1,ncphi
  DO  j=1,nrphi
    IF (iblock(j,i) == 1) THEN
      DO  k=2,nzphi
        dzu = gz(k,j,i)
        dzc = 0.5*(gz(k,j,i)+gz(k-1,j,i))
        dzl = gz(k-1,j,i)
        
!        --flux below PHI
        
        ved = vtemp(k,j,i)
        simpa(k,j,i) = simpa(k,j,i) - dthimp*del3*ved/(dzu*dzc)
        simpb(k,j,i) = simpb(k,j,i) + dthimp*del3*ved/(dzu*dzc)
        sexp(k,j,i) = sexp(k,j,i) - (1.0-dthimp)*del3*ved*  &
            (phi(k,j,i)-phi(k-1,j,i))/(dzu*dzc)
!        --flux above PHI
        
        simpb(k-1,j,i) = simpb(k-1,j,i) + dthimp*del3*ved/(dzc*dzl)
        simpc(k-1,j,i) = simpc(k-1,j,i) - dthimp*del3*ved/(dzc*dzl)
        sexp(k-1,j,i) = sexp(k-1,j,i) + (1.0-dthimp)*del3*ved*  &
            (phi(k,j,i)-phi(k-1,j,i))/(dzc*dzl)
      END DO
    END IF
  END DO
END DO


!       2.2 BOTTOM AND SURFACE BOUNDARY FLUXES
!       --------------------------------------

DO  i=1,ncphi
  DO  j=1,nrphi
    IF (iblock(j,i) == 1) THEN
!       ---surface
      IF (isur == 0) THEN
!        --Neumann
        sexp(nzphi,j,i) = sexp(nzphi,j,i) + sflux(j,i)*del3/gz(nzphi,j,i) +  &
            stran(j,i)*sphi(j,i)*del3/gz(nzphi,j,i)
        simpb(nzphi,j,i) = simpb(nzphi,j,i) + del3*stran(j,i)/gz(nzphi,j,i)
      ELSE IF (isur == 1) THEN
!        --Dirichlet
        simpa(nzphi,j,i) = 0.0
        simpb(nzphi,j,i) = 1.0
        simpc(nzphi,j,i) = 0.0
        sexp(nzphi,j,i)  = sphi(j,i)
      ELSE IF (isur == 2) THEN
!        --Dirichlet applied at first point below surface
!          (for kl and eps equation)
        simpa(nzphi,j,i) = 0.0
        simpb(nzphi,j,i) = 1.0
        simpc(nzphi,j,i) = 0.0
        sexp(nzphi,j,i)  = 0.0
        simpa(nzphi-1,j,i) = 0.0
        simpb(nzphi-1,j,i) = 1.0
        simpc(nzphi-1,j,i) = 0.0
        sexp(nzphi-1,j,i)  = sphi(j,i)
      END IF
!       ---bottom
      IF (ibot == 0) THEN
!        --Neumann
        sexp(1,j,i) =  sexp(1,j,i) + bflux(j,i)*del3/gz(1,j,i) +  &
            btran(j,i)*bphi(j,i)*del3/gz(1,j,i)
        simpb(1,j,i) = simpb(1,j,i)  + del3*btran(j,i)/gz(1,j,i)
      ELSE IF (ibot == 1) THEN
!        --Dirichlet
        simpa(1,j,i) = 0.0
        simpb(1,j,i) = 1.0
        simpc(1,j,i) = 0.0
        sexp(1,j,i)  = bphi(j,i)
      ELSE IF (ibot == 2) THEN
!        --Dirichlet applied at first point above bottom
!          (for kl and eps equation)
        simpa(1,j,i) = 0.0
        simpb(1,j,i) = 1.0
        simpc(1,j,i) = 0.0
        sexp(1,j,i)  = 0.0
        simpa(2,j,i) = 0.0
        simpb(2,j,i) = 1.0
        simpc(2,j,i) = 0.0
        sexp(2,j,i)  = bphi(j,i)
      END IF
    END IF
  END DO
END DO


RETURN

END SUBROUTINE zadvdis

!************************************************************************

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:36:10

!    *bcsin.f*    READ OPEN BOUNDARY PARAMETERS AND DATA

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 11 May 1998       @(COHERENS)bcsin.f 8.4

!       DESCRIPTION - READ OPEN BOUNDARY PARAMETER FILE
!                   - CHECK POSSIBLE ERRORS IN OPEN BOUNDARY PARAMETERS
!                   - DETERMINE ORIENTATION OF OPEN BOUNDARY FACES
!                   - READ OPEN BOUNDARY DATA AT SELECTED INTERVALS

!       REFERENCE -

!       CALLING PROGRAM -

!       SUBROUTINES - BBCIN, BCSIN, BC2IN, CBCIN, HBCIN, PBCIN

!************************************************************************


!

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:36:31

SUBROUTINE boundc
!************************************************************************

!    *BOUNDC*     BOUNDARY CONDITIONS FOR THE 2-D MODE

!       AUTHOR - KEVIN RUDDICK

!       LAST UPDATE - 13 Dec 1999        @(COHERENS)boundc.f 8.4

!       DESCRIPTION - USES THE METHOD OF RIEMANN CHARACTERISTICS
!                   - VARIOUS FORMULATIONS (SELECTED BY THE ARRAYS
!                     ITYPOBU AND ITYPOBV) ARE AVAILABLE

!       REFERENCE - Sections III-1.6.3a (general description) and III-4.3.11b
!                   of the User Documentation

!       CALLING PROGRAM - CRRNT2

!       EXTERNALS -

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, icon, ii, j, jj
REAL :: crad, crad0, duu, dvv, dxx, dyy, obah, obcor, obdep
REAL :: obdrag, obdr2x, obdr2y, obdudx, obdvdy, obfs, obgs
REAL :: obpr, obqden, obrhs, obsph
!REAL :: cud2atv, cvd2atu, gx2u, gy2v, h1atc, h1atu, h1atv
!REAL :: ud2atc, vd2atc


!WRITE (*,*) 'TRACE: open sea boundary conditions for UD2/VD2'


!     1. UPDATE PHASES OF TIDAL COMPONENTS
!     ------------------------------------


DO  icon=1,ncon
  phase0(icon) = phase0(icon) + delt*sigma(icon)
  phase0(icon) = MOD(phase0(icon),2.0*pi)
END DO


!     2. COMPUTE UD2 AT EAST/WEST BOUNDARIES
!     --------------------------------------

DO  ii=1,nobu
  i = iobu(ii)
  j = jobu(ii)
  obpr = 0.0
  
  IF (westob(ii)) THEN
    obdep = h1atc(j,i)
    crad = SQRT(g*obdep)
    crad0= SQRT(g*dep(j,i))
    dvv = cosphiv(j+1)*vd2(j+1,i)-cosphiv(j)*vd2(j,i)
    dxx = gx2(j,i)
    dyy = gy2(j)*cosphi(j)
    IF (nwd(j,i+1) == 1) obpr =  &
        -(p2(j,i+1)/r0ref-p2(j,i)/r0ref)*h1atu(j,i+1)/gx2u(j,i+1)
    obfs  = fs(j,i)
    obcor = coriol(j)*cvd2atu(j,i)
    obdvdy= crad*dvv/dyy
    obdrag= fb(j,i+1)
    obdr2x= 2.0*crad*((ud2atc(j,i) - crad*zeta2(j,i)) - r2obu(ii))/dxx
    obqden= udqden(j,i+1)
    obah  = uadhdev(j,i+1) - uah2d(j,i+1) + udh2d(j,i+1)
  ELSE
    obdep = h1atc(j,i-1)
    crad =-SQRT(g*obdep)
    crad0=-SQRT(g*dep(j,i-1))
    dvv = cosphiv(j+1)*vd2(j+1,i-1)-cosphiv(j)*vd2(j,i-1)
    dxx = gx2(j,i-1)
    dyy = gy2(j)*cosphi(j)
    IF (nwd(j,i-2) == 1) obpr = -(p2(j,i-1)/r0ref-p2(j,i-2)/r0ref)  &
        *h1atu(j,i-1)/gx2u(j,i-1)
    obfs  = fs(j,i-1)
    obcor = coriol(j)*cvd2atu(j,i)
    obdvdy= crad*dvv/dyy
    obdrag= fb(j,i-1)
    obdr2x= 2.0*crad*(r2obu(ii) - (ud2atc(j,i-1) - crad*zeta2(j,i-1)))/dxx
    obqden= udqden(j,i-1)
    obah  = uadhdev(j,i-1) - uah2d(j,i-1) + udh2d(j,i-1)
  END IF
  
  r2obu(ii) = r2obu(ii) + delt*(obdr2x + obcor + obpr  &
      + obqden + obfs + obdvdy + obdrag + obah)
  IF (itypobu(ii) == 1) THEN
    r1obu(ii) = 2.0*crad*ampobu(ii,0)
    DO  icon=1,ncon
      r1obu(ii) = r1obu(ii) + 2.0*crad*ampobu(ii,icon)*  &
          COS(phase0(icon)-phaobu(ii,icon))
    END DO
  ELSE IF (itypobu(ii) == 2) THEN
    r1obu(ii) = r1obu(ii) + delt*(obcor + obpr + obqden  &
        + obfs - obdvdy + obdrag + obah)
  ELSE IF (itypobu(ii) == 3) THEN
    obrhs = ampobu(ii,0)
    DO  icon=1,ncon
      obrhs = obrhs + ampobu(ii,icon)* COS(phase0(icon)-phaobu(ii,icon))
    END DO
    r1obu(ii) = r2obu(ii) + 2.0*crad*obrhs
  ELSE IF (itypobu(ii) == 4) THEN
    obrhs = crad0*ampobu(ii,0)
    DO  icon=1,ncon
      obrhs = obrhs + crad0*ampobu(ii,icon)* COS(phase0(icon)-phaobu(ii,icon))
    END DO
    r1obu(ii) = 2.0*obrhs - r2obu(ii)
  END IF
  
END DO

DO  ii=1,nobu
  i = iobu(ii)
  j = jobu(ii)

if(westob(ii))then
 ud2(j,i) = ud2(j,i+1)
else
 ud2(j,i) = ud2(j,i-1)
end if

!  IF (itypobu(ii) > 0) THEN
!    ud2(j,i) = 0.5*(r1obu(ii)+r2obu(ii))
!  ELSE
!    ud2(j,i) = 0.0
!  END IF

END DO


!     3. COMPUTE VD2 AT NORTH/SOUTH BOUNDARIES
!     ----------------------------------------

DO  jj=1,nobv
  i = iobv(jj)
  j = jobv(jj)
  obpr = 0.0
  
  IF (soutob(jj)) THEN
    obdep = h1atc(j,i)
    crad = SQRT(g*obdep)
    crad0= SQRT(g*dep(j,i))
    duu = ud2(j,i+1)-ud2(j,i)
    dxx = gx2(j,i)
    dyy = gy2(j)
    IF (nwd(j+1,i) == 1) obpr =  &
        -(p2(j+1,i)/r0ref-p2(j,i)/r0ref)*h1atv(j+1,i)/gy2v(j+1)
    obgs  = gs(j,i)
    obcor =-coriolv(j)*cud2atv(j,i)
    obdudx= crad*duu/dxx
    obsph = crad*sphcurv(j)*vd2(j,i)
    obdrag= gb(j+1,i)
    obdr2y= 2.0*crad*((vd2atc(j,i) - crad*zeta2(j,i)) - r2obv(jj))/dyy
    obqden= vdqden(j+1,i)
    obah  = vadhdev(j+1,i) - vah2d(j+1,i)+ vdh2d(j+1,i)
  ELSE
    obdep = h1atc(j-1,i)
    crad =-SQRT(g*obdep)
    crad0=-SQRT(g*dep(j-1,i))
    duu = ud2(j-1,i+1)-ud2(j-1,i)
    dxx = gx2(j-1,i)
    dyy = gy2(j-1)
    IF (nwd(j-2,i) == 1) obpr =  &
        -(p2(j-1,i)/r0ref-p2(j-2,i)/r0ref)*h1atv(j-1,i)/gy2v(j-1)
    obgs  = gs(j-1,i)
    obcor =-coriolv(j)*cud2atv(j,i)
    obdudx= crad*duu/dxx
    obsph = crad*sphcurv(j)*vd2(j,i)
    obdrag= gb(j-1,i)
    obdr2y= 2.0*crad*(r2obv(jj) - (vd2atc(j-1,i) - crad*zeta2(j-1,i)))/dyy
    obqden= vdqden(j-1,i)
    obah  = vadhdev(j-1,i) - vah2d(j-1,i) + vdh2d(j-1,i)
  END IF
  
  r2obv(jj) = r2obv(jj) + delt*(obdr2y + obcor + obpr  &
      + obqden + obgs + obdudx - obsph + obdrag + obah)
  
  IF (itypobv(jj) == 1) THEN
    r1obv(jj) = 2.0*crad*ampobv(jj,0)
    DO  icon=1,ncon
      r1obv(jj) = r1obv(jj) + 2.0*crad*ampobv(jj,icon)*  &
          COS(phase0(icon)-phaobv(jj,icon))
    END DO
  ELSE IF (itypobv(jj) == 2) THEN
    r1obv(jj) = r1obv(jj) + delt*(obcor + obpr + obqden  &
        + obgs - obdudx + obsph + obdrag + obah)
  ELSE IF (itypobv(jj) == 3) THEN
    obrhs = ampobv(jj,0)
    DO  icon=1,ncon
      obrhs = obrhs + ampobv(jj,icon)* COS(phase0(icon)-phaobv(jj,icon))
    END DO
    r1obv(jj) = r2obv(jj) + 2.0*crad*obrhs
  ELSE IF (itypobv(jj) == 4) THEN
    obrhs = crad0*ampobv(jj,0)
    DO  icon=1,ncon
      obrhs = obrhs + crad0*ampobv(jj,icon)* COS(phase0(icon)-phaobv(jj,icon))
    END DO
    r1obv(jj) = 2.0*obrhs - r2obv(jj)
  END IF
  
END DO

DO  jj=1,nobv
  i = iobv(jj)
  j = jobv(jj)

if(soutob(jj))then
 vd2(j,i) = vd2(j+1,i)
else
 vd2(j,i) = vd2(j-1,i)
end if

!  IF (itypobv(jj) > 0) THEN
!    vd2(j,i) = 0.5*(r1obv(jj)+r2obv(jj))
!  ELSE
!    vd2(j,i) = 0.0
!  END IF

END DO


RETURN


END SUBROUTINE boundc

!

!

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:37:24

SUBROUTINE contny
!************************************************************************

!    *CONTNY*      SOLVE THE 2-D CONTINUITY EQUATION FOR THE SURFACE
!                  ELEVATION

!       AUTHOR - KEVIN RUDDICK AND PATRICK LUYTEN

!       LAST UPDATE - 26 Mar 1999 @(COHERENS)contny.f 8.4

!       DESCRIPTION - CHECKS THE CFL CONDITION ON FIRST CALL AND STOPS
!                     THE PROGRAM IF THE CONDITION IS NOT SATISFIED
!                   - SOLVE 2-D CONTINUITY EQUATION FOR ZETA2

!       REFERENCE - Sections III-1.1, III-4.3 of the User Documentation

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS - ERROR

!************************************************************************
!*    LOCAL VARIABLES

LOGICAL :: call1
INTEGER :: i, j
SAVE    call1
DATA    call1 /.true./
!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL

!------------------------------------------------------------------------

!     1. CHECK CFL-LIMIT ON FIRST CALL
!     --------------------------------

IF (call1) THEN
  IF (delt > dtmax) THEN
    nerrs = 1
    WRITE (0,'(A)') 'Time step for 2-D mode DELT exceeds '// 'CFL-limit'
    WRITE (0,8001) delt, dtmax
  END IF
END IF

!WRITE (*,*) 'TRACE: ZETA2'


!     2. COMPUTE NEW ELEVATION FROM DEPTH INTEGRATED FLUXES
!     -----------------------------------------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      zeta2(j,i) = zeta2(j,i) - delt*((ud2(j,i+1)-ud2(j,i))/gx2(j,i) +  &
          (vd2(j+1,i)*cosphiv(j+1)-vd2(j,i)*cosphiv(j)) /(gy2(j)*cosphi(j)))
    END IF
  END DO
END DO

RETURN

8001 FORMAT(2X,'2-D time step : ',1PE14.7,/, 2X,'CLF-limit     : ',1PE14.7)
8002 FORMAT ('Invalid value for element (',i3,',',i3,') of array ',  &
    a,' : ',1PE14.7)

END SUBROUTINE contny

!************************************************************************

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:37:42

!    *crrnt2.f*      ENSEMBLE OF ROUTINES TO SOLVE THE 2-D MOMENTUM
!                    EQUATIONS

!       AUTHOR - KEVIN RUDDICK AND PATRICK LUYTEN

!       LAST UPDATE - 13 Dec 1999 @(COHERENS)crrnt2.f 8.4

!       DESCRIPTION -

!       REFERENCE - Section III-4.3 of the User Documentation

!       SUBROUTINES - CRRNT2, UDCALC, VDCALC

!*************************************************************************

!=========================================================================

SUBROUTINE crrnt2
!************************************************************************

!    *CRRNT2*      SOLVE THE 2-D MOMENTUM EQUATIONS (MAIN UNIT)

!       AUTHOR - KEVIN RUDDICK AND PATRICK LUYTEN

!       LAST UPDATE - 26 Jan 1996 @(COHERENS)crrnt2.f 8.0

!       DESCRIPTION - EVALUATE THE 2-D ADVECTIVE AND DIFFUSION TERMS FOR
!                     THE 2-D CURRENTS
!                   - SOLVE THE 2-D MOMENTUM EQUATIONS
!                   - UPDATE THE 2-D CURRENT AT THE OPEN BOUNDARIES
!                   - UPDATE THE "FILTERED" VELOCITIES UD2F, VD2F

!       REFERENCE - Section III-4.3 of the User Documentation

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS - BOUNDC, HAD2DU, HAD2DV, UDCALC, VDCALC

!************************************************************************
!*    LOCAL VARIABLES

LOGICAL :: udfirst
REAL :: eee
INTEGER :: i, j
SAVE udfirst
DATA udfirst /.true./

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *UDFIRST*   LOGICAL   FLAG TO DETERMINE WHETHER UD2 IS CALCULATED
!                          BEFORE VD2 (.TRUE. AT ODD 2-D TIME STEPS)

!------------------------------------------------------------------------

!     1. UPDATE ADVECTIVE TERMS
!     -------------------------

CALL had2du
CALL had2dv


!     2. ALTERNATE BETWEEN COMPUTING UD OR VD FIRST
!     ---------------------------------------------

IF (udfirst) THEN
  CALL udcalc
  CALL vdcalc
ELSE
  CALL vdcalc
  CALL udcalc
END IF
udfirst = .NOT.udfirst

!     3. APPLY OPEN SEA BOUNDARY CONDITIONS
!     -------------------------------------

CALL boundc

!     4. UPDATE 2D MODE AVERAGE (OVER 3D TIME STEP)
!     ---------------------------------------------


DO  i=1,nc+1
  DO  j=1,nr
    ud2f(j,i) = ud2f(j,i) + ud2(j,i)*delt/del3
  END DO
END DO

DO  i=1,nc
  DO  j=1,nr+1
    vd2f(j,i) = vd2f(j,i) + vd2(j,i)*delt/del3
  END DO
END DO


RETURN

END SUBROUTINE crrnt2

!=======================================================================

SUBROUTINE udcalc
!************************************************************************

!    *UDCALC*      SOLVE THE 2-D U-MOMENTUM EQUATION

!       AUTHOR - KEVIN RUDDICK

!       LAST UPDATE - 13 Dec 1999 @(COHERENS)crrnt2.f 8.4

!       DESCRIPTION -

!       REFERENCE - Section III-4.3 of the User Documentation

!       CALLING PROGRAM - CRRNT2

!       EXTERNALS -

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, j
REAL :: ulhs(nr,nc+1), urhs(nr,nc+1)
REAL :: udcor, udfs, udqse, udqsp,umean, adjust
!REAL :: cvd2atu, gx2u, h1atu

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ULHS*     REAL      LEFT  HAND SIDE OF MOMENTUM EQUATION/UD2
!    *URHS*     REAL      RIGHT HAND SIDE OF MOMENTUM EQUATION

!------------------------------------------------------------------------

!     1. INITIALISE TEMPORARY ARRAYS
!     ------------------------------

!WRITE (*,*) 'TRACE: UD2 - internal points'

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      urhs(j,i) = ud2(j,i)
      ulhs(j,i) = 1.0
    END IF
  END DO
END DO

!     2. ADD CORIOLIS, WIND AND BOTTOM STRESS, SURFACE ELEVATION,
!        AND ATMOSPHERIC PRESSURE GRADIENT TERMS
!     -----------------------------------------------------------

umean = 0.0
DO  j=1,nr-1
 umean=umean+ud2(j,2)/h2atu(j,2)
END DO
 umean = umean/real(nr-1)
 
adjust = (umean+0.2)/0.2
adjust = max(adjust,-1.0)
adjust = min(adjust,1.0)
adjust = sign(1.0,adjust)*sqrt(abs(adjust))

WRITE(6,*)"umean,adjust",umean,adjust

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      udcor =  coriol(j)*cvd2atu(j,i)
      udqse = -g*h1atu(j,i)*(zeta2(j,i)-zeta2(j,i-1))/gx2u(j,i)
! forcing JK

        if(i<10) udqse = udqse -g*h1atu(j,i)*etam/gx2u(j,i)*adjust

      udqsp = -(p2(j,i)/r0ref-p2(j,i-1)/r0ref)*h1atu(j,i) /gx2u(j,i)
      udfs  = 0.5*(fs(j,i) + fs(j,i-1))
      urhs(j,i) = urhs(j,i) + delt*(udcor + udqse + udqsp + udfs)
!     quasi-implicit bottom stress
!     normal condition is that FB  is opposite sign to U2 and UD2 (and UDP)
!     and FBK is positive
      urhs(j,i) = urhs(j,i) + delt*(fb(j,i)+fbk(j,i)*udp(j,i)/h1atu(j,i))
      ulhs(j,i) = ulhs(j,i) + delt*fbk(j,i)/h1atu(j,i)
    END IF
  END DO
END DO


!     3. ADD HORIZONTAL ADVECTION AND DIFFUSION TERMS
!     -----------------------------------------------

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      urhs(j,i) = urhs(j,i) + delt* (uadhdev(j,i)-uah2d(j,i)+udh2d(j,i))
    END IF
  END DO
END DO

!     4. ADD DEPTH-INTEGRATED BAROCLINIC PRESSURE GRADIENT
!     ----------------------------------------------------

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      urhs(j,i) = urhs(j,i) + delt*udqden(j,i)
    END IF
  END DO
END DO

!     5. UPDATE MASS TRANSPORT
!  ------------------------

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      ud2(j,i) = urhs(j,i)/ulhs(j,i)
    END IF
  END DO
END DO


RETURN

END SUBROUTINE udcalc
!=======================================================================

SUBROUTINE vdcalc
!************************************************************************

!    *VDCALC*      SOLVE THE 2-D V-MOMENTUM EQUATION

!       AUTHOR - KEVIN RUDDICK

!       LAST UPDATE - 13 Dec 1999 @(COHERENS)crrnt2.f 8.4

!       DESCRIPTION -

!       REFERENCE - Section III-4.3 of the User Documentation

!       CALLING PROGRAM - CRRNT2

!       EXTERNALS -

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, j
REAL :: vlhs(nr+1,nc), vrhs(nr+1,nc)
REAL :: vdcor, vdgs, vdqse, vdqsp
!REAL :: cud2atv, gy2v, h1atv

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *VLHS*     REAL      LEFT  HAND SIDE OF VD MOMENTUM EQUATION/VD2
!    *VRHS*     REAL      RIGHT HAND SIDE OF VD MOMENTUM EQUATION

!------------------------------------------------------------------------

!     1. INITIALISE TEMPORARY ARRAYS
!     ------------------------------

!WRITE (*,*) 'TRACE: VD2 - internal points'

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      vrhs(j,i) = vd2(j,i)
      vlhs(j,i) = 1.0
    END IF
  END DO
END DO

!     2. ADD CORIOLIS, WIND AND BOTTOM STRESS, SURFACE ELEVATION,
!        AND ATMOSPHERIC PRESSURE GRADIENT TERMS
!     ----------------------------------------------------------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN

!      vdcor = -coriolv(j)*(cud2atv(j,i)-ugeo*h1atv(j,i))
      vdcor = -coriolv(j)*cud2atv(j,i)
      vdqse = -g*h1atv(j,i)*(zeta2(j,i)-zeta2(j-1,i))/gy2v(j)
      vdqsp = -(p2(j,i)/r0ref-p2(j-1,i)/r0ref)*h1atv(j,i)/gy2v(j)

! forcing JK
!      IF((h1atv(j,i)>500.0).AND.(h1atv(j,i)<1200.0)) vdqse = vdqse - 1.e-4*0.1*h1atv(j,i)
      vdgs  = 0.5*(gs(j,i) + gs(j-1,i))
      vrhs(j,i) = vrhs(j,i) + delt*(vdcor + vdqse + vdqsp + vdgs)
!     quasi-implicit bottom stress
!     normal condition is that GB  is opposite sign to V2 and VD2 (and VDP)
!     and GBK is positive
      vrhs(j,i) = vrhs(j,i) + delt*(gb(j,i)+gbk(j,i)*vdp(j,i)/h1atv(j,i))
      vlhs(j,i) = vlhs(j,i) + delt*gbk(j,i)/h1atv(j,i)
    END IF
  END DO
END DO

!     3. ADD HORIZONTAL ADVECTION AND DIFFUSION TERMS
!     -----------------------------------------------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      vrhs(j,i) = vrhs(j,i) + delt* (vadhdev(j,i)-vah2d(j,i)+vdh2d(j,i))
    END IF
  END DO
END DO


!     4. ADD DEPTH-INTEGRATED BAROCLINIC PRESSURE GRADIENT
!     ----------------------------------------------------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      vrhs(j,i) = vrhs(j,i) + delt*vdqden(j,i)
    END IF
  END DO
END DO

!     5. UPDATE MASS TRANSPORT
!     ------------------------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      vd2(j,i) = vrhs(j,i)/vlhs(j,i)
    END IF
  END DO
END DO


RETURN

END SUBROUTINE vdcalc

!************************************************************************

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:37:46

!    *crrnt3.f*      ENSEMBLE OF ROUTINES TO SOLVE THE 3-D MOMENTUM
!                    EQUATIONS

!       AUTHOR - KEVIN RUDDICK AND PATRICK LUYTEN

!       LAST UPDATE - 13 Dec 1999 @(COHERENS)crrnt3.f 8.4

!       DESCRIPTION -

!       REFERENCE - Section III-4.3 of the User Documentation

!       SUBROUTINES - CRRNT3P, CRRNT3C, CRRNT1C, UCALC, VCALC,
!                     WCALC, WPCALC

!*************************************************************************

!=======================================================================

SUBROUTINE crrnt3p
!************************************************************************

!    *CRRNT3P*      UPDATE 3-D HORIZONTAL CURRENT (PREDICTOR STEP)

!       AUTHOR - KEVIN RUDDICK AND PATRICK LUYTEN

!       LAST UPDATE - 27 Mar 1998 @(COHERENS)crrnt3.f 8.4

!       DESCRIPTION - STORE OLD VALUES (SURFACE ELEVATION, ...)
!                   - UPDATE U2 AND V2 AT PREDICTOR STEP BY CALLING UCALC,
!                     VCALC
!                   - EVALUATE DEPTH-INTEGRATED CURRENT AT PREDICTOR STEP
!                   - EVALUATE BOTTOM STRESS AT PREDICTOR STEP

!       REFERENCE - Section III-4.3 of the User Documentation

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS - BSTRES, UCALC, VCALC

!************************************************************************
!*    LOCAL VARIABLES

LOGICAL :: ufirst
INTEGER :: i, j, k
!REAL :: gz1u, gz1v
SAVE ufirst
DATA ufirst /.true./


!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *UFIRST*   LOGICAL   FLAG TO DETERMINE WHETHER U2 IS CALCULATED
!                         BEFORE V2 (.TRUE. AT ODD 3-D TIME STEPS)

!------------------------------------------------------------------------

!     1. STORE OLD VALUES
!     -------------------
!     1.1 SURFACE ELEVATION
!     ---------------------

DO  i=1,nc
  DO  j=1,nr
    zeta1(j,i) = zeta2(j,i)
  END DO
END DO



!     1.2 VERTICAL GRID SPACING
!     -------------------------

DO  i=1,nc
  DO  j=1,nr
    DO  k=1,nz
      gz1(k,j,i) = gz2(k,j,i)
    END DO
  END DO
END DO

!     1.3 CURRENTS
!     ------------

!     ---U-nodes
DO  i=1,nc+1
  DO  j=1,nr
    DO  k=1,nz
      u1(k,j,i) = u2(k,j,i)
    END DO
  END DO
END DO

!     ---V-nodes
DO  i=1,nc
  DO  j=1,nr+1
    DO  k=1,nz
      v1(k,j,i) = v2(k,j,i)
    END DO
  END DO
END DO

!     2. INITIALISE 2-D CURRENT AVERAGES AND TEMPORARY ARRAYS
!     -------------------------------------------------------

CALL zero33(ud2f,1,nr,nc+1)
CALL zero33(vd2f,1,nr+1,nc)

!     3. ALTERNATE BETWEEN COMPUTING U OR V FIRST
!     -------------------------------------------

IF (ufirst) THEN
  CALL ucalc
  CALL vcalc
ELSE
  CALL vcalc
  CALL ucalc
END IF
ufirst = .NOT.ufirst

!     4. 2-D TRANSPORT (PREDICTOR STEP)
!     ---------------------------------

DO  i=2,nc
  DO  j=1,nr
    udp(j,i) = 0.0
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        udp(j,i) = udp(j,i) + u2(k,j,i)*gz1u(k,j,i)
      END DO
    END IF
  END DO
END DO

DO  i=1,nc
  DO  j=2,nr
    vdp(j,i) = 0.0
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        vdp(j,i) = vdp(j,i) + v2(k,j,i)*gz1v(k,j,i)
      END DO
    END IF
  END DO
END DO

!     5. UPDATE BOTTOM STRESS
!     -----------------------

IF (igrdim == 3) CALL bstres


RETURN

END SUBROUTINE crrnt3p

!========================================================================

SUBROUTINE crrnt3c
!************************************************************************

!    *CRRNT3C*      UPDATE THE HORIZONTAL CURRENT AT THE CORRECTOR STEP
!                   (3-D APPLICATIONS)

!       AUTHOR - KEVIN RUDDICK AND PATRICK LUYTEN

!       LAST UPDATE - 26 Mar 1999 @(COHERENS)crrnt3.f 8.4

!       DESCRIPTION - UPDATE U2,V2 AND U2F,V2F AT INTERIOR POINTS
!                     USING THE DISCRETISED FORMS OF Section III-4.3.1c
!                   - UPDATE U2,V2 AND U2F,V2F AT OPEN BOUNDARIES USING
!                     THE PROCEDURES DESCRIBED IN Section III-4.3.10c
!                   - UPDATE THE BOTTOM STRESS AT THE NEW (CORRECTED) TIME

!       REFERENCE - Section III-4.3 of the User Documentation
!                 - E.Deleersnijder, 1993 : "Numerical mass conservation
!                   in a free-surface sigma coordinate marine model with
!                   mode splitting", J.Mar.Sys. Vol 4 no 5, pp 365-170.

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS - BSTRES

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, ii, iucorr, ivcorr, ivp, j, jj, jucorr, jvcorr, k
REAL :: ucorr, udint, uerr1, uerr2, vcorr, vdint, verr1, verr2
!REAL :: gz1u, gz2u, gz1v, gz2v, h1atu, h1atv, h2atc

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *IUCORR*   INTEGER   X-INDEX OF THE (INTERIOR) U-NODE WHERE THE
!                         DEVIATION OF THE DEPTH-MEAN CURRENT (UD2/H)
!                         FROM ITS PREDICTED VALUE (UDP/H) HAS ITS
!                         LARGEST MAGNITUDE
!    *JVCORR*   INTEGER   Y-INDEX OF THE (INTERIOR) V-NODE WHERE THE
!                         DEVIATION OF THE DEPTH-MEAN CURRENT (VD2/H)
!                         FROM ITS PREDICTED VALUE (VDP/H) HAS ITS
!                         LARGEST MAGNITUDE
!    *UCORR*    REAL      MAXIMUM DEVIATION (IN MAGNITUDE) OF THE DEPTH-MEAN
!                         U-CURRENT FROM ITS PREDICTED VALUE    [m/s]
!    *VCORR*    REAL      MAXIMUM DEVIATION (IN MAGNITUDE) OF THE DEPTH-MEAN
!                         V-CURRENT FROM ITS PREDICTED VALUE    [m/s]

!------------------------------------------------------------------------


!WRITE (*,*) 'TRACE: correcting u and v - internal points'


!     1. CORRECT U AND V AT INTERNAL POINTS FOR EXTERNAL MODE TRANSPORT
!     -----------------------------------------------------------------

!     1.1 U-NODES
!     -----------


ucorr = 0.0
DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      uerr1 = (ud2(j,i)  - udp(j,i))/h1atu(j,i)
      uerr2 = (ud2f(j,i) - udp(j,i))/h1atu(j,i)
      DO  k=1,nz
!        ---NB careful to correct U2F before U2 to avoid overwriting U2
        u2f(k,j,i) = (u2(k,j,i) + uerr2)*gz1u(k,j,i)/gz2u(k,j,i)
        u2(k,j,i)  = (u2(k,j,i) + uerr1)*gz1u(k,j,i)/gz2u(k,j,i)
      END DO
      IF (ABS(uerr1) > ucorr) THEN
        ucorr = ABS(uerr1)
        iucorr= i
        jucorr= j
      END IF
    END IF
  END DO
END DO
111   CONTINUE

!     1.2 V-NODES
!     -----------


vcorr = 0.0
DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      verr1 = (vd2(j,i)  - vdp(j,i))/h1atv(j,i)
      verr2 = (vd2f(j,i) - vdp(j,i))/h1atv(j,i)
      DO  k=1,nz
!        ---NB careful to correct V2F before V2 to avoid overwriting V2
        v2f(k,j,i) = (v2(k,j,i) + verr2)*gz1v(k,j,i)/gz2v(k,j,i)
        v2(k,j,i)  = (v2(k,j,i) + verr1)*gz1v(k,j,i)/gz2v(k,j,i)
      END DO
      IF (ABS(verr1) > vcorr) THEN
        vcorr = ABS(verr1)
        ivcorr= i
        jvcorr= j
      END IF
    END IF
  END DO
END DO
121   CONTINUE

!     ----output corrections
!WRITE (*,9101) 'U',ucorr,' m/s',iucorr,jucorr
!WRITE (*,9101) 'V',vcorr,' m/s',ivcorr,jvcorr


!     2. APPLY VELOCITY DEVIATION OPEN BOUNDARY CONDITIONS
!     ---------------------------------------------------------

!     2.1 U-OPEN BOUNDARIES
!     ---------------------


!WRITE (*,*) 'TRACE: velocity deviation open boundary conditions'

DO  ii=1,nobu
  DO  k=1,nz
    i = iobu(ii)
    j = jobu(ii)
!    ivp = ivpobu(ii)
!    IF (uvp(k,ivp) <= uvp(0,ivp)) THEN
!        --normal zero gradient boundary condition
      IF (westob(ii)) THEN
        u2(k,j,i)  = u2(k,j,i+1)*gz2u(k,j,i+1)/gz2(k,j,i)
      ELSE
        u2(k,j,i)  = u2(k,j,i-1)*gz2u(k,j,i-1)/gz2(k,j,i-1)
      END IF
!    ELSE
!        --imposed profile of velocity deviation
!      u2(k,j,i)  = uvp(k,ivp)
!    END IF
  END DO
END DO

!     2.2 V-OPEN BOUNDARIES
!     ---------------------


DO  jj=1,nobv
  DO  k=1,nz
    i = iobv(jj)
    j = jobv(jj)
!    ivp = ivpobv(jj)
!    IF (uvp(k,ivp) <= uvp(0,ivp)) THEN
!        --normal zero gradient boundary condition
      IF (soutob(jj)) THEN
        v2(k,j,i)  = v2(k,j+1,i)*gz2v(k,j+1,i)*cosphiv(j+1)  &
            /(gz2(k,j,i)*cosphi(j))
      ELSE
        v2(k,j,i)  = v2(k,j-1,i)*gz2v(k,j-1,i)*cosphiv(j-1)  &
            /(gz2(k,j-1,i)*cosphiv(j))
      END IF
!    ELSE
!        --imposed profile of velocity deviation
!      v2(k,j,i)  = uvp(k,ivp)
!    END IF
  END DO
END DO


!     3. CORRECT U AND V AT OPEN BOUNDARIES
!     -------------------------------------


!     3.1 U-OPEN BOUNDARIES
!     ---------------------


DO  ii=1,nobu
  i = iobu(ii)
  j = jobu(ii)
  udint = 0.0
  IF (westob(ii)) THEN
    DO  k=1,nz
      udint = udint + u2(k,j,i)*gz2(k,j,i)
    END DO
    DO  k=1,nz
!        ---NB careful to correct U2F before U2 to avoid overwriting U2
      u2f(k,j,i) = u2(k,j,i) + (ud2f(j,i)- udint)/h2atc(j,i)
      u2(k,j,i)  = u2(k,j,i) + (ud2(j,i) - udint)/h2atc(j,i)
    END DO
  ELSE
    DO  k=1,nz
      udint = udint + u2(k,j,i)*gz2(k,j,i-1)
    END DO
    DO  k=1,nz
!        ---NB careful to correct U2F before U2 to avoid overwriting U2
      u2f(k,j,i) = u2(k,j,i) + (ud2f(j,i)- udint)/h2atc(j,i-1)
      u2(k,j,i)  = u2(k,j,i) + (ud2(j,i) - udint)/h2atc(j,i-1)
    END DO
  END IF
END DO

!     3.2 V-OPEN BOUNDARIES
!     ---------------------


DO  jj=1,nobv
  i = iobv(jj)
  j = jobv(jj)
  vdint = 0.0
  IF (soutob(jj)) THEN
    DO  k=1,nz
      vdint = vdint + v2(k,j,i)*gz2(k,j,i)
    END DO
    DO  k=1,nz
!        ---NB careful to correct V2F before V2 to avoid overwriting V2
      v2f(k,j,i) = v2(k,j,i) + (vd2f(j,i)- vdint)/h2atc(j,i)
      v2(k,j,i)  = v2(k,j,i) + (vd2(j,i) - vdint)/h2atc(j,i)
    END DO
  ELSE
    DO  k=1,nz
      vdint = vdint + v2(k,j,i)*gz2(k,j-1,i)
    END DO
    DO  k=1,nz
!        --NB careful to correct V2F before V2 to avoid overwriting V2
      v2f(k,j,i) = v2(k,j,i) + (vd2f(j,i)- vdint)/h2atc(j-1,i)
      v2(k,j,i)  = v2(k,j,i) + (vd2(j,i) - vdint)/h2atc(j-1,i)
    END DO
  END IF
END DO

!     4. UPDATE BOTTOM STRESS
!     -----------------------

CALL bstres


RETURN

9101 FORMAT (a1,'ERR=',1PE14.7,a5,' (',i3,',',i3,')')

END SUBROUTINE crrnt3c

!========================================================================

SUBROUTINE crrnt1c
!************************************************************************

!    *CRRNT1C*      UPDATE THE HORIZONTAL CURRENT AT THE CORRECTOR STEP
!                   (1-D APPLICATION)

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 5 May 1999 @(COHERENS)crrnt3.f 8.4

!       DESCRIPTION - UPDATE HORIZONTAL CURRENT AT ALL INTERIOR POINTS
!                   - EVALUATE DEPTH-INTEGRATED CURRENT
!                   - UPDATE CURRENTS AT OPEN BOUNDARIES
!                     (EQUAL TO THE INTERIOR VALUES)
!                   - UPDATE THE BOTTOM STRESS

!       REFERENCE - Sections III-1.1.4 and III-4.7 of the User Documentation

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS - BSTRES

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, ii, j, jj, k
!REAL :: gz1u, gz1v, gz2u, gz2v


!     1. CORRECTOR STEP
!     -----------------

!     ---U-nodes
DO  i=2,nc
  DO  j=1,nr
    DO  k=1,nz
      u2(k,j,i) = u2(k,j,i)*gz1u(k,j,i)/gz2u(k,j,i)
    END DO
  END DO
END DO
!     ---V-nodes
DO  i=1,nc
  DO  j=2,nr
    DO  k=1,nz
      v2(k,j,i) = v2(k,j,i)*gz1v(k,j,i)/gz2v(k,j,i)
    END DO
  END DO
END DO


!     2. DEPTH INTEGRATED CURRENTS
!     ----------------------------

!     ---U-nodes
DO  i=2,nc
  DO  j=1,nr
    ud2(j,i) = 0.0
    DO  k=1,nz
      ud2(j,i) = ud2(j,i)+gz2u(k,j,i)*u2(k,j,i)
    END DO
  END DO
END DO

!     ---V-nodes
DO  i=1,nc
  DO  j=2,nr
    vd2(j,i) = 0.0
    DO  k=1,nz
      vd2(j,i) = vd2(j,i)+gz2v(k,j,i)*v2(k,j,i)
    END DO
  END DO
END DO


!     3. U-OPEN BOUNDARIES
!     ---------------------

DO  ii=1,nobu
  i = iobu(ii)
  j = jobu(ii)
  IF (westob(ii)) THEN
    DO  k=1,nz
      u2(k,j,i) = u2(k,j,i+1)
    END DO
    ud2(j,i) = ud2(j,i+1)
  ELSE
    DO  k=1,nz
      u2(k,j,i) = u2(k,j,i-1)
    END DO
    ud2(j,i) = ud2(j,i-1)
  END IF
END DO


!     4. V-OPEN BOUNDARIES
!     --------------------

DO  jj=1,nobv
  i = iobv(jj)
  j = jobv(jj)
  IF (soutob(jj)) THEN
    DO  k=1,nz
      v2(k,j,i) = v2(k,j+1,i)
    END DO
    vd2(j,i) = vd2(j+1,i)
  ELSE
    DO  k=1,nz
      v2(k,j,i) = v2(k,j-1,i)
    END DO
    vd2(j,i) = vd2(j-1,i)
  END IF
END DO


!     5. UPDATE BOTTOM STRESS
!     -----------------------

CALL bstres


RETURN

END SUBROUTINE crrnt1c

!=======================================================================

SUBROUTINE ucalc
!************************************************************************

!    *UCALC*      SOLVE U-COMPONENT OF THE 3-D MOMENTUM EQUATION
!                 (PREDICTOR STEP)

!       AUTHOR - PATRICK LUYTEN AND KEVIN RUDDICK

!       LAST UPDATE - 13 Dec 1999 @(COHERENS)crrnt3.f 8.4

!       DESCRIPTION - SOLVES U2-EQUATION AT U-NODES
!                   - OPERATOR SPLITTING IS DISABLED IN THE ABSENCE OF
!                     ADVECTION OR FOR THE UPWIND SCHEME, AND IS ENABLED
!                     OTHERWISE

!       REFERENCE - Section III-4.3.1, 4.3.2 of the User Documentation

!       CALLING PROGRAM - CRRNT3P

!       EXTERNALS - THOMV, XHAD3DU, YHAD3DU, ZADVDIS, ZERO

!************************************************************************
!*    LOCAL VARIABLES

REAL :: fsu(nr,nc+1), gzu(nz,nr,nc+1)
REAL :: ulhsa(nz,nr,nc+1), ulhsb(nz,nr,nc+1)
REAL :: ulhsc(nz,nr,nc+1), urhs(nz,nr,nc+1)
REAL :: u2a(nz,nr,nc+1), u2b(nz,nr,nc+1)
REAL :: vedu(nz+1,nr,nc+1), w2u(nz+1,nr,nc+1)
REAL :: xintu(nr,nc+1), yintu(nr,nc+1)
REAL :: zerou1(nr,nc+1), zerou2(nr,nc+1)
REAL :: zerou3(nr,nc+1), zerou4(nr,nc+1)
INTEGER :: i, j, k
REAL :: ucor, uqse, uqsp
!REAL :: cv2atu, gx2u, gz2u

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *FSU*      REAL      X-COMPONENT OF SURFACE STRESS AT U-NODE [m2/s2]
!    *GZU*      REAL      VERTICAL GRID SPACING AT U-NODE
!    *ULHSA*    REAL      LEFT HAND SIDE OF U2-EQUATION /U2(I,J,K-1)
!    *ULHSB*    REAL      LEFT HAND SIDE OF U2-EQUATION /U2(I,J,K)
!    *ULHSC*    REAL      LEFT HAND SIDE OF U2-EQUATION /U2(I,J,K+1)
!    *URHS*     REAL      RIGHT HAND SIDE OF U2-EQUATION [m/s]
!    *U2A*      REAL      NEW VALUE OF U2 AFTER FIRST FRACTIONAL TIME STEP
!                                                                   [m/s]
!    *U2B*      REAL      NEW VALUE OF U2 AFTER SECOND FRACTIONAL TIME STEP
!                                                                   [m/s]
!    *VEDU*     REAL      VERTICAL EDDY COEFFICIENT BELOW U-NODE   [m2/s]
!    *W2U*      REAL      VERTICAL ADVECTIVE VELOCITY BELOW U-NODE  [m/s]
!    *XINTU*    REAL      DEPTH-INTEGRATED HORIZONTAL ADVECTION TERM FOR U2
!                         IN X-DIRECTION                          [m2/s2]
!    *YINTU*    REAL      DEPTH-INTEGRATED HORIZONTAL ADVECTION TERM FOR U2
!                         IN Y-DIRECTION                          [m2/s2]
!    *ZEROU1*   REAL      =0 DUMMY ARRAY FOR ZADVDIS ROUTINE
!    *ZEROU2*   REAL      =0 DUMMY ARRAY FOR ZADVDIS ROUTINE
!    *ZEROU3*   REAL      =0 DUMMY ARRAY FOR ZADVDIS ROUTINE
!    *ZEROU4*   REAL      =0 DUMMY ARRAY FOR ZADVDIS ROUTINE

!------------------------------------------------------------------------

!     1. INITIALISE ARRAYS
!     --------------------


!WRITE (*,*) 'TRACE: u(predicted) - internal points'

CALL zero33(zerou1,1,nr,nc+1)
CALL zero33(zerou2,1,nr,nc+1)
CALL zero33(zerou3,1,nr,nc+1)
CALL zero33(zerou4,1,nr,nc+1)

CALL zero33(urhs,nz,nr,nc+1)

DO  i=1,nc+1
  DO  j=1,nr
    DO  k=1,nz
      u2a(k,j,i) = u1(k,j,i)
      u2b(k,j,i) = u1(k,j,i)
    END DO
  END DO
END DO

!     2. INTERPOLATE ARRAYS AT U-NODES
!     --------------------------------

!     2.1 VERTICAL VELOCITY AND EDDY VISCOSITY
!     ----------------------------------------

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz+1
        w2u(k,j,i) = 0.5*(w2(k,j,i-1)+w2(k,j,i))
        vedu(k,j,i) = 0.5*(veddyv(k,j,i-1)+veddyv(k,j,i))
      END DO
    END IF
  END DO
END DO

!     2.2 GRID SPACING
!     ----------------


DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        gzu(k,j,i) = gz2u(k,j,i)
      END DO
    END IF
  END DO
END DO

!     2.3 SURFACE STRESS
!     ------------------

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      fsu(j,i) = 0.5*(fs(j,i-1)+fs(j,i))
    END IF
  END DO
END DO

!     3. DEPTH-INTEGRATED CORRECTION TERMS FOR 2-D MODE
!        ADVECTIVE/DIFFUSIVE TERMS FOR UPWIND SCHEME
!     -------------------------------------------------


CALL xhad3du(u1,u1,urhs,xintu,iodif)
CALL yhad3du(u1,v1,urhs,yintu,iodif)
CALL had2du
DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      uadhdev(j,i) = xintu(j,i) + yintu(j,i) + uah2d(j,i) - udh2d(j,i)
    END IF
  END DO
END DO

!     4. UPDATE U2 WITHOUT OPERATOR SPLITTING
!     ---------------------------------------


IF (iadvc <= 1) THEN
  
  
!     4.1 ASSEMBLE INTERIA TERMS
!     --------------------------
  
  DO  i=1,nc+1
    DO  j=1,nr
      DO  k=1,nz
        ulhsa(k,j,i) = 0.0
        ulhsb(k,j,i) = 1.0
        ulhsc(k,j,i) = 0.0
        IF (npix(j,i) == 0) THEN
          urhs(k,j,i) = 0.0
        ELSE
          urhs(k,j,i) = urhs(k,j,i) + u2(k,j,i)
        END IF
      END DO
    END DO
  END DO


!     4.2 ADD EXPLICIT TERMS
!     ----------------------
  
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) THEN
        DO  k=1,nz
          ucor = coriol(j)*cv2atu(k,j,i)
          uqse = -g*(zeta2(j,i)-zeta2(j,i-1))/gx2u(j,i)
          uqsp = -(p2(j,i)/r0ref-p2(j,i-1)/r0ref)/gx2u(j,i)
          urhs(k,j,i) = urhs(k,j,i) + del3*(ucor+uqse +uqsp+uqden(k,j,i))
        END DO
      END IF
    END DO
  END DO
  421   CONTINUE

!     4.3 ADD VERTICAL ADVECTION/DIFFUSION TERMS
!     ------------------------------------------
  
  CALL zadvdis(u2,vedu,w2u,gzu,npix,0,0,  &
      nz,nr,nc+1,fsu,zerou1,zerou2,zerou3,fbk,zerou4,  &
      iadvc,ulhsa,ulhsb,ulhsc,urhs)
 
!     4.4 INVERT TRIDIAGONAL SYSTEM TO OBTAIN U2
!     ------------------------------------------
  
  
  CALL thomv(ulhsa,ulhsb,ulhsc,urhs,u2,nz,nr,nc+1)
  
  GO TO 3000
  
END IF


!     5. EXPLICIT X-DIRECTION ADVECTION AND DIFFUSION - STEP A
!     --------------------------------------------------------
!     5.1 ASSEMBLE INTERIA TERMS - STEP A
!     -----------------------------------


DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        urhs(k,j,i) = u2(k,j,i)
      END DO
    END IF
  END DO
END DO
511   CONTINUE

!     5.2 ASSEMBLE X-ADVECTION/DIFFUSION TERMS - STEP A
!     -------------------------------------------------


CALL xhad3du(u2,u1,urhs,xintu,iodif)


!     5.3 UPDATE TO GET U2A-X
!     -----------------------


DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        u2a(k,j,i) = urhs(k,j,i)
      END DO
    END IF
  END DO
END DO
531   CONTINUE

!     6. EXPLICIT Y-DIRECTION ADVECTION AND DIFFUSION - STEP A
!     --------------------------------------------------------

!     6.1 ASSEMBLE Y-ADVECTION/DIFFUSION TERMS - STEP A
!     -------------------------------------------------


CALL yhad3du(u2a,v1,urhs,yintu,iodif)


!     6.2 UPDATE TO GET U2A - XY
!     --------------------------

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        u2a(k,j,i) = urhs(k,j,i)
      END DO
    END IF
  END DO
END DO
621   CONTINUE

!     7. EXPLICIT TERMS, IMPLICIT Z-DIR. ADVECTION/DIFFUSION TERMS - STEP A
!     ---------------------------------------------------------------------

!     7.1 PREPARE ARGUMENTS FOR CALL TO ZADVDIS - STEP A
!     --------------------------------------------------

DO  i=1,nc+1
  DO  j=1,nr
    DO  k=1,nz
      ulhsa(k,j,i) = 0.0
      ulhsb(k,j,i) = 1.0
      ulhsc(k,j,i) = 0.0
    END DO
  END DO
END DO


!     7.2 ADD EXPLICIT TERMS (CORIOLIS AND PRESSURE GRADIENT) - STEP A
!     ----------------------------------------------------------------

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        ucor = coriol(j)*cv2atu(k,j,i)
        uqse = -g*(zeta2(j,i)-zeta2(j,i-1))/gx2u(j,i)
        uqsp = -(p2(j,i)/r0ref-p2(j,i-1)/r0ref)/gx2u(j,i)
        urhs(k,j,i) = urhs(k,j,i) + del3*(ucor+uqse +uqsp+uqden(k,j,i))
      END DO
    END IF
  END DO
END DO
721   CONTINUE

!     7.3 EVALUATE VERTICAL ADVECTION AND DIFFUSION TERMS - STEP A
!     ------------------------------------------------------------


CALL zadvdis(u2a,vedu,w2u,gzu,npix,0,0,  &
    nz,nr,nc+1,fsu,zerou1,zerou2,zerou3,fbk,zerou4, iadvc,ulhsa,ulhsb,ulhsc,urhs)


!     7.4 INVERT TRIDIAGONAL SYSTEM TO OBTAIN U2A
!     -------------------------------------------


CALL thomv(ulhsa,ulhsb,ulhsc,urhs,u2a,nz,nr,nc+1)


!     8. EXPLICIT TERMS, IMPLICIT Z-DIR. ADVECTION/DIFFUSION TERMS - STEP B
!     ---------------------------------------------------------------------

!     8.1 PREPARE ARGUMENTS FOR CALL TO ZADVDIS - STEP B
!     --------------------------------------------------

DO  i=1,nc+1
  DO  j=1,nr
    DO  k=1,nz
      ulhsa(k,j,i) = 0.0
      ulhsb(k,j,i) = 1.0
      ulhsc(k,j,i) = 0.0
      IF (npix(j,i) == 1) THEN
        urhs(k,j,i) = u2(k,j,i)
      ELSE
        urhs(k,j,i) = 0.0
      END IF
    END DO
  END DO
END DO


!     8.2 ADD EXPLICIT TERMS (CORIOLIS AND PRESSURE GRADIENT) - STEP B
!     ----------------------------------------------------------------

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        ucor = coriol(j)*cv2atu(k,j,i)
        uqse = -g*(zeta2(j,i)-zeta2(j,i-1))/gx2u(j,i)
        uqsp = -(p2(j,i)/r0ref-p2(j,i-1)/r0ref)/gx2u(j,i)
        urhs(k,j,i) = urhs(k,j,i) + del3*(ucor+uqse +uqsp+uqden(k,j,i))
      END DO
    END IF
  END DO
END DO
821   CONTINUE

!     8.3 EVALUATE VERTICAL ADVECTION AND DIFFUSION TERMS - STEP B
!     ------------------------------------------------------------


CALL zadvdis(u2,vedu,w2u,gzu,npix,0,0,  &
    nz,nr,nc+1,fsu,zerou1,zerou2,zerou3,fbk,zerou4, iadvc,ulhsa,ulhsb,ulhsc,urhs)


!     8.4 INVERT TRIDIAGONAL SYSTEM TO OBTAIN U2B
!     -------------------------------------------


CALL thomv(ulhsa,ulhsb,ulhsc,urhs,u2b,nz,nr,nc+1)

!     9. EXPLICIT Y-DIRECTION ADVECTION AND DIFFUSION - STEP B
!     --------------------------------------------------------

!     9.1 ASSEMBLE INERTIA TERMS - STEP B
!     -----------------------------------


DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        urhs(k,j,i) = u2b(k,j,i)
      END DO
    END IF
  END DO
END DO
911   CONTINUE


!     9.2 ASSEMBLE Y-ADVECTION/DIFFUSION TERMS - STEP B
!     -------------------------------------------------


CALL yhad3du(u2b,v1,urhs,yintu,iodif)


!     9.3 UPDATE TO GET U2B - XY
!     --------------------------


DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        u2b(k,j,i) = urhs(k,j,i)
      END DO
    END IF
  END DO
END DO
931   CONTINUE

!     10. EXPLICIT X-DIRECTION ADVECTION AND DIFFUSION - STEP B
!     ---------------------------------------------------------

!     10.1 ASSEMBLE X-ADVECTION/DIFFUSION TERMS - STEP B
!     --------------------------------------------------


CALL xhad3du(u2b,u1,urhs,xintu,iodif)


!     10.2 UPDATE TO GET U2B
!     ----------------------

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        u2b(k,j,i) = urhs(k,j,i)
      END DO
    END IF
  END DO
END DO

!     11. TAKE AVERAGE OF U2A AND U2B
!     -------------------------------


DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        u2(k,j,i) = 0.5*(u2a(k,j,i)+u2b(k,j,i))
      END DO
    END IF
  END DO
END DO

3000 CONTINUE


RETURN

END SUBROUTINE ucalc

!========================================================================

SUBROUTINE vcalc
!************************************************************************

!    *VCALC*      SOLVE V-COMPONENT OF THE 3-D MOMENTUM EQUATION
!                 (PREDICTOR STEP)

!       AUTHOR - PATRICK LUYTEN AND KEVIN RUDDICK

!       LAST UPDATE - 13 Dec 1999 @(COHERENS)crrnt3.f 8.4

!       DESCRIPTION - SOLVES V2-EQUATION AT V-NODES
!                   - OPERATOR SPLITTING IS DISABLED IN THE ABSENCE OF
!                     ADVECTION OR FOR THE UPWIND SCHEME, AND IS ENABLED
!                     OTHERWISE

!       REFERENCE - Section III-4.3.1, 4.3.2 of the User Documentation

!       CALLING PROGRAM - CRRNT3P

!       EXTERNALS - THOMV, XHAD3DV, YHAD3DV, ZADVDIS, ZERO

!************************************************************************
!*    LOCAL VARIABLES

REAL :: gsv(nr+1,nc), gzv(nz,nr+1,nc), vedv(nz+1,nr+1,nc)
REAL :: vlhsa(nz,nr+1,nc), vlhsb(nz,nr+1,nc)
REAL :: vlhsc(nz,nr+1,nc), vrhs(nz,nr+1,nc)
REAL :: v2a(nz,nr+1,nc), v2b(nz,nr+1,nc)
REAL :: w2v(nz+1,nr+1,nc)
REAL :: xintv(nr+1,nc), yintv(nr+1,nc)
REAL :: zerov1(nr+1,nc), zerov2(nr+1,nc)
REAL :: zerov3(nr+1,nc), zerov4(nr+1,nc)
INTEGER :: i, j, k
REAL :: vcor, vqse, vqsp
!REAL :: cu2atv, gy2v, gz2v

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *GSV*      REAL      Y-COMPONENT OF SURFACE STRESS AT V-NODE [m2/s2]
!    *GZV*      REAL      VERTICAL GRID SPACING AT V-NODE
!    *VEDV*     REAL      VERTICAL EDDY COEFFICIENT BELOW V-NODE   [m2/s]
!    *VLHSA*    REAL      LEFT HAND SIDE OF V2-EQUATION /V2(I,J,K-1)
!    *VLHSB*    REAL      LEFT HAND SIDE OF V2-EQUATION /V2(I,J,K)
!    *VLHSC*    REAL      LEFT HAND SIDE OF V2-EQUATION /V2(I,J,K+1)
!    *VRHS*     REAL      RIGHT HAND SIDE OF V2-EQUATION            [m/s]
!    *V2A*      REAL      NEW VALUE OF V2 AFTER FIRST FRACTIONAL TIME STEP
!                                                                   [m/s]
!    *V2B*      REAL      NEW VALUE OF V2 AFTER SECOND FRACTIONAL TIME STEP
!                                                                   [m/s]
!    *W2V*      REAL      VERTICAL ADVECTIVE VELOCITY BELOW V-NODE  [m/s]
!    *XINTV*    REAL      DEPTH-INTEGRATED HORIZONTAL ADVECTION TERM FOR V2
!                         IN X-DIRECTION                          [m2/s2]
!    *YINTV*    REAL      DEPTH-INTEGRATED HORIZONTAL ADVECTION TERM FOR V2
!                         IN Y-DIRECTION                          [m2/s2]
!    *ZEROV1*   REAL      =0 DUMMY ARRAY FOR ZADVDIS ROUTINE
!    *ZEROV2*   REAL      =0 DUMMY ARRAY FOR ZADVDIS ROUTINE
!    *ZEROV3*   REAL      =0 DUMMY ARRAY FOR ZADVDIS ROUTINE
!    *ZEROV4*   REAL      =0 DUMMY ARRAY FOR ZADVDIS ROUTINE

!------------------------------------------------------------------------

!     1. INITIALISE ARRAYS
!     --------------------


!WRITE (*,*) 'TRACE: v(predicted) - internal points'

CALL zero33(zerov1,1,nr+1,nc)
CALL zero33(zerov2,1,nr+1,nc)
CALL zero33(zerov3,1,nr+1,nc)
CALL zero33(zerov4,1,nr+1,nc)

CALL zero33(vrhs,nz,nr+1,nc)

DO  i=1,nc
  DO  j=1,nr+1
    DO  k=1,nz
      v2a(k,j,i) = v1(k,j,i)
      v2b(k,j,i) = v1(k,j,i)
    END DO
  END DO
END DO

!     2. INTERPOLATE ARRAYS AT V-NODES
!     --------------------------------

!     2.1 VERTICAL VELOCITY AND EDDY VISCOSITY
!     ----------------------------------------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz+1
        w2v(k,j,i) = 0.5*(w2(k,j-1,i)+w2(k,j,i))
        vedv(k,j,i) = 0.5*(veddyv(k,j-1,i)+veddyv(k,j,i))
      END DO
    END IF
  END DO
END DO

!     2.2 GRID SPACING
!     ----------------


DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        gzv(k,j,i) = gz2v(k,j,i)
      END DO
    END IF
  END DO
END DO

!     2.3 SURFACE STRESS
!     ------------------


DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      gsv(j,i) = 0.5*(gs(j-1,i)+gs(j,i))
    END IF
  END DO
END DO

!     3. DEPTH-INTEGRATED CORRECTION TERMS FOR 2-D MODE
!        ADVECTIVE/DIFFUSIVE TERMS FOR UPWIND SCHEME
!     -------------------------------------------------


CALL xhad3dv(v1,u1,vrhs,xintv,iodif)

CALL yhad3dv(v1,v1,vrhs,yintv,iodif)

CALL had2dv

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      vadhdev(j,i) = xintv(j,i) + yintv(j,i) + vah2d(j,i) - vdh2d(j,i)
    END IF
  END DO
END DO

!     4. UPDATE V2 WITHOUT OPERATOR SPLITTING
!     ---------------------------------------

IF (iadvc <= 1) THEN
  
  
!     4.1 ASSEMBLE INTERIA TERMS
!     --------------------------
  
  DO  i=1,nc
    DO  j=1,nr+1
      DO  k=1,nz
        vlhsa(k,j,i) = 0.0
        vlhsb(k,j,i) = 1.0
        vlhsc(k,j,i) = 0.0
        IF (npiy(j,i) == 0) THEN
          vrhs(k,j,i) = 0.0
        ELSE
          vrhs(k,j,i) = vrhs(k,j,i) + v2(k,j,i)
        END IF
      END DO
    END DO
  END DO
!     4.2 ADD EXPLICIT TERMS
!     ----------------------
  
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        DO  k=1,nz
  !        vcor = -coriolv(j)*(cu2atv(k,j,i)-ugeo)
           vcor = -coriolv(j)*cu2atv(k,j,i)
  !        zzz = h2atv(j,i)*(1.0-gz0(k,j,i))
  !        if((zzz>200.0).and.(zzz<500.0).and.(i<1000)) vcor = vcor-ad*0.3e-5

          vqse = -g*(zeta2(j,i)-zeta2(j-1,i))/gy2v(j)
          vqsp = -(p2(j,i)/r0ref-p2(j-1,i)/r0ref)/gy2v(j)
          vrhs(k,j,i) = vrhs(k,j,i) + del3*(vcor+vqse +vqsp+vqden(k,j,i))
        END DO
      END IF
    END DO
  END DO

  421   CONTINUE
  
!     4.3 ADD VERTICAL ADVECTION/DIFFUSION TERMS
!     ------------------------------------------
  
  CALL zadvdis(v2,vedv,w2v,gzv,npiy,0,0,  &
      nz,nr+1,nc,gsv,zerov1,zerov2,zerov3,gbk,zerov4,  &
      iadvc,vlhsa,vlhsb,vlhsc,vrhs)
  
  
!     4.4 INVERT TRIDIAGONAL SYSTEM TO OBTAIN U2
!     ------------------------------------------

  CALL thomv(vlhsa,vlhsb,vlhsc,vrhs,v2,nz,nr+1,nc)

  GO TO 3000
  
END IF

!     5. EXPLICIT X-DIRECTION ADVECTION AND DIFFUSION - STEP A
!     --------------------------------------------------------

!     5.1 ASSEMBLE INTERIA TERMS - STEP A
!     -----------------------------------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        vrhs(k,j,i) = v2(k,j,i)
      END DO
    END IF
  END DO
END DO
511   CONTINUE

!     5.2 ASSEMBLE X-ADVECTION/DIFFUSION TERMS - STEP A
!     -------------------------------------------------


CALL xhad3dv(v2,u1,vrhs,xintv,iodif)


!     5.3 UPDATE TO GET V2A-X
!     -----------------------


DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        v2a(k,j,i) = vrhs(k,j,i)
      END DO
    END IF
  END DO
END DO
531   CONTINUE


!     6. EXPLICIT Y-DIRECTION ADVECTION AND DIFFUSION - STEP A
!     --------------------------------------------------------

!     6.1 ASSEMBLE Y-ADVECTION/DIFFUSION TERMS - STEP A
!     -------------------------------------------------

CALL yhad3dv(v2a,v1,vrhs,yintv,iodif)


!     6.2 UPDATE TO GET V2A - XY
!     --------------------------


DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        v2a(k,j,i) = vrhs(k,j,i)
      END DO
    END IF
  END DO
END DO
621   CONTINUE


!     7. EXPLICIT TERMS, IMPLICIT Z-DIR. ADVECTION/DIFFUSION TERMS - STEP A
!     ---------------------------------------------------------------------

!     7.1 PREPARE ARGUMENTS FOR CALL TO ZADVDIS - STEP A
!     --------------------------------------------------


DO  i=1,nc
  DO  j=1,nr+1
    DO  k=1,nz
      vlhsa(k,j,i) = 0.0
      vlhsb(k,j,i) = 1.0
      vlhsc(k,j,i) = 0.0
    END DO
  END DO
END DO


!     7.2 ADD EXPLICIT TERMS (CORIOLIS AND PRESSURE GRADIENT) - STEP A
!     ----------------------------------------------------------------


DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
!        vcor = -coriolv(j)*(cu2atv(k,j,i)-ugeo)
           vcor = -coriolv(j)*cu2atv(k,j,i)
 !         zzz = h2atv(j,i)*(1.0-gz0(k,j,i))
 !         if((zzz>200.0).and.(zzz<500.0).and.(i<1000)) vcor = vcor-0.3e-5

        vqse = -g*(zeta2(j,i)-zeta2(j-1,i))/gy2v(j)
        vqsp = -(p2(j,i)/r0ref-p2(j-1,i)/r0ref)/gy2v(j)
        vrhs(k,j,i) = vrhs(k,j,i) + del3*(vcor+vqse +vqsp+vqden(k,j,i))
      END DO
    END IF
  END DO
END DO
721   CONTINUE


!     7.3 EVALUATE VERTICAL ADVECTION AND DIFFUSION TERMS - STEP A
!     ------------------------------------------------------------


CALL zadvdis(v2a,vedv,w2v,gzv,npiy,0,0,  &
    nz,nr+1,nc,gsv,zerov1,zerov2,zerov3,gbk,zerov4, iadvc,vlhsa,vlhsb,vlhsc,vrhs)


!     7.4 INVERT TRIDIAGONAL SYSTEM TO OBTAIN V2A
!     -------------------------------------------

CALL thomv(vlhsa,vlhsb,vlhsc,vrhs,v2a,nz,nr+1,nc)


!     8. EXPLICIT TERMS, IMPLICIT Z-DIR. ADVECTION/DIFFUSION TERMS - STEP B
!     ---------------------------------------------------------------------

!     8.1 PREPARE ARGUMENTS FOR CALL TO ZADVDIS - STEP B
!     --------------------------------------------------

DO  i=1,nc
  DO  j=1,nr+1
    DO  k=1,nz
      vlhsa(k,j,i) = 0.0
      vlhsb(k,j,i) = 1.0
      vlhsc(k,j,i) = 0.0
      IF (npiy(j,i) == 1) THEN
        vrhs(k,j,i) = v2(k,j,i)
      ELSE
        vrhs(k,j,i) = 0.0
      END IF
    END DO
  END DO
END DO

!     8.2 ADD EXPLICIT TERMS (CORIOLIS AND PRESSURE GRADIENT) - STEP B
!     ----------------------------------------------------------------


DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
!        vcor = -coriolv(j)*(cu2atv(k,j,i)-ugeo)
        vcor = -coriolv(j)*cu2atv(k,j,i)
!          zzz = h2atv(j,i)*(1.0-gz0(k,j,i))
!          if((zzz>200.0).and.(zzz<500.0).and.(i<1000)) vcor = vcor-ad*0.3e-5

        vqse = -g*(zeta2(j,i)-zeta2(j-1,i))/gy2v(j)
        vqsp = -(p2(j,i)/r0ref-p2(j-1,i)/r0ref)/gy2v(j)
        vrhs(k,j,i) = vrhs(k,j,i) + del3*(vcor+vqse +vqsp+vqden(k,j,i))
      END DO
    END IF
  END DO
END DO
821   CONTINUE


!     8.3 EVALUATE VERTICAL ADVECTION AND DIFFUSION TERMS - STEP B
!     ------------------------------------------------------------


CALL zadvdis(v2,vedv,w2v,gzv,npiy,0,0,  &
    nz,nr+1,nc,gsv,zerov1,zerov2,zerov3,gbk,zerov4, iadvc,vlhsa,vlhsb,vlhsc,vrhs)


!     8.4 INVERT TRIDIAGONAL SYSTEM TO OBTAIN V2B
!     -------------------------------------------


CALL thomv(vlhsa,vlhsb,vlhsc,vrhs,v2b,nz,nr+1,nc)

!     9. EXPLICIT Y-DIRECTION ADVECTION AND DIFFUSION - STEP B
!     --------------------------------------------------------

!     9.1 ASSEMBLE INERTIA TERMS - STEP B
!     -----------------------------------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        vrhs(k,j,i) = v2b(k,j,i)
      END DO
    END IF
  END DO
END DO
911   CONTINUE


!     9.2 ASSEMBLE Y-ADVECTION/DIFFUSION TERMS - STEP B
!     -------------------------------------------------


CALL yhad3dv(v2b,v1,vrhs,yintv,iodif)


!     9.3 UPDATE TO GET V2B - XY
!     --------------------------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        v2b(k,j,i) = vrhs(k,j,i)
      END DO
    END IF
  END DO
END DO
931   CONTINUE

!     10. EXPLICIT X-DIRECTION ADVECTION AND DIFFUSION - STEP B
!     ---------------------------------------------------------

!     10.1 ASSEMBLE X-ADVECTION/DIFFUSION TERMS - STEP B
!     --------------------------------------------------


CALL xhad3dv(v2b,u1,vrhs,xintv,iodif)


!     10.2 UPDATE TO GET V2B
!     ----------------------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        v2b(k,j,i) = vrhs(k,j,i)
      END DO
    END IF
  END DO
END DO


!     11. TAKE AVERAGE OF V2A AND V2B
!     -------------------------------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        v2(k,j,i) = 0.5*(v2a(k,j,i)+v2b(k,j,i))
      END DO
    END IF
  END DO
END DO


3000  CONTINUE

RETURN

END SUBROUTINE vcalc

!=======================================================================

SUBROUTINE wcalc
!************************************************************************

!    *WCALC*      SOLVE THE 3-D CONTINUITY EQUATION FOR THE TRANSFORMED
!                 VERTICAL VELOCITY

!       AUTHOR - KEVIN RUDDICK

!       LAST UPDATE - 26 Mar 1999 @(COHERENS)crrnt3.f 8.4

!       DESCRIPTION -

!       REFERENCE - Section III-4.3.1d of the User Documentation

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS - WPCALC

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, ii, iwerr, j, jj, jwerr, k
REAL :: fluxu(nz,nr,nc+1), fluxv(nz,nr+1,nc)
REAL :: werr, w2divu, w2divv
!REAL :: gz2u, gz2v, h2atu, h2atv, h2atc

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *FLUXU*     REAL     DEVIATION OF MASS FLUX AT U FACES [m2/s]
!    *FLUXV*     REAL     DEVIATION OF MASS FLUX AT U FACES [m2/s]

!------------------------------------------------------------------------

!     1. INITIALISE DEVIATION OF MASS FLUX ARRAYS
!     -------------------------------------------


!WRITE (*,*) 'TRACE: calculating wJ'

CALL zero33(fluxu,nz,nr,nc+1)
CALL zero33(fluxv,nz,nr+1,nc)


!     2. CALCULATE DEVIATION OF MASS FLUXES AT U FACES
!     ------------------------------------------------

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        fluxu(k,j,i) = gz2u(k,j,i)*(u2f(k,j,i)-ud2f(j,i) /h2atu(j,i))
      END DO
    END IF
  END DO
END DO


!     3. CALCULATE DEVIATION OF MASS FLUXES AT V FACES
!     ------------------------------------------------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        fluxv(k,j,i) = gz2v(k,j,i)*cosphiv(j)*  &
            (v2f(k,j,i)-vd2f(j,i)/h2atv(j,i))
      END DO
    END IF
  END DO
END DO


!     4. CALCULATE DEVIATION OF MASS FLUXES AT U OPEN BOUNDARIES
!     ----------------------------------------------------------

DO  ii=1,nobu
  i=iobu(ii)
  j=jobu(ii)
  IF (westob(ii)) THEN
    DO  k=1,nz
      fluxu(k,j,i) = gz2(k,j,i)*(u2f(k,j,i)-ud2f(j,i) /h2atc(j,i))
    END DO
  ELSE
    DO  k=1,nz
      fluxu(k,j,i) = gz2(k,j,i-1)*(u2f(k,j,i) - ud2f(j,i)/h2atc(j,i-1))
    END DO
  END IF
END DO


!     5. CALCULATE DEVIATION OF MASS FLUXES AT V OPEN BOUNDARIES
!     ----------------------------------------------------------

DO  jj=1,nobv
  i=iobv(jj)
  j=jobv(jj)
  IF (soutob(jj)) THEN
    DO  k=1,nz
      fluxv(k,j,i) = gz2(k,j,i)*cosphiv(j)* (v2f(k,j,i)-vd2f(j,i)/h2atc(j,i))
    END DO
  ELSE
    DO  k=1,nz
      fluxv(k,j,i) = gz2(k,j-1,i)*cosphiv(j)*  &
          (v2f(k,j,i)-vd2f(j,i)/h2atc(j-1,i))
    END DO
  END IF
END DO


!     6. COMBINE ALL MASS FLUXES TO GIVE TRANSFORMED VERTICAL VELOCITY
!     ----------------------------------------------------------------

werr = 0.0
DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      w2(1,j,i) = 0.0
      DO  k=1,nz
        w2divu = (fluxu(k,j,i+1)-fluxu(k,j,i))/gx2(j,i)
        w2divv = (fluxv(k,j+1,i)-fluxv(k,j,i)) /(gy2(j)*cosphi(j))
        w2(k+1,j,i) = w2(k,j,i) - w2divu - w2divv
      END DO
!         --note error in vertical velocity at surface
      IF (ABS(w2(nz+1,j,i)) > werr) THEN
        werr = ABS(w2(nz+1,j,i))
        iwerr= i
        jwerr= j
      END IF
      w2(nz+1,j,i) = 0.0
    END IF
  END DO
END DO


!     7. UPDATE "PHYSICAL" VELOCITY
!     -----------------------------


CALL wpcalc

!     ---- output vertical velocity non-conservativity
!WRITE (*,9101) 'W',werr*100,' cm/s',iwerr,jwerr


RETURN

9101 FORMAT (a1,'ERR=',1PE14.7,a5,' (',i3,',',i3,')')

END SUBROUTINE wcalc

!=======================================================================

SUBROUTINE wpcalc
!************************************************************************

!    *WPCALC*      CALCULATE "PHYSICAL" VERTICAL VELOCITY

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 26 Mar 1999 @(COHERENS)crrnt3.f 8.4

!       DESCRIPTION -

!       REFERENCE - Section III-4.3.1d of the User Documentation

!       CALLING PROGRAM - WCALC

!       EXTERNALS -

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, ii, j, jj, k
REAL :: u2jz(nz,nr,nc+1), v2jz(nz,nr+1,nc), w2jz(nz+1,nr,nc)
REAL :: zjnew(nz,nr,nc), zjold(nz,nr,nc)
!REAL :: depatu, depatv, gz0c, gz2u, gz2v
!REAL :: h1atc, h2atc, h2atu, h2atv

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *U2JZ*     REAL      VERTICAL COORDINATE TIMES J AND U2 AT U-NODES
!    *V2JZ*     REAL      VERTICAL COORDINATE TIMES J AND V2 AT V-NODES
!    *W2JZ*     REAL      VERTICAL COORDINATE TIMES J AND TRANSFORMED
!                         VERTICAL VELOCITY AT W-NODES
!    *ZJNEW*    REAL      VERTICAL COORDINATE TIMES J AT NEW TIME STEP
!    *ZJOLD*    REAL      VERTICAL COORDINATE TIMES J AT OLD TIME STEP

!-----------------------------------------------------------------------

!     1. INITIALISE ARRAYS
!     --------------------

CALL zero33(zjold,nz,nr,nc)
CALL zero33(zjnew,nz,nr,nc)
CALL zero33(u2jz,nz,nr,nc+1)
CALL zero33(v2jz,nz,nr+1,nc)
CALL zero33(w2jz,nz+1,nr,nc)


!     2. Z TIMES J
!     ------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        zjold(k,j,i) = (gz0c(k,j,i)*h1atc(j,i)-dep(j,i))*gz1(k,j,i)
        zjnew(k,j,i) = (gz0c(k,j,i)*h2atc(j,i)-dep(j,i))*gz2(k,j,i)
      END DO
    END IF
  END DO
END DO
210   CONTINUE


!     3. Z TIMES U2 AND J
!     -------------------

!     3.1 INTERIOR POINTS
!     -------------------

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        u2jz(k,j,i) = (gz0c(k,j,i)*h2atu(j,i)-depatu(j,i)) *gz2u(k,j,i)*u2(k,j,i)
      END DO
    END IF
  END DO
END DO
310   CONTINUE


!     3.2 U-OPEN BOUNDARIES
!     ---------------------

DO  ii=1,nobu
  i = iobu(ii)
  j = jobu(ii)
  DO  k=1,nz
    IF (westob(ii)) THEN
      u2jz(k,j,i) = (gz0c(k,j,i)*h2atc(j,i)-dep(j,i)) *gz2(k,j,i)*u2(k,j,i)
    ELSE
      u2jz(k,j,i) = (gz0c(k,j,i)*h2atc(j,i-1)-dep(j,i-1)) *gz2(k,j,i-1)*u2(k,j,i)
    END IF
  END DO
END DO


!     4. Z TIMES V2 AND J
!     -------------------

!     4.1 INTERIOR POINTS
!     -------------------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        v2jz(k,j,i) = (gz0c(k,j,i)*h2atv(j,i)-depatv(j,i))  &
            *gz2v(k,j,i)*cosphiv(j)*v2(k,j,i)
      END DO
    END IF
  END DO
END DO
410   CONTINUE


!     4.2 V-OPEN BOUNDARIES
!     ---------------------

DO  jj=1,nobv
  i = iobv(jj)
  j = jobv(jj)
  DO  k=1,nz
    IF (soutob(jj)) THEN
      v2jz(k,j,i) = (gz0c(k,j,i)*h2atc(j,i)-dep(j,i))  &
          *gz2(k,j,i)*cosphiv(j)*v2(k,j,i)
    ELSE
      v2jz(k,j,i) = (gz0c(k,j,i)*h2atc(j-1,i)-dep(j-1,i))  &
          *gz2(k,j-1,i)*cosphiv(j)*v2(k,j,i)
    END IF
  END DO
END DO


!     5. Z TIMES W2
!     -------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        w2jz(k,j,i) = (gz0(k,j,i)*h2atc(j,i)-dep(j,i))*w2(k,j,i)
      END DO
    END IF
  END DO
END DO
510   CONTINUE


!     6. EVALUATE VERTICAL PHYSICAL VELOCITY
!     --------------------------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        w2phys(k,j,i) = ((zjnew(k,j,i)  - zjold(k,j,i))/del3  &
            + (u2jz(k,j,i+1) - u2jz(k,j,i))/gx2(j,i)  &
            + (v2jz(k,j+1,i) - v2jz(k,j,i))/(gy2(j)*cosphi(j))  &
            +  w2jz(k+1,j,i) - w2jz(k,j,i))/gz2(k,j,i)
      END DO
    END IF
  END DO
END DO
610   CONTINUE


RETURN

END SUBROUTINE wpcalc

SUBROUTINE dissip
!************************************************************************

!    *DISSIP*      SOLVE TRANSPORT EQUATION FOR DISSIPATION RATE

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 6 Nov 1998        @(COHERENS)dissip.f 8.4

!       DESCRIPTION - APPLY SURFACE AND BOTTOM BOUNDARY CONDITIONS
!                   - EVALUATE SOURCE/SINK TERMS
!                   - UPDATE DISSIPATION BY CALLING TRANSPW

!       REFERENCE - Section III-1.2.2b and III-4.5.3 of the User Documentation

!       CALLING PROGRAM - VEDDY2

!       EXTERNALS - TRANSPW

!************************************************************************
!*    LOCAL VARIABLES

LOGICAL :: call1, flag
INTEGER :: i, iadvt, ibot, idift, isur, j, k
REAL :: bdiss(nr,nc), sdiss(nr,nc)
REAL :: sink(nz+1,nr,nc), source(nz+1,nr,nc)
REAL :: zl1, zl2
!REAL :: h2atc
SAVE call1
DATA call1 /.true./

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL
!    *FLAG*     LOGICAL   .TRUE. IF TRANSPW IS CALLED FOR THE FIRST TIME
!    *IADVT*    INTEGER   TYPE OF SCHEME FOR ADVECTION (0/1/2/3/4)
!    *IBOT*     INTEGER   =2 (DIRICHLET NEAR BOTTOM B.C.)
!    *IDIFT*    INTEGER   TYPE OF SCHEME FOR HORIZONTAL DIFFUSION (0/1/2)
!    *ISUR*     INTEGER   =2 (DIRICHLET NEAR SURFACE B.C.)
!    *BDISS*    REAL      NEAR BOTTOM VALUE OF DISSW
!    *SDISS*    REAL      NEAR SURFACE VALUE OF DISSW
!    *SINK*     REAL      SINK TERMS IN DISSW EQUATION
!    *SOURCE*   REAL      SOURCE TERMS IN DISSW EQUATION

!------------------------------------------------------------------------

!     1. INITIALISE PARAMETERS AND ARRAYS
!     -----------------------------------

!WRITE (*,*) 'TRACE: dissip'

isur = 2
ibot = 2
IF (iahdht == 0) THEN
  iadvt = 0
  idift = 0
ELSE
  iadvt = iadvs
  idift = iodif
END IF

flag = call1
call1 = .false.


!     2. BOUNDARY CONDITIONS (SURFACE, BOTTOM)
!     ----------------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      zl2 = ckar*(h2atc(j,i)*(1.0-gz0(nz,j,i))+z0sur)
      zl1 = ckar*(h2atc(j,i)*gz0(2,j,i)+z0bot)
      sdiss(j,i) = eps0*tkew(nz,j,i)**1.5/zl2
      bdiss(j,i) = eps0*tkew(2,j,i)**1.5/zl1
    END IF
  END DO
END DO


!     3. SOURCE/SINK TERMS
!     --------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      source(2,j,i) = 0.0
      sink(2,j,i) = 0.0
      DO  k=3,nz-1
        IF (tkeold(k,j,i) > 0.0) THEN
          source(k,j,i) = c1e*shprod(k,j,i) *(dissw(k,j,i)/tkeold(k,j,i))
          sink(k,j,i) = c2e*dissw(k,j,i)/tkeold(k,j,i)
          IF (buprod(k,j,i) > 0.0) THEN
            source(k,j,i) = source(k,j,i) +c1e*c3e2*buprod(k,j,i)  &
                *dissw(k,j,i)/tkeold(k,j,i)
          ELSE
            sink(k,j,i) = sink(k,j,i)-c1e*c3e1*buprod(k,j,i) /tkeold(k,j,i)
          END IF
        ELSE
          source(k,j,i) = 0.0
          sink(k,j,i)   = 0.0
        END IF
      END DO
      source(nz,j,i) = 0.0
      sink(nz,j,i) = 0.0
    END IF
  END DO
END DO


!     4. UPDATE DISSW
!     ---------------


CALL transpw(dissw,veddye,source,sink, flag,isur,ibot,iadvt,iadvt,idift,  &
    ivpobu,ivpobv,turvp, zeros1,zeros2,sdiss,  &
    zeros4,zeros5,bdiss)


RETURN

END SUBROUTINE dissip

!=======================================================================


SUBROUTINE xhad3du (u2phi,uadv,sexp,xintu,jdifc)

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:39:32

!************************************************************************

!    *XHAD3DU*    EVALUATE X-DIRECTION ADVECTION-DIFFUSION TERMS IN
!                 U-MOMENTUM EQUATION INCLUDING DEPTH INTEGRALS FOR
!                 UD2-EQUATION

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 26 Mar 1999 @(COHERENS)hadvdis.f 8.4

!       DESCRIPTION - COMPUTES X-DIRECTION HORIZONTAL ADVECTION/DIFFUSION
!                     TERMS FOR U2
!                   - ADVECTIVE/DIFFUSIVE TERMS ARE EVALUATED EXPLICITLY
!                   - ADVECTION SCHEME SELECTED BY IADVC
!                      = 1 => FIRST ORDER UPWIND
!                      = 2 => LAX-WENDROFF
!                      = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                      = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!                   - DIFFUSION SELECTED BY JDIFC
!                      = 0 => OFF
!                      = 1 => UNIFORM
!                      = 2 => SMAGORINSKY

!       REFERENCES - Sections III-4.3.3a and 4.3.5a of the User Documentation

!       CALLING PROGRAM - UCALC

!       EXTERNALS - FNLIM, ZERO

!************************************************************************
!*    ARGUMENTS


REAL, INTENT(IN)                     :: u2phi(nz,nr,nc+1)
REAL, INTENT(IN)                         :: uadv(nz,nr,nc+1)
REAL, INTENT(IN OUT)                        :: sexp(nz,nr,nc+1)
REAL, INTENT(IN OUT)                        :: xintu(nr,nc+1)
INTEGER, INTENT(IN)                  :: jdifc


!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *U2PHI*    REAL      CURRENT TO BE ADVECTED/DIFFUSED [m/s]        (IN)
!    *UADV*     REAL      ADVECTING VELOCITY [m/s]                     (IN)
!    *SEXP*     REAL      RIGHT HAND SIDE OF U2-EQUATION (EXPLICIT TERMS)
!                                                        [m/s]      (IN/OUT)
!    *XINTU*    REAL      DEPTH-INTEGRATED HORIZONTAL ADVECTION TERM FOR U2
!                         IN X-DIRECTION                 [m2/s2]       (OUT)
!    *JDIFC*    INTEGER   SWITCH TO SELECT HORIZONTAL DIFFUSION (0/1/2) (IN)

!------------------------------------------------------------------------

!*    LOCAL VARIABLES

INTEGER :: i, j, k
REAL :: advfl(nz,nr,nc), advflw(nz,nr,0:nc+1),uadv2(nz,nr,0:nc+1)
REAL :: advfup(nz,nr,0:nc+1), xdifluj(nz,nr,nc)
REAL :: cfl, psi, uadvc, uwt, xadvu, xdifu
!REAL :: fnlim, gx2u, gz2u, v1atc

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ADVFL*    REAL      TVD ADVECTIVE FLUX OF U2 IN X-DIRECTION AT CELL
!                         CENTRE                                     [m2/s2]
!    *ADVFLW*   REAL      LW  ADVECTIVE FLUX OF U2 IN X-DIRECTION AT CELL
!                         CENTRE                                     [m2/s2]
!    *ADVFUP*   REAL      UPWIND ADVECTIVE FLUX OF U2 IN X-DIRECTION AT CELL
!                         CENTRE                                     [m2/s2]
!    *XDIFLUJ*  REAL      (1,1)-COMPONENT OF DIFFUSIVE STRESS AT CELL CENTRE
!                         TIMES J                                    [m2/s2]
!    *CFL*      REAL      COURANT-FRIEDRICH-LEVY NUMBER
!    *PSI*      REAL      FLUX LIMITER

!------------------------------------------------------------------------

!     1. INITIALISATION OF TEMPORARY ARRAYS
!     -------------------------------------


CALL zero33(advfl,nz,nr,nc)
CALL zero33(advflw,(nc+2)*nr*nz,1,1)
CALL zero33(advfup,(nc+2)*nr*nz,1,1)
CALL zero33(xdifluj,nz,nr,nc)
CALL zero33(xintu,1,nr,nc+1)

DO  i=1,nc+1
DO  j=1,nr
DO  k=1,nz
  uadv2(k,j,i) = uadv(k,j,i)+ugeo
END DO
END DO
END DO

!     2. ADVECTIVE FLUXES AT CENTRES
!     ------------------------------

IF (iadvc > 0) THEN
  
  
!     2.1 UPWIND AND LW FLUXES
!     ------------------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
          uadvc = 0.5*(uadv2(k,j,i)+uadv2(k,j,i+1))
          uwt = SIGN(1.0,uadvc)
          cfl = uadvc*del3/gx2(j,i)
          advfup(k,j,i) = 0.5*uadvc  &
              *((1.0+uwt)*u2phi(k,j,i)+(1.0-uwt)*u2phi(k,j,i+1))
          advflw(k,j,i) = 0.5*uadvc  &
              *((1.0+cfl)*u2phi(k,j,i)+(1.0-cfl)*u2phi(k,j,i+1))
        END DO
      END IF
    END DO
  END DO
  211   CONTINUE
  
  
!     2.2 ADVECTIVE FLUXES
!     --------------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
          psi = 0.0
          IF (iadvc == 1) THEN
!          ---fully upwind
            psi = 0.0
          ELSE IF (iadvc == 2) THEN
!          ---Lax-Wendroff
            psi = 1.0
          ELSE IF ((iadvc == 3).OR.(iadvc == 4)) THEN
!          ---TVD scheme
            uadvc = 0.5*(uadv2(k,j,i)+uadv2(k,j,i+1))
            IF (uadvc > 0.0) THEN
              psi = fnlim(advflw(k,j,i),  advfup(k,j,i),  &
                  advflw(k,j,i-1),advfup(k,j,i-1),iadvc)
            ELSE IF (uadvc < 0.0) THEN
              psi = fnlim(advflw(k,j,i),  advfup(k,j,i),  &
                  advflw(k,j,i+1),advfup(k,j,i+1),iadvc)
            END IF
          END IF
          advfl(k,j,i) = (advfup(k,j,i) +  &
              psi*(advflw(k,j,i)-advfup(k,j,i)))*gz2(k,j,i)
        END DO
      END IF
    END DO
  END DO
  221   CONTINUE
  
END IF


!     3. DIFFUSIVE FLUXES
!     -------------------

IF (jdifc > 0) THEN
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
          xdifluj(k,j,i) = 2.0*gz2(k,j,i)*heddyvc(k,j,i)*  &
              ((u1(k,j,i+1)-u1(k,j,i))/gx2(j,i)- sphcur(j)*v1atc(k,j,i))
        END DO
      END IF
    END DO
  END DO
  301   CONTINUE
  
END IF


!     4. ADVECTIVE AND DIFFUSIVE TERMS AT U-NODES
!     -------------------------------------------

IF ((iadvc > 0).OR.(jdifc > 0)) THEN
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) THEN
        DO  k=1,nz
          xadvu = (advfl(k,j,i)-advfl(k,j,i-1))/gx2u(j,i)
          xdifu = (xdifluj(k,j,i)-xdifluj(k,j,i-1))/gx2u(j,i)
          xintu(j,i) = xintu(j,i) - xadvu + xdifu
          sexp(k,j,i) = sexp(k,j,i) - (xadvu-xdifu)*del3/gz2u(k,j,i)
        END DO
      END IF
    END DO
  END DO
  401   CONTINUE
END IF


RETURN

END SUBROUTINE xhad3du

!========================================================================

SUBROUTINE yhad3du (u2phi,vadv,sexp,yintu,jdifc)
!************************************************************************

!    *YHAD3DU*    EVALUATE Y-DIRECTION ADVECTION-DIFFUSION TERMS IN
!                 U-MOMENTUM EQUATION INCLUDING DEPTH INTEGRALS FOR
!                 UD2-EQUATION

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 19 Apr 1999 @(COHERENS)hadvdis.f 8.4

!       DESCRIPTION - COMPUTES Y-DIRECTION HORIZONTAL ADVECTION/DIFFUSION
!                     TERMS FOR U2
!                   - ADVECTIVE/DIFFUSIVE TERMS ARE EVALUATED EXPLICITLY
!                   - ADVECTIVE FLUXES ARE EVALUATED AT V-NODES USING EITHER
!                     U2-VALUES WEST OR EAST OF V-NODE
!                   - ADVECTIVE TERMS ARE OBTAINED BY THE AVERAGE OF
!                     THE WEST/EAST FLUXES
!                   - DIFFUSIVE FLUXES USING X-(Y-)DERIVATIVES ARE EVALUATED
!                     AT THE U-(V-)NODES
!                   - ADVECTION SCHEME SELECTED BY IADVC
!                      = 1 => FIRST ORDER UPWIND
!                      = 2 => LAX-WENDROFF
!                      = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                      = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!                   - DIFFUSION SELECTED BY JDIFC
!                      = 0 => OFF
!                      = 1 => UNIFORM
!                      = 2 => SMAGORINSKY

!       REFERENCES - Sections III-4.3.3c and 4.3.5a of the User Documentation

!       CALLING PROGRAM - UCALC

!       EXTERNALS - FNLIM, ZERO

!************************************************************************
!*    ARGUMENTS


REAL, INTENT(IN)                         :: u2phi(nz,nr,nc+1)
REAL, INTENT(IN)                         :: vadv(nz,nr+1,nc)
REAL, INTENT(IN OUT)                        :: sexp(nz,nr,nc+1)
REAL, INTENT(IN OUT)                        :: yintu(nr,nc+1)
INTEGER, INTENT(IN)                  :: jdifc




!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *U2PHI*    REAL      CURRENT TO BE ADVECTED/DIFFUSED [m/s]         (IN)
!    *VADV*     REAL      ADVECTING VELOCITY [m/s]                      (IN)
!    *SEXP*     REAL      RIGHT HAND SIDE OF U2-EQUATION (EXPLICIT TERMS)
!                                                        [m/s]      (IN/OUT)
!    *YINTU*  REAL        DEPTH-INTEGRATED HORIZONTAL ADVECTION TERM FOR U2
!                         IN Y-DIRECTION                 [m2/s2]       (OUT)
!    *JDIFC*    INTEGER   SWITCH TO SELECT HORIZONTAL DIFFUSION (0/1/2) (IN)

!------------------------------------------------------------------------

!*    LOCAL VARIABLES

INTEGER :: i, ii, j, jj, k

REAL :: advfle(nz,nr+1,nc), advflw(nz,nr+1,nc)
REAL :: advflwe(nz,nr+1,nc),advflww(nz,nr+1,nc)
REAL :: advfupe(nz,nr+1,nc),advfupw(nz,nr+1,nc)
REAL :: xdiflv(nz,nr,nc+1), xdiflvj(nz,nr,nc+1)
REAL :: ydiflu(nz,nr+1,nc), ydifluj(nz,nr+1,nc)

REAL :: cfl, difcor, psi, vadvv, vwt
REAL :: xdifv, yadvu, yadvus, ydifu
!REAL :: fnlim, gx2u, gy2u, gy2v, gz2u, gz2v, u1atc, v1atc, v1atu

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ADVFLE*   REAL      TVD ADVECTIVE FLUX OF U2 IN Y-DIRECTION AT V-NODE
!                         USING U2-VALUES EAST OF V-NODE         [m2/s2]
!    *ADVFLW*   REAL      TVD ADVECTIVE FLUX OF U2 IN Y-DIRECTION AT V-NODE
!                         USING U2-VALUES WEST OF V-NODE         [m2/s2]
!    *ADVFLWE*  REAL      LW ADVECTIVE FLUX OF U2 IN Y-DIRECTION AT V-NODE
!                         USING U2-VALUES EAST OF V-NODE         [m2/s2]
!    *ADVFLWW*  REAL      LW ADVECTIVE FLUX OF U2 IN Y-DIRECTION AT V-NODE
!                         USING U2-VALUES WEST OF V-NODE         [m2/s2]
!    *ADVFUPE*  REAL      UPWIND ADVECTIVE FLUX OF U2 IN Y-DIRECTION AT V-NODE
!                         USING U2-VALUES EAST OF V-NODE         [m2/s2]
!    *ADVFUPW*  REAL      UPWIND ADVECTIVE FLUX OF U2 IN Y-DIRECTION AT V-NODE
!                         USING U2-VALUES WEST OF V-NODE         [m2/s2]
!    *XDIFLV*   REAL      PART OF (2,1)-COMPONENT OF DIFFUSIVE STRESS
!                         EVALUATED AT U-NODE                    [m2/s2]
!    *YDIFLU*   REAL      PART OF (2,1)-COMPONENT OF DIFFUSIVE STRESS
!                         EVALUATED AT V-NODE                    [m2/s2]
!    *XDIFLVJ*  REAL      PART OF (2,1)-COMPONENT OF DIFFUSIVE STRESS
!                         EVALUATED AT U-NODE TIMES J            [m2/s2]
!    *YDIFLUJ*  REAL      PART OF (2,1)-COMPONENT OF DIFFUSIVE STRESS
!                         EVALUATED AT V-NODE TIMES J            [m2/s2]
!    *CFL*      REAL      COURANT-FRIEDRICH-LEVY NUMBER
!    *PSI*      REAL      FLUX LIMITER

!------------------------------------------------------------------------

!     1. INITIALISATION
!     -----------------


CALL zero33(advfupw,nz,nr+1,nc)
CALL zero33(advflww,nz,nr+1,nc)
CALL zero33(advflw,nz,nr+1,nc)
CALL zero33(advfupe,nz,nr+1,nc)
CALL zero33(advflwe,nz,nr+1,nc)
CALL zero33(advfle,nz,nr+1,nc)
CALL zero33(ydiflu,nz,nr+1,nc)
CALL zero33(xdiflv,nz,nr,nc+1)
CALL zero33(ydifluj,nz,nr+1,nc)
CALL zero33(xdiflvj,nz,nr,nc+1)
CALL zero33(yintu,1,nr,nc+1)


!     2. ADVECTIVE FLUXES AT V-NODES
!     ------------------------------

IF (iadvc > 0) THEN
  
  
!     2.1 UPWIND AND LW-FLUXES AT INTERIOR POINTS
!     -------------------------------------------
  
!      ---using U2-values west of V-node
  DO  i=2,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        DO  k=1,nz
          vwt = SIGN(1.0,vadv(k,j,i))
          cfl = vadv(k,j,i)*del3/gy2v(j)
          advfupw(k,j,i) = 0.5*vadv(k,j,i)  &
              *((1.0+vwt)*u2phi(k,j-1,i)+(1.0-vwt)*u2phi(k,j,i))
          advflww(k,j,i) = 0.5*vadv(k,j,i)  &
              *((1.0+cfl)*u2phi(k,j-1,i)+(1.0-cfl)*u2phi(k,j,i))
        END DO
      END IF
    END DO
  END DO
  211   CONTINUE
  
!      ---using U2-values east of V-node
  DO  i=1,nc-1
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        DO  k=1,nz
          vwt = SIGN(1.0,vadv(k,j,i))
          cfl = vadv(k,j,i)*del3/gy2v(j)
          advfupe(k,j,i) = 0.5*vadv(k,j,i)  &
              *((1.0+vwt)*u2phi(k,j-1,i+1)+(1.0-vwt)*u2phi(k,j,i+1))
          advflwe(k,j,i) = 0.5*vadv(k,j,i)  &
              *((1.0+cfl)*u2phi(k,j-1,i+1)+(1.0-cfl)*u2phi(k,j,i+1))
        END DO
      END IF
    END DO
  END DO
  212   CONTINUE
  
  
!     2.2 ADVECTIVE FLUXES AT INTERIOR POINTS
!     ---------------------------------------
  
!      ---using U2-values west of V-node
  DO  i=2,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        DO  k=1,nz
          psi = 0.0
          IF (iadvc == 1) THEN
!          ---fully upwind
            psi = 0.0
          ELSE IF (iadvc == 2) THEN
!          ---Lax-Wendroff
            psi = 1.0
          ELSE IF ((iadvc == 3).OR.(iadvc == 4)) THEN
!          ---TVD scheme
            vadvv = vadv(k,j,i)
            IF (vadvv > 0.0) THEN
              psi = fnlim(advflww(k,j,i),  advfupw(k,j,i),  &
                  advflww(k,j-1,i),advfupw(k,j-1,i),iadvc)
            ELSE IF (vadvv < 0.0) THEN
              psi = fnlim(advflww(k,j,i),  advfupw(k,j,i),  &
                  advflww(k,j+1,i),advfupw(k,j+1,i),iadvc)
            END IF
          END IF
          advflw(k,j,i) = (advfupw(k,j,i) +  &
              psi*(advflww(k,j,i)-advfupw(k,j,i))) *gz2v(k,j,i)*cosphiv(j)
        END DO
      END IF
    END DO
  END DO
  221   CONTINUE
  
!      ---using U2-values east of V-node
  DO  i=1,nc-1
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        DO  k=1,nz
          psi = 0.0
          IF (iadvc == 1) THEN
!          ---fully upwind
            psi = 0.0
          ELSE IF (iadvc == 2) THEN
!          ---Lax-Wendroff
            psi = 1.0
          ELSE IF ((iadvc == 3).OR.(iadvc == 4)) THEN
!         ---TVD scheme
            vadvv = vadv(k,j,i)
            IF (vadvv > 0.0) THEN
              psi = fnlim(advflwe(k,j,i),  advfupe(k,j,i),  &
                  advflwe(k,j-1,i),advfupe(k,j-1,i),iadvc)
            ELSE IF (vadvv < 0.0) THEN
              psi = fnlim(advflwe(k,j,i),  advfupe(k,j,i),  &
                  advflwe(k,j+1,i),advfupe(k,j+1,i),iadvc)
            END IF
          END IF
          advfle(k,j,i) = (advfupe(k,j,i) +  &
              psi*(advflwe(k,j,i)-advfupe(k,j,i))) *gz2v(k,j,i)*cosphiv(j)
        END DO
      END IF
    END DO
  END DO
  222   CONTINUE
  
  
!     2.3 ADVECTIVE FLUXES AT V-OPEN BOUNDARIES
!     -----------------------------------------
  
  DO  jj=1,nobv
    i = iobv(jj)
    j = jobv(jj)
    DO  k=1,nz
      IF (soutob(jj)) THEN
        IF (vadv(k,j,i) <= 0.0) THEN
          advflw(k,j,i) = vadv(k,j,i)*u2phi(k,j,i) *gz2(k,j,i)*cosphiv(j)
          advfle(k,j,i) = vadv(k,j,i)*u2phi(k,j,i+1) *gz2(k,j,i)*cosphiv(j)
        ELSE
          advflw(k,j,i) = advflw(k,j+1,i)
          advfle(k,j,i) = advfle(k,j+1,i)
        END IF
      ELSE
        IF (vadv(k,j,i) >= 0.0) THEN
          advflw(k,j,i) = vadv(k,j,i)*u2phi(k,j-1,i) *gz2(k,j-1,i)*cosphiv(j)
          advfle(k,j,i) = vadv(k,j,i)*u2phi(k,j-1,i+1)  &
              *gz2(k,j-1,i)*cosphiv(j)
        ELSE
          advflw(k,j,i) = advflw(k,j-1,i)
          advfle(k,j,i) = advfle(k,j-1,i)
        END IF
      END IF
    END DO
  END DO
  231   CONTINUE
  
END IF


!     3. DIFFUSIVE FLUXES
!     -------------------

IF (jdifc > 0) THEN
  
  
!     3.1 INTERIOR POINTS
!     -------------------
  
!      ---(2,1)-component at V-nodes
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        DO  k=1,nz
          ydiflu(k,j,i) = heddyvv(k,j,i)  &
              *(u1atc(k,j,i) - u1atc(k,j-1,i))/gy2v(j)
          ydifluj(k,j,i) = gz2v(k,j,i)*ydiflu(k,j,i)*cosphiv(j)
        END DO
      END IF
    END DO
  END DO
  311   CONTINUE
  
!      ---(2,1)-component at U-nodes
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) THEN
        DO  k=1,nz
          xdiflv(k,j,i) = heddyvu(k,j,i)  &
              *((v1atc(k,j,i)-v1atc(k,j,i-1))/gx2u(j,i) +sphcur(j)*u1(k,j,i))
          xdiflvj(k,j,i) = gz2u(k,j,i)*xdiflv(k,j,i)*cosphi(j)
        END DO
      END IF
    END DO
  END DO
  312   CONTINUE
  
  
!     3.2 DIFFUSIVE FLUXES AT OPEN BOUNDARIES
!     ---------------------------------------
  
!      ---(2,1)-component at V-nodes
  DO  jj=1,nobv
    i = iobv(jj)
    j = jobv(jj)
    DO  k=1,nz
      IF (soutob(jj)) THEN
        ydiflu(k,j,i)  = ydiflu(k,j+1,i)
        ydifluj(k,j,i) = ydifluj(k,j+1,i)
      ELSE
        ydiflu(k,j,i)  = ydiflu(k,j-1,i)
        ydifluj(k,j,i) = ydifluj(k,j-1,i)
      END IF
    END DO
  END DO
  321   CONTINUE
  
!      ---(2,1)-component at U-nodes
  DO  ii=1,nobu
    i = iobu(ii)
    j = jobu(ii)
    DO  k=1,nz
      IF (westob(ii)) THEN
        xdiflv(k,j,i)  = xdiflv(k,j,i+1)
        xdiflvj(k,j,i) = xdiflvj(k,j,i+1)
      ELSE
        xdiflv(k,j,i)  = xdiflv(k,j,i-1)
        xdiflvj(k,j,i) = xdiflvj(k,j,i-1)
      END IF
    END DO
  END DO
  322   CONTINUE
  
END IF


!     4. ADVECTIVE AND DIFFUSIVE TERMS
!     --------------------------------

!     4.1 ADVECTION
!     -------------

IF (iadvc > 0) THEN
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) THEN
        DO  k=1,nz
          yadvu = 0.5*(advflw(k,j+1,i) + advfle(k,j+1,i-1) -  &
              advflw(k,j,i)   - advfle(k,j,i-1)) /(gy2u(j)*cosphi(j))
          yadvus = - sphcur(j)*u1(k,j,i)*v1atu(k,j,i)
          yintu(j,i) = yintu(j,i) - yadvu - yadvus*gz2u(k,j,i)
          sexp(k,j,i) = sexp(k,j,i) - del3*(yadvu/gz2u(k,j,i)+yadvus)
        END DO
      END IF
    END DO
  END DO
  411   CONTINUE
END IF


!     4.2 DIFFUSION
!     -------------

IF (jdifc > 0) THEN
  
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) THEN
        DO  k=1,nz
          ydifu = 0.5*(ydifluj(k,j+1,i) + ydifluj(k,j+1,i-1) -  &
              ydifluj(k,j,i)   - ydifluj(k,j,i-1))/gy2u(j)
          IF (j == 1) THEN
            xdifv = (xdiflvj(k,j+1,i) - xdiflvj(k,j,i))  &
                /(0.5*gy2u(j+1)+1.5*gy2u(j))
          ELSE IF (j == nr) THEN
            xdifv = (xdiflvj(k,j,i) - xdiflvj(k,j-1,i))  &
                /(0.5*gy2u(j-1)+1.5*gy2u(j))
          ELSE
            xdifv = (xdiflvj(k,j+1,i) - xdiflvj(k,j-1,i))  &
                /(0.5*(gy2u(j-1)+gy2u(j+1))+gy2u(j))
          END IF
          difcor = -sphcur(j)*(0.25*(ydiflu(k,j,i)+ydiflu(k,j+1,i)+  &
              ydiflu(k,j,i-1)+ydiflu(k,j+1,i-1))+xdiflv(k,j,i))
          sexp(k,j,i) = sexp(k,j,i) + del3*  &
              ((ydifu+xdifv)/(gz2u(k,j,i)*cosphi(j))+difcor)
          yintu(j,i) = yintu(j,i) + (ydifu+xdifv)/cosphi(j)  &
              + difcor*gz2u(k,j,i)
        END DO
      END IF
    END DO
  END DO
  421   CONTINUE
  
END IF


RETURN

END SUBROUTINE yhad3du

!========================================================================

SUBROUTINE xhad3dv (v2phi,uadv,sexp,xintv,jdifc)
!************************************************************************

!    *XHAD3DV*    EVALUATE X-DIRECTION ADVECTION-DIFFUSION TERMS IN
!                 V-MOMENTUM EQUATION INCLUDING DEPTH INTEGRALS FOR
!                 VD2-EQUATION

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 19 Apr 1999 @(COHERENS)hadvdis.f 8.4

!       DESCRIPTION - COMPUTES X-DIRECTION HORIZONTAL ADVECTION/DIFFUSION
!                     TERMS FOR V2
!                   - ADVECTIVE/DIFFUSIVE TERMS ARE EVALUATED EXPLICITLY
!                   - ADVECTIVE FLUXES ARE EVALUATED AT U-NODES USING EITHER
!                     V2-VALUES SOUTH OR NORTH OF U-NODE
!                   - ADVECTIVE TERMS ARE OBTAINED BY THE AVERAGE OF
!                     THE SOUTH/NORTH FLUXES
!                   - DIFFUSIVE FLUXES USING X-(Y-)DERIVATIVES ARE EVALUATED
!                     AT THE U-(V-)NODES
!                   - ADVECTION SCHEME SELECTED BY IADVC
!                      = 1 => FIRST ORDER UPWIND
!                      = 2 => LAX-WENDROFF
!                      = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                      = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!                   - DIFFUSION SELECTED BY JDIFC
!                      = 0 => OFF
!                      = 1 => UNIFORM
!                      = 2 => SMAGORINSKY

!       REFERENCES - Sections III-4.3.3d and 4.3.5b of the User Documentation

!       CALLING PROGRAM - VCALC

!       EXTERNALS - FNLIM, ZERO

!************************************************************************
!*    ARGUMENTS


REAL, INTENT(IN)                         :: v2phi(nz,nr+1,nc)
REAL, INTENT(IN)                         :: uadv(nz,nr,nc+1)
REAL, INTENT(IN OUT)                        :: sexp(nz,nr+1,nc)
REAL, INTENT(IN OUT)                        :: xintv(nr+1,nc)
INTEGER, INTENT(IN)                  :: jdifc

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *V2PHI*    REAL      CURRENT TO BE ADVECTED/DIFFUSED [m/s]         (IN)
!    *UADV*     REAL      ADVECTING VELOCITY [m/s]                      (IN)
!    *SEXP*     REAL      RIGHT HAND SIDE OF V2-EQUATION (EXPLICIT TERMS)
!                                                        [m/s]      (IN/OUT)
!    *XINTV*    REAL      DEPTH-INTEGRATED HORIZONTAL ADVECTION TERM FOR V2
!                         IN X-DIRECTION                 [m2/s2]       (OUT)
!    *JDIFC*    INTEGER   SWITCH TO SELECT HORIZONTAL DIFFUSION (0/1/2) (IN)

!------------------------------------------------------------------------

!*    LOCAL VARIABLES

INTEGER :: i, ii, j, jj, k

REAL :: advfln(nz,nr,nc+1), advfls(nz,nr,nc+1),uadv2(nz,nr,nc+1)
REAL :: advflwn(nz,nr,nc+1),advflws(nz,nr,nc+1)
REAL :: advfupn(nz,nr,nc+1),advfups(nz,nr,nc+1)
REAL :: xdiflu(nz,nr,nc),   xdiflvj(nz,nr,nc+1)
REAL :: ydifluj(nz,nr+1,nc)

REAL :: cfl, difcor, psi, uadvu, uwt
REAL :: xadvv, xadvvs, xdifv, ydifu
!REAL :: fnlim, gx2u, gx2v, gy2v, gz2u, gz2v, u1atc, u1atv, v1atc

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ADVFLN*   REAL      TVD ADVECTIVE FLUX OF V2 IN X-DIRECTION AT U-NODE
!                         USING V2-VALUES NORTH OF U-NODE         [m2/s2]
!    *ADVFLS*   REAL      TVD ADVECTIVE FLUX OF V2 IN X-DIRECTION AT U-NODE
!                         USING V2-VALUES SOUTH OF U-NODE         [m2/s2]
!    *ADVFLWN*  REAL      LW ADVECTIVE FLUX OF V2 IN X-DIRECTION AT U-NODE
!                         USING V2-VALUES NORTH OF U-NODE         [m2/s2]
!    *ADVFLWS*  REAL      LW ADVECTIVE FLUX OF V2 IN X-DIRECTION AT U-NODE
!                         USING V2-VALUES SOUTH OF U-NODE         [m2/s2]
!    *ADVFUPN*  REAL      UPWIND ADVECTIVE FLUX OF V2 IN X-DIRECTION AT U-NODE
!                         USING V2-VALUES NORTH OF U-NODE         [m2/s2]
!    *ADVFUPS*  REAL      UPWIND ADVECTIVE FLUX OF V2 IN X-DIRECTION AT U-NODE
!                         USING V2-VALUES SOUTH OF U-NODE         [m2/s2]
!    *XDIFLVJ*  REAL      PART OF (1,2)-COMPONENT OF DIFFUSIVE STRESS
!                         EVALUATED AT U-NODE TIMES J             [m2/s2]
!    *XDIFLU*   REAL      (1,1)-COMPONENT OF DIFFUSIVE STRESS AT CELL CENTRE
!                         (WITHOUT CURVATURE TERM)                [m2/s2]
!    *YDIFLUJ*  REAL      PART OF (1,2)-COMPONENT OF DIFFUSIVE STRESS
!                         EVALUATED AT V-NODE TIMES J             [m2/s2]
!    *CFL*      REAL      COURANT-FRIEDRICH-LEVY NUMBER
!    *PSI*      REAL      FLUX LIMITER

!------------------------------------------------------------------------

!     1. INITIALISATION
!     -----------------


CALL zero33(advfups,nz,nr,nc+1)
CALL zero33(advflws,nz,nr,nc+1)
CALL zero33(advfls,nz,nr,nc+1)
CALL zero33(advfupn,nz,nr,nc+1)
CALL zero33(advflwn,nz,nr,nc+1)
CALL zero33(advfln,nz,nr,nc+1)
CALL zero33(xdiflvj,nz,nr,nc+1)
CALL zero33(ydifluj,nz,nr+1,nc)
CALL zero33(xdiflu,nz,nr,nc)
CALL zero33(xintv,1,nr+1,nc)

DO  i=1,nc+1
DO  j=1,nr
DO  k=1,nz
  uadv2(k,j,i) = uadv(k,j,i)+ugeo
END DO
END DO
END DO

!     2. ADVECTIVE FLUXES AT U-NODES
!     ------------------------------

IF (iadvc > 0) THEN
  
  
!     2.1 UPWIND AND LW-FLUXES AT INTERIOR POINTS
!     -------------------------------------------
  
!      ---using V2-values south of U-node
  DO  i=2,nc
    DO  j=2,nr
      IF (npix(j,i) == 1) THEN
        DO  k=1,nz
          uwt = SIGN(1.0,uadv2(k,j,i))
          cfl = uadv2(k,j,i)*del3/gx2u(j,i)
          advfups(k,j,i) = 0.5*uadv2(k,j,i)  &
              *((1.0+uwt)*v2phi(k,j,i-1)+(1.0-uwt)*v2phi(k,j,i))
          advflws(k,j,i) = 0.5*uadv2(k,j,i)  &
              *((1.0+cfl)*v2phi(k,j,i-1)+(1.0-cfl)*v2phi(k,j,i))
        END DO
      END IF
    END DO
  END DO
  211   CONTINUE
  
!      ---using V2-values north of U-node
  DO  i=2,nc
    DO  j=1,nr-1
      IF (npix(j,i) == 1) THEN
        DO  k=1,nz
          uwt = SIGN(1.0,uadv2(k,j,i))
          cfl = uadv2(k,j,i)*del3/gx2u(j,i)
          advfupn(k,j,i) = 0.5*uadv2(k,j,i)  &
              *((1.0+uwt)*v2phi(k,j+1,i-1)+(1.0-uwt)*v2phi(k,j+1,i))
          advflwn(k,j,i) = 0.5*uadv2(k,j,i)  &
              *((1.0+cfl)*v2phi(k,j+1,i-1)+(1.0-cfl)*v2phi(k,j+1,i))
        END DO
      END IF
    END DO
  END DO
  212   CONTINUE
  

!     2.2 ADVECTIVE FLUXES AT INTERIOR POINTS
!     ---------------------------------------
  
!      ---using V2-values south of U-node
  DO  i=2,nc
    DO  j=2,nr
      IF (npix(j,i) == 1) THEN
        DO  k=1,nz
          psi = 0.0
          IF (iadvc == 1) THEN
!          ---fully upwind
            psi = 0.0
          ELSE IF (iadvc == 2) THEN
!          ---Lax-Wendroff
            psi = 1.0
          ELSE IF ((iadvc == 3).OR.(iadvc == 4)) THEN
!          ---TVD scheme
            uadvu = uadv2(k,j,i)
            IF (uadvu > 0.0) THEN
              psi = fnlim(advflws(k,j,i),  advfups(k,j,i),  &
                  advflws(k,j,i-1),advfups(k,j,i-1),iadvc)
            ELSE IF (uadvu < 0.0) THEN
              psi = fnlim(advflws(k,j,i),  advfups(k,j,i),  &
                  advflws(k,j,i+1),advfups(k,j,i+1),iadvc)
            END IF
          END IF
          advfls(k,j,i) = (advfups(k,j,i) +  &
              psi*(advflws(k,j,i)-advfups(k,j,i))) *gz2u(k,j,i)
        END DO
      END IF
    END DO
  END DO
  221   CONTINUE
  
!      ---using V2-values north of U-node
  DO  i=2,nc
    DO  j=1,nr-1
      IF (npix(j,i) == 1) THEN
        DO  k=1,nz
          psi = 0.0
          IF (iadvc == 1) THEN
!          ---fully upwind
            psi = 0.0
          ELSE IF (iadvc == 2) THEN
!          ---Lax-Wendroff
            psi = 1.0
          ELSE IF ((iadvc == 3).OR.(iadvc == 4)) THEN
!          ---TVD scheme
            uadvu = uadv2(k,j,i)
            IF (uadvu > 0.0) THEN
              psi = fnlim(advflwn(k,j,i),  advfupn(k,j,i),  &
                  advflwn(k,j,i-1),advfupn(k,j,i-1),iadvc)
            ELSE IF (uadvu < 0.0) THEN
              psi = fnlim(advflwn(k,j,i),  advfupn(k,j,i),  &
                  advflwn(k,j,i+1),advfupn(k,j,i+1),iadvc)
            END IF
          END IF
          advfln(k,j,i) = (advfupn(k,j,i) +  &
              psi*(advflwn(k,j,i)-advfupn(k,j,i))) *gz2u(k,j,i)
        END DO
      END IF
    END DO
  END DO
  222   CONTINUE


!     2.3 ADVECTIVE FLUXES AT U-OPEN BOUNDARIES
!     -----------------------------------------
  
  DO  ii=1,nobu
    i = iobu(ii)
    j = jobu(ii)

    DO  k=1,nz
      IF (westob(ii)) THEN
        IF (uadv2(k,j,i) <= 0.0) THEN
          advfls(k,j,i) = uadv2(k,j,i)*v2phi(k,j,i)*gz2(k,j,i)
          advfln(k,j,i) = uadv2(k,j,i)*v2phi(k,j+1,i)*gz2(k,j,i)
        ELSE
          advfls(k,j,i) = advfls(k,j,i+1)
          advfln(k,j,i) = advfln(k,j,i+1)
        END IF
      ELSE
        IF (uadv2(k,j,i) >= 0.0) THEN
          advfls(k,j,i) = uadv2(k,j,i)*v2phi(k,j,i-1)*gz2(k,j,i-1)
          advfln(k,j,i) = uadv2(k,j,i)*v2phi(k,j+1,i-1)*gz2(k,j,i-1)
        ELSE
          advfls(k,j,i) = advfls(k,j,i-1)
          advfln(k,j,i) = advfln(k,j,i-1)
        END IF
      END IF
    END DO
  END DO
  231   CONTINUE
  
END IF

!     3. DIFFUSIVE FLUXES
!     -------------------

IF (jdifc > 0) THEN
  
  
!     3.1 INTERIOR POINTS
!     -------------------
  
!      ---(1,2)-component at U-nodes
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) THEN
        DO  k=1,nz
          xdiflvj(k,j,i) = gz2u(k,j,i)*heddyvu(k,j,i)  &
              *((v1atc(k,j,i) - v1atc(k,j,i-1))/gx2u(j,i) +sphcur(j)*u1(k,j,i))
        END DO
      END IF
    END DO
  END DO
  311   CONTINUE
  
!      ---(1,2)-component at V-nodes
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        DO  k=1,nz
          ydifluj(k,j,i) = gz2v(k,j,i)*heddyvv(k,j,i)  &
              *(u1atc(k,j,i)-u1atc(k,j-1,i))/gy2v(j)
        END DO
      END IF
    END DO
  END DO
  312   CONTINUE
  
!      ---(1,1)-component at centres
  IF (igtrh == 1) THEN
    DO  i=1,nc
      DO  j=1,nr
        IF (nwd(j,i) == 1) THEN
          DO  k=1,nz
            xdiflu(k,j,i) = 2.0*heddyvc(k,j,i)  &
                *(u1(k,j,i+1)-u1(k,j,i))/gx2(j,i)
          END DO
        END IF
      END DO
    END DO
    313   CONTINUE
  END IF
  
!     3.2 DIFFUSIVE FLUXES AT OPEN BOUNDARIES
!     ---------------------------------------
  
!      ---(1,2) component at U-nodes
  DO  ii=1,nobu
    i = iobu(ii)
    j = jobu(ii)
    DO  k=1,nz
      IF (westob(ii)) THEN
        xdiflvj(k,j,i) = xdiflvj(k,j,i+1)
      ELSE
        xdiflvj(k,j,i) = xdiflvj(k,j,i-1)
      END IF
    END DO
  END DO
  321   CONTINUE
  
!      ---(1,2)-component at V-nodes
  DO  jj=1,nobv
    i = iobv(jj)
    j = jobv(jj)
    DO  k=1,nz
      IF (soutob(jj)) THEN
        ydifluj(k,j,i) = ydifluj(k,j+1,i)
      ELSE
        ydifluj(k,j,i) = ydifluj(k,j-1,i)
      END IF
    END DO
  END DO
  322   CONTINUE
  
END IF

!     4. ADVECTIVE AND DIFFUSIVE TERMS
!     --------------------------------


!     4.1 ADVECTION
!     -------------


IF (iadvc > 0) THEN
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        DO  k=1,nz
          xadvv = 0.5*(advfls(k,j,i+1) + advfln(k,j-1,i+1) -  &
              advfls(k,j,i)   - advfln(k,j-1,i))/gx2v(j,i)
          xadvvs =  sphcurv(j)*u1atv(k,j,i)**2
          xintv(j,i) = xintv(j,i) - xadvv - xadvvs*gz2v(k,j,i)
          sexp(k,j,i) = sexp(k,j,i) - del3*(xadvv/gz2v(k,j,i)+xadvvs)
        END DO
      END IF
    END DO
  END DO
  411   CONTINUE
END IF

!     4.2 DIFFUSION
!     -------------

IF (jdifc > 0) THEN
  
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        DO  k=1,nz
          xdifv = 0.5*(xdiflvj(k,j,i+1) + xdiflvj(k,j-1,i+1) -  &
              xdiflvj(k,j,i)   - xdiflvj(k,j-1,i))/gx2v(j,i)
          IF (i == 1) THEN
            ydifu = (ydifluj(k,j,i+1) - ydifluj(k,j,i))  &
                /(0.5*gx2v(j,i+1)+1.5*gx2v(j,i))
          ELSE IF (i == nc) THEN
            ydifu = (ydifluj(k,j,i) - ydifluj(k,j,i-1))  &
                / (0.5*gx2v(j,i-1)+1.5*gx2v(j,i))
          ELSE
            ydifu = (ydifluj(k,j,i+1) - ydifluj(k,j,i-1))  &
                /(0.5*(gx2v(j,i-1)+gx2v(j,i+1))+gx2v(j,i))
          END IF
          difcor = sphcurv(j)*(0.5*(xdiflu(k,j-1,i)+xdiflu(k,j,i))  &
              - 2.0*sphcurv(j)*heddyvv(k,j,i)*v1(k,j,i))
          sexp(k,j,i) = sexp(k,j,i) + del3*((xdifv+ydifu)/gz2v(k,j,i)+difcor)
          xintv(j,i) = xintv(j,i) + xdifv + ydifu + difcor*gz2v(k,j,i)
        END DO
      END IF
    END DO
  END DO
  421   CONTINUE
  
END IF

RETURN

END SUBROUTINE xhad3dv

!========================================================================

SUBROUTINE yhad3dv (v2phi,vadv,sexp,yintv,jdifc)
!************************************************************************

!    *YHAD3DV*    EVALUATE Y-DIRECTION ADVECTION-DIFFUSION TERMS IN
!                 V-MOMENTUM EQUATION INCLUDING DEPTH INTEGRALS FOR
!                 VD2-EQUATION

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 26 Mar 1999 @(COHERENS)hadvdis.f 8.4

!       DESCRIPTION - COMPUTES Y-DIRECTION HORIZONTAL ADVECTION/DIFFUSION
!                     TERMS FOR V2
!                   - ADVECTIVE/DIFFUSIVE TERMS ARE EVALUATED EXPLICITLY
!                   - ADVECTION SCHEME SELECTED BY IADVC
!                      = 1 => FIRST ORDER UPWIND
!                      = 2 => LAX-WENDROFF
!                      = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                      = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!                   - DIFFUSION SELECTED BY JDIFC
!                      = 0 => OFF
!                      = 1 => UNIFORM
!                      = 2 => SMAGORINSKY

!       REFERENCES - Sections III-4.3.3b and 4.3.5b of the User Documentation

!       CALLING PROGRAM - VCALC

!       EXTERNALS - FNLIM, ZERO

!************************************************************************
!*    ARGUMENTS


REAL, INTENT(IN)                     :: v2phi(nz,nr+1,nc)
REAL, INTENT(IN)                         :: vadv(nz,nr+1,nc)
REAL, INTENT(IN OUT)                        :: sexp(nz,nr+1,nc)
REAL, INTENT(IN OUT)                        :: yintv(nr+1,nc)
INTEGER, INTENT(IN)                  :: jdifc




!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *V2PHI*    REAL      CURRENT TO BE ADVECTED/DIFFUSED [m/s]      (IN)
!    *VADV*     REAL      ADVECTING VELOCITY              [m/s]      (IN)
!    *SEXP*     REAL      RIGHT HAND SIDE OF V2-EQUATION (EXPLICIT TERMS)
!                                                         [m/s]  (IN/OUT)
!    *YINTV*    REAL      DEPTH INTEGRATED HORIZONTAL ADVECTION TERM FOR V2
!                         IN Y-DIRECTION                  [m2/s2]   (OUT)
!    *JDIFC*    INTEGER   SWITCH TO SELECT HORIZONTAL DIFFUSION (0/1/2) (IN)

!------------------------------------------------------------------------

!*    LOCAL VARIABLES

INTEGER :: i, j, k
REAL :: advfl(nz,nr,nc), advflw(nz,0:nr+1,nc)
REAL :: advfup(nz,0:nr+1,nc), ydiflvj(nz,nr,nc)
REAL :: cfl, psi, vadvc, vwt, yadvv, ydifv
!REAL :: fnlim, gy2v, gz2v

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ADVFL*    REAL      TVD ADVECTIVE FLUX OF V2 IN Y-DIRECTION AT CELL
!                         CENTRE                                       [m2/s2]
!    *ADVFLW*   REAL      LW  FLUX OF V2 IN Y-DIRECTION AT CELL CENTRE [m2/s2]
!    *ADVFUP*   REAL      UPWIND FLUX OF V2 IN Y-DIRECTION AT CELL CENTRE
!                                                                      [m2/s2]
!    *YDIFLVJ*  REAL      (2,2)-COMPONENT OF DIFFUSIVE STRESS AT CELL CENTRE
!                         TIMES J                                      [m2/s2]
!    *CFL*      REAL      COURANT-FRIEDRICH-LEVY NUMBER
!    *PSI*      REAL      FLUX LIMITER

!------------------------------------------------------------------------

!     1. INITIALISATION OF TEMPORARY ARRAYS
!     -------------------------------------


CALL zero33(advfl,nz,nr,nc)
CALL zero33(advflw,nc*(nr+2)*nz,1,1)
CALL zero33(advfup,nc*(nr+2)*nz,1,1)
CALL zero33(ydiflvj,nz,nr,nc)
CALL zero33(yintv,1,nr+1,nc)


!     2. ADVECTIVE FLUXES AT CENTRES
!     ------------------------------

IF (iadvc > 0) THEN
  
  
!     2.1 UPWIND AND LW FLUXES
!     ------------------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
          vadvc = 0.5*(vadv(k,j,i)+vadv(k,j+1,i))
          vwt = SIGN(1.0,vadvc)
          cfl = vadvc*del3/gy2(j)
          advfup(k,j,i) = 0.5*vadvc  &
              *((1.0+vwt)*v2phi(k,j,i)+(1.0-vwt)*v2phi(k,j+1,i))
          advflw(k,j,i) = 0.5*vadvc  &
              *((1.0+cfl)*v2phi(k,j,i)+(1.0-cfl)*v2phi(k,j+1,i))
        END DO
      END IF
    END DO
  END DO
  211   CONTINUE
  
  
!     2.2 ADVECTIVE FLUXES
!     --------------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
          psi = 0.0
          IF (iadvc == 1) THEN
!          ---fully upwind
            psi = 0.0
          ELSE IF (iadvc == 2) THEN
!          ---Lax-Wendroff
            psi = 1.0
          ELSE IF ((iadvc == 3).OR.(iadvc == 4)) THEN
!          ---TVD scheme
            vadvc = 0.5*(vadv(k,j,i)+vadv(k,j+1,i))
            IF (vadvc > 0.0) THEN
              psi = fnlim(advflw(k,j,i),  advfup(k,j,i),  &
                  advflw(k,j-1,i),advfup(k,j-1,i),iadvc)
            ELSE IF (vadvc < 0.0) THEN
              psi = fnlim(advflw(k,j,i),  advfup(k,j,i),  &
                  advflw(k,j+1,i),advfup(k,j+1,i),iadvc)
            END IF
          END IF
          advfl(k,j,i) = (advfup(k,j,i) +  &
              psi*(advflw(k,j,i)-advfup(k,j,i))) *gz2(k,j,i)*cosphi(j)
        END DO
      END IF
    END DO
  END DO
  221   CONTINUE
  
END IF


!     3. DIFFUSIVE FLUXES
!     -------------------

IF (jdifc > 0) THEN
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
          ydiflvj(k,j,i) = 2.0*heddyvc(k,j,i)  &
              *(v1(k,j+1,i)-v1(k,j,i))/gy2(j) *gz2(k,j,i)*cosphi(j)
        END DO
      END IF
    END DO
  END DO
  301   CONTINUE
  
END IF


!     4. ADVECTIVE AND DIFFUSIVE TERMS AT V-NODES
!     -------------------------------------------

IF ((iadvc > 0).OR.(jdifc > 0)) THEN
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        DO  k=1,nz
          yadvv = (advfl(k,j,i)-advfl(k,j-1,i))/(gy2v(j)*cosphiv(j))
          ydifv = (ydiflvj(k,j,i)-ydiflvj(k,j-1,i)) /(gy2v(j)*cosphiv(j))
          yintv(j,i) = yintv(j,i) - yadvv + ydifv
          sexp(k,j,i) = sexp(k,j,i) - del3*(yadvv-ydifv) /gz2v(k,j,i)
        END DO
      END IF
    END DO
  END DO
  401   CONTINUE
END IF

RETURN

END SUBROUTINE yhad3dv

!========================================================================

SUBROUTINE had2du
!************************************************************************

!    *HAD2DU*    2-D HORIZONTAL ADVECTION AND DIFFUSION IN THE UD2-MOMENTUM
!                EQUATION

!       AUTHOR - PATRICK LUYTEN AND KEVIN RUDDICK

!       LAST UPDATE - 26 Mar 1999 @(COHERENS)hadvdis.f 8.4

!       DESCRIPTION - EVALUATES HORIZONTAL ADVECTION/DIFFUSION
!                     TERMS FOR UD2
!                   - ADVECTIVE/DIFFUSIVE TERMS ARE EVALUATED EXPLICITLY
!                   - PROCEDURES ARE SIMILAR TO THE ONES USED IN
!                     XHAD3DU, YHAD3DU
!                   - ADVECTION SCHEME SELECTED BY IADVC
!                      = 1 => FIRST ORDER UPWIND
!                      = 2 => LAX-WENDROFF
!                      = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                      = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!                   - DIFFUSION SELECTED BY IODIF
!                      = 0 => OFF
!                      = 1 => UNIFORM
!                      = 2 => SMAGORINSKYI

!       REFERENCES - Sections III-4.3.3,4.3.4,4.3.5,4.3.6 and V-1.8 of the
!                    User Documentation

!       CALLING PROGRAM - CRRNT2, UCALC

!       EXTERNALS - FNLIM, ZERO

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, ii, j, jj

REAL :: advflc(nr,nc),     advfle(nr+1,nc), advflw(nr+1,nc)
REAL :: advflwc(nr,0:nc+1),advflwe(nr+1,nc),advflww(nr+1,nc)
REAL :: advfupc(nr,0:nc+1),advfupe(nr+1,nc),advfupw(nr+1,nc)
REAL :: xdiflu(nr,nc),     xdiflv(nr,nc+1), ydiflu(nr+1,nc)

REAL :: cfl, psi, udadv, uwt, vdadv, vwt
REAL :: xadvu, xdifu, xdifv, yadvu, ydifu
!REAL :: fnlim, gx2u, gy2u, gy2v, h2atc, h2atu, h2atv
!REAL :: ud2atc, vd2atc, vd2atu

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ADVFLC*   REAL      TVD ADVECTIVE FLUX OF UD2 IN X-DIRECTION AT
!                         CELL CENTRE                          [m3/s2]
!    *ADVFLE*   REAL      TVD ADVECTIVE FLUX IN Y-DIRECTION AT V-NODE
!                         USING UD2-VALUES EAST OF V-NODE      [m3/s2]
!    *ADVFLW*   REAL      TVD ADVECTIVE FLUX IN Y-DIRECTION AT V-NODE
!                         USING UD2-VALUES WEST OF V-NODE      [m3/s2]
!    *ADVFLWC*  REAL      LW ADVECTIVE FLUX OF UD2 IN X-DIRECTION AT
!                         CELL CENTRE                          [m3/s2]
!    *ADVFLWE*  REAL      LW ADVECTIVE FLUX IN Y-DIRECTION AT V-NODE
!                         USING UD2-VALUES EAST OF V-NODE      [m3/s2]
!    *ADVFLWW*  REAL      LW ADVECTIVE FLUX IN Y-DIRECTION AT V-NODE
!                         USING UD2-VALUES WEST OF V-NODE      [m3/s2]
!    *ADVFUPC* REAL       UPWIND ADVECTIVE FLUX OF UD2 IN X-DIRECTION AT
!                         CELL CENTRE                          [m3/s2]
!    *ADVFUPE*  REAL      UPWIND ADVECTIVE FLUX IN Y-DIRECTION AT V-NODE
!                         USING UD2-VALUES EAST OF V-NODE      [m3/s2]
!    *ADVFUPW*  REAL      UPWIND ADVECTIVE FLUX IN Y-DIRECTION AT V-NODE
!                         USING UD2-VALUES WEST OF V-NODE      [m3/s2]
!    *XDIFLU*   REAL      (1,1)-COMPONENT OF DIFFUSIVE STRESS AT CELL CENTRE
!                                                              [m3/s2]
!    *XDIFLV*   REAL      PART OF (2,1)-COMPONENT OF DIFFUSIVE STRESS
!                         EVALUATED AT U-NODE                  [m3/s2]
!    *YDIFLU*   REAL      PART OF (2,1)-COMPONENT OF DIFFUSIVE STRESS
!                         EVALUATED AT V-NODE                  [m3/s2]
!    *CFL*      REAL      COURANT-FRIEDRICH-LEVY NUMBER
!    *PSI*      REAL      FLUX LIMITER

!------------------------------------------------------------------------

!     1. INITIALISATION OF TEMPORARY ARRAYS
!     -------------------------------------


CALL zero33(advfupc,(nc+2)*nr,1,1)
CALL zero33(advflwc,(nc+2)*nr,1,1)
CALL zero33(advflc,1,nr,nc)

CALL zero33(advfupw,1,nr+1,nc)
CALL zero33(advflww,1,nr+1,nc)
CALL zero33(advflw,1,nr+1,nc)

CALL zero33(advfupe,1,nr+1,nc)
CALL zero33(advflwe,1,nr+1,nc)
CALL zero33(advfle,1,nr+1,nc)

CALL zero33(xdiflu,1,nr,nc)
CALL zero33(ydiflu,1,nr+1,nc)
CALL zero33(xdiflv,1,nr,nc+1)


!     2. ADVECTIVE FLUXES AT CENTRES
!     ------------------------------

IF (iadvc > 0) THEN
  
  
!     2.1 UPWIND AND LW FLUXES
!     ------------------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        udadv = ud2atc(j,i)
        uwt = SIGN(1.0,udadv)
        cfl = udadv*delt/(gx2(j,i)*h2atc(j,i))
        advfupc(j,i) = 0.5*udadv *((1.0+uwt)*ud2(j,i)+(1.0-uwt)*ud2(j,i+1))
        advflwc(j,i) = 0.5*udadv *((1.0+cfl)*ud2(j,i)+(1.0-cfl)*ud2(j,i+1))
      END IF
    END DO
  END DO
  
  
!     2.2 ADVECTIVE FLUXES
!     --------------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        psi = 0.0
        IF (iadvc == 1) THEN
!         ---fully upwind
          psi = 0.0
        ELSE IF (iadvc == 2) THEN
!         ---Lax-Wendroff
          psi = 1.0
        ELSE IF ((iadvc == 3).OR.(iadvc == 4)) THEN
!         ---TVD scheme
          udadv = ud2atc(j,i)
          IF (udadv > 0.0) THEN
            psi = fnlim(advflwc(j,i),  advfupc(j,i),  &
                advflwc(j,i-1),advfupc(j,i-1),iadvc)
          ELSE IF (udadv < 0.0) THEN
            psi = fnlim(advflwc(j,i),  advfupc(j,i),  &
                advflwc(j,i+1),advfupc(j,i+1),iadvc)
          END IF
        END IF
        advflc(j,i) = (advfupc(j,i) +  &
            psi*(advflwc(j,i)-advfupc(j,i)))/h2atc(j,i)
      END IF
    END DO
  END DO
  
END IF


!     3. ADVECTIVE FLUXES AT V-NODES
!     ------------------------------

IF (iadvc > 0) THEN
  
  
!     3.1 UPWIND AND LW-FLUXES AT INTERIOR POINTS
!     -------------------------------------------
  
!      ---using UD2-values west of V-node
  DO  i=2,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        vwt = SIGN(1.0,vd2(j,i))
        cfl = vd2(j,i)*delt/(gy2v(j)*h2atv(j,i))
        advfupw(j,i) = 0.5*vd2(j,i) *((1.0+vwt)*ud2(j-1,i)+(1.0-vwt)*ud2(j,i))
        advflww(j,i) = 0.5*vd2(j,i) *((1.0+cfl)*ud2(j-1,i)+(1.0-cfl)*ud2(j,i))
      END IF
    END DO
  END DO
  
!      ---using UD2-values east of V-node
  DO  i=1,nc-1
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        vwt = SIGN(1.0,vd2(j,i))
        cfl = vd2(j,i)*delt/(gy2v(j)*h2atv(j,i))
        advfupe(j,i) = 0.5*vd2(j,i)  &
            *((1.0+vwt)*ud2(j-1,i+1)+(1.0-vwt)*ud2(j,i+1))
        advflwe(j,i) = 0.5*vd2(j,i)  &
            *((1.0+cfl)*ud2(j-1,i+1)+(1.0-cfl)*ud2(j,i+1))
      END IF
    END DO
  END DO
  
  
!     3.2 ADVECTIVE FLUXES AT INTERIOR POINTS
!     ---------------------------------------
  
!      ---using UD2-values west of V-node
  DO  i=2,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        psi = 0.0
        IF (iadvc == 1) THEN
!         ---fully upwind
          psi = 0.0
        ELSE IF (iadvc == 2) THEN
!         ---Lax-Wendroff
          psi = 1.0
        ELSE IF ((iadvc == 3).OR.(iadvc == 4)) THEN
!         ---TVD-scheme
          vdadv = vd2(j,i)
          IF (vdadv > 0.0) THEN
            psi = fnlim(advflww(j,i),  advfupw(j,i),  &
                advflww(j-1,i),advfupw(j-1,i),iadvc)
          ELSE IF (vdadv < 0.0) THEN
            psi = fnlim(advflww(j,i),  advfupw(j,i),  &
                advflww(j+1,i),advfupw(j+1,i),iadvc)
          END IF
        END IF
        advflw(j,i) = (advfupw(j,i) + psi*(advflww(j,i)-advfupw(j,i)))  &
            *cosphiv(j)/h2atv(j,i)
      END IF
    END DO
  END DO
  
!      ---using UD2-values east of V-node
  DO  i=1,nc-1
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        psi = 0.0
        IF (iadvc == 1) THEN
!         ---fully upwind
          psi = 0.0
        ELSE IF (iadvc == 2) THEN
!         ---Lax-Wendroff
          psi = 1.0
        ELSE IF ((iadvc == 3).OR.(iadvc == 4)) THEN
!         ---TVD-scheme
          vdadv = vd2(j,i)
          IF (vdadv > 0.0) THEN
            psi = fnlim(advflwe(j,i),  advfupe(j,i),  &
                advflwe(j-1,i),advfupe(j-1,i),iadvc)
          ELSE IF (vdadv < 0.0) THEN
            psi = fnlim(advflwe(j,i),  advfupe(j,i),  &
                advflwe(j+1,i),advfupe(j+1,i),iadvc)
          END IF
        END IF
        advfle(j,i) = (advfupe(j,i) + psi*(advflwe(j,i)-advfupe(j,i)))  &
            *cosphiv(j)/h2atv(j,i)
      END IF
    END DO
  END DO
  
  
!     3.3 ADVECTIVE FLUXES AT V-OPEN BOUNDARIES
!     -----------------------------------------
  
  DO  jj=1,nobv
    i = iobv(jj)
    j = jobv(jj)
    IF (soutob(jj)) THEN
      IF (vd2(j,i) <= 0.0) THEN
        advflw(j,i) = vd2(j,i)*ud2(j,i)*cosphiv(j)/h2atc(j,i)
        advfle(j,i) = vd2(j,i)*ud2(j,i+1)*cosphiv(j)/h2atc(j,i)
      ELSE
        advflw(j,i) = advflw(j+1,i)
        advfle(j,i) = advfle(j+1,i)
      END IF
    ELSE
      IF (vd2(j,i) >= 0.0) THEN
        advflw(j,i) = vd2(j,i)*ud2(j-1,i)*cosphiv(j)/h2atc(j-1,i)
        advfle(j,i) = vd2(j,i)*ud2(j-1,i+1)*cosphiv(j)/h2atc(j-1,i)
      ELSE
        advflw(j,i) = advflw(j-1,i)
        advfle(j,i) = advfle(j-1,i)
      END IF
    END IF
  END DO
  
END IF


!     4. DIFFUSIVE FLUXES
!     -------------------

IF (iodif > 0) THEN
  
  
!     4.1 INTERIOR POINTS
!     -------------------
  
!      ---(1,1)-component at centres
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        IF (npix(j,i) == 1) THEN
          xdiflu(j,i) = xdiflu(j,i) - ud2(j,i)/h2atu(j,i)
        ELSE IF (npix(j,i) > 1) THEN
          xdiflu(j,i) = xdiflu(j,i) - ud2(j,i)/h2atc(j,i)
        END IF
        IF (npix(j,i+1) == 1) THEN
          xdiflu(j,i) = xdiflu(j,i) + ud2(j,i+1)/h2atu(j,i+1)
        ELSE IF (npix(j,i+1) > 1) THEN
          xdiflu(j,i) = xdiflu(j,i) + ud2(j,i+1)/h2atc(j,i)
        END IF
        xdiflu(j,i) = 2.0*xdiflu(j,i)*dheddyvc(j,i)/gx2(j,i)  &
            - 2.0*dheddyvc(j,i)*sphcur(j)*vd2atc(j,i)/h2atc(j,i)
      END IF
    END DO
  END DO
  
!      ---(2,1)-component at V-nodes
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        ydiflu(j,i) = dheddyvv(j,i)*cosphiv(j)  &
            *(ud2atc(j,i)/h2atc(j,i)-ud2atc(j-1,i)/h2atc(j-1,i)) /gy2v(j)
      END IF
    END DO
  END DO
  
!      ---(2,1)-component at U-nodes
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) THEN
        xdiflv(j,i) = dheddyvu(j,i)*cosphi(j)  &
            *((vd2atc(j,i)/h2atc(j,i)-vd2atc(j,i-1)/h2atc(j,i-1))  &
            /gx2u(j,i) + sphcur(j)*ud2(j,i)/h2atu(j,i))
      END IF
    END DO
  END DO
  
  
!     4.2 DIFFUSIVE FLUXES AT OPEN BOUNDARIES
!     ---------------------------------------
  
!      ---(2,1)-component at V-nodes
  DO  jj=1,nobv
    i = iobv(jj)
    j = jobv(jj)
    IF (soutob(jj)) THEN
      ydiflu(j,i) = ydiflu(j+1,i)
    ELSE
      ydiflu(j,i) = ydiflu(j-1,i)
    END IF
  END DO
  
!      ---(2,1)-component at U-nodes
  DO  ii=1,nobu
    i = iobu(ii)
    j = jobu(ii)
    IF (westob(ii)) THEN
      xdiflv(j,i) = xdiflv(j,i+1)
    ELSE
      xdiflv(j,i) = xdiflv(j,i-1)
    END IF
  END DO
  
END IF


!     5. ADVECTIVE AND DIFFUSIVE TERMS
!     --------------------------------

!     5.1 ADVECTION
!     -------------

IF (iadvc > 0) THEN
  
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) THEN
        xadvu = (advflc(j,i)-advflc(j,i-1))/gx2u(j,i)
        yadvu = 0.5*(advflw(j+1,i)+ advfle(j+1,i-1) -  &
            advflw(j,i)  - advfle(j,i-1))/(gy2u(j)*cosphi(j))  &
            - sphcur(j)*ud2(j,i)*vd2atu(j,i)/h2atu(j,i)
        uah2d(j,i) = xadvu + yadvu
      END IF
    END DO
  END DO
  
END IF


!     5.2 DIFFUSION
!     -------------

IF (iodif > 0) THEN
  
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) THEN
        xdifu = (xdiflu(j,i)-xdiflu(j,i-1))/gx2u(j,i)
        ydifu = 0.5*(ydiflu(j+1,i)+ ydiflu(j+1,i-1) -  &
            ydiflu(j,i)  - ydiflu(j,i-1))/gy2u(j)
        IF (j == 1) THEN
          xdifv = (xdiflv(j+1,i) - xdiflv(j,i)) /(0.5*gy2u(j+1)+1.5*gy2u(j))
        ELSE IF (j == nr) THEN
          xdifv = (xdiflv(j,i) - xdiflv(j-1,i)) /(0.5*gy2u(j-1)+1.5*gy2u(j))
        ELSE
          xdifv = (xdiflv(j+1,i) - xdiflv(j-1,i))  &
              /(0.5*(gy2u(j-1)+gy2u(j+1))+gy2u(j))
        END IF
        udh2d(j,i) = xdifu + (ydifu+xdifv)/cosphi(j)
        udh2d(j,i) = udh2d(j,i) - sphcur(j)*  &
            (0.25*(ydiflu(j,i)+ydiflu(j+1,i)+ydiflu(j,i-1)+  &
            ydiflu(j+1,i-1))+xdiflv(j,i))
      END IF
    END DO
  END DO
  
END IF


RETURN

END SUBROUTINE had2du

!========================================================================

SUBROUTINE had2dv
!************************************************************************

!    *HAD2DV*    2-D HORIZONTAL ADVECTION AND DIFFUSION IN THE VD2-MOMENTUM
!                EQUATION

!       AUTHOR - PATRICK LUYTEN AND KEVIN RUDDICK

!       LAST UPDATE - 26 Mar 1999 @(COHERENS)hadvdis.f 8.4

!       DESCRIPTION - EVALUATES HORIZONTAL ADVECTION/DIFFUSION
!                     TERMS FOR VD2
!                   - ADVECTIVE/DIFFUSIVE TERMS ARE EVALUATED EXPLICITLY
!                   - PROCEDURES ARE SIMILAR TO THE ONES USED IN
!                     XHAD3DV, YHAD3DV
!                   - ADVECTION SCHEME SELECTED BY IADVC
!                      = 1 => FIRST ORDER UPWIND
!                      = 2 => LAX-WENDROFF
!                      = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                      = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!                   - DIFFUSION SELECTED BY IODIF
!                      = 0 => OFF
!                      = 1 => UNIFORM
!                      = 2 => SMAGORINSKYI

!       REFERENCES - Sections III-4.3.3,4.3.4,4.3.5,4.3.6 and V-1.8 of the
!                    User Documentation

!       CALLING PROGRAM - CRRNT2, VCALC

!       EXTERNALS - FNLIM, ZERO

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, ii, j, jj

REAL :: advflc(nr,nc),     advfln(nr,nc+1), advfls(nr,nc+1)
REAL :: advflwc(0:nr+1,nc),advflwn(nr,nc+1),advflws(nr,nc+1)
REAL :: advfupc(0:nr+1,nc),advfupn(nr,nc+1),advfups(nr,nc+1)
REAL :: xdiflu(nr,nc),     xdiflv(nr,nc+1)
REAL :: ydiflu(nr+1,nc),   ydiflv(nr,nc)

REAL :: cfl, psi, udadv, uwt, vdadv, vwt
REAL :: xadvv, xdifv, yadvv, ydifu, ydifv
!REAL :: fnlim, gx2u, gx2v, gy2v, h2atc, h2atu, h2atv
!REAL :: ud2atc, ud2at v, vd2atc

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ADVFLC*   REAL      TVD ADVECTIVE FLUX OF VD2 IN Y-DIRECTION AT
!                         CELL CENTRE                          [m3/s2]
!    *ADVFLN*   REAL      TVD ADVECTIVE FLUX IN X-DIRECTION AT U-NODE
!                         USING VD2-VALUES NORTH OF U-NODE     [m3/s2]
!    *ADVFLS*   REAL      TVD ADVECTIVE FLUX IN X-DIRECTION AT U-NODE
!                         USING VD2-VALUES SOUTH OF U-NODE     [m3/s2]
!    *ADVFLWC*  REAL      LW ADVECTIVE FLUX OF VD2 IN Y-DIRECTION AT
!                         CELL CENTRE                          [m3/s2]
!    *ADVFLWN*  REAL      LW ADVECTIVE FLUX IN X-DIRECTION AT U-NODE
!                         USING VD2-VALUES NORTH OF U-NODE     [m3/s2]
!    *ADVFLWS*  REAL      LW ADVECTIVE FLUX IN X-DIRECTION AT U-NODE
!                         USING VD2-VALUES SOUTH OF U-NODE     [m3/s2]
!    *ADVFUPC*  REAL      UPWIND ADVECTIVE FLUX OF VD2 IN Y-DIRECTION AT
!                         CELL CENTRE                          [m3/s2]
!    *ADVFUPN*  REAL      UPWIND ADVECTIVE FLUX IN X-DIRECTION AT U-NODE
!                         USING VD2-VALUES NORTH OF U-NODE     [m3/s2]
!    *ADVFUPS*  REAL      UPWIND ADVECTIVE FLUX IN X-DIRECTION AT U-NODE
!                         USING VD2-VALUES SOUTH OF U-NODE     [m3/s2]
!    *XDIFLU*   REAL      (1,1)-COMPONENT OF DIFFUSIVE STRESS AT CELL CENTRE
!                         (WITHOUT CURVATURE TERM)             [m3/s2]
!    *XDIFLV*   REAL      PART OF (1,2)-COMPONENT OF DIFFUSIVE STRESS
!                         EVALUATED AT U-NODE                  [m3/s2]
!    *YDIFLU*   REAL      PART OF (1,2)-COMPONENT OF DIFFUSIVE STRESS
!                         EVALUATED AT V-NODE                  [m3/s2]
!    *YDIFLV*   REAL      (2,2)-COMPONENT OF DIFFUSIVE STRESS AT CELL CENTRE
!                                                              [m3/s2]
!    *CFL*      REAL      COURANT-FRIEDRICH-LEVY NUMBER
!    *PSI*      REAL      FLUX LIMITER

!------------------------------------------------------------------------

!     1. INITIALISATION OF TEMPORARY ARRAYS
!     -------------------------------------


CALL zero33(advfupc,nc*(nr+2),1,1)
CALL zero33(advflwc,nc*(nr+2),1,1)
CALL zero33(advflc,1,nr,nc)

CALL zero33(advfups,1,nr,nc+1)
CALL zero33(advflws,1,nr,nc+1)
CALL zero33(advfls,1,nr,nc+1)

CALL zero33(advfupn,1,nr,nc+1)
CALL zero33(advflwn,1,nr,nc+1)
CALL zero33(advfln,1,nr,nc+1)

CALL zero33(ydiflv,1,nr,nc)
CALL zero33(xdiflv,1,nr,nc+1)
CALL zero33(ydiflu,1,nr+1,nc)
CALL zero33(xdiflu,1,nr,nc)


!     2. ADVECTIVE FLUXES AT CENTRES
!     ------------------------------

IF (iadvc > 0) THEN
  
  
!     2.1 UPWIND AND LW FLUXES
!     ------------------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        vdadv = vd2atc(j,i)
        vwt = SIGN(1.0,vdadv)
        cfl = vdadv*delt/(gy2(j)*h2atc(j,i))
        advfupc(j,i) = 0.5*vdadv *((1.0+vwt)*vd2(j,i)+(1.0-vwt)*vd2(j+1,i))
        advflwc(j,i) = 0.5*vdadv *((1.0+cfl)*vd2(j,i)+(1.0-cfl)*vd2(j+1,i))
      END IF
    END DO
  END DO
  
!     2.2 ADVECTIVE FLUXES
!     --------------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        psi = 0.0
        IF (iadvc == 1) THEN
!         ---fully upwind
          psi = 0.0
        ELSE IF (iadvc == 2) THEN
!         ---Lax-Wendroff
          psi = 1.0
        ELSE IF ((iadvc == 3).OR.(iadvc == 4)) THEN
!         ---TVD scheme
          vdadv = vd2atc(j,i)
          IF (vdadv > 0.0) THEN
            psi = fnlim(advflwc(j,i),  advfupc(j,i),  &
                advflwc(j-1,i),advfupc(j-1,i),iadvc)
          ELSE IF (vdadv < 0.0) THEN
            psi = fnlim(advflwc(j,i),  advfupc(j,i),  &
                advflwc(j+1,i),advfupc(j+1,i),iadvc)
          END IF
        END IF
        advflc(j,i) = (advfupc(j,i) + psi*(advflwc(j,i)-advfupc(j,i)))  &
            *cosphi(j)/h2atc(j,i)
      END IF
    END DO
  END DO
  
END IF

!     3. ADVECTIVE FLUXES AT U-NODES
!     ------------------------------

IF (iadvc > 0) THEN
  
  
!     3.1 UPWIND AND LW-FLUXES AT INTERIOR POINTS
!     -------------------------------------------
  
!      ---using VD2-values south of U-node
  DO  i=2,nc
    DO  j=2,nr
      IF (npix(j,i) == 1) THEN
        uwt = SIGN(1.0,ud2(j,i))
        cfl = ud2(j,i)*delt/(gx2u(j,i)*h2atu(j,i))
        advfups(j,i) = 0.5*ud2(j,i) *((1.0+uwt)*vd2(j,i-1)+(1.0-uwt)*vd2(j,i))
        advflws(j,i) = 0.5*ud2(j,i) *((1.0+cfl)*vd2(j,i-1)+(1.0-cfl)*vd2(j,i))
      END IF
    END DO
  END DO
  
!      ---using VD2-values north of U-node
  DO  i=2,nc
    DO  j=1,nr-1
      IF (npix(j,i) == 1) THEN
        uwt = SIGN(1.0,ud2(j,i))
        cfl = ud2(j,i)*delt/(gx2u(j,i)*h2atu(j,i))
        advfupn(j,i) = 0.5*ud2(j,i)  &
            *((1.0+uwt)*vd2(j+1,i-1)+(1.0-uwt)*vd2(j+1,i))
        advflwn(j,i) = 0.5*ud2(j,i)  &
            *((1.0+cfl)*vd2(j+1,i-1)+(1.0-cfl)*vd2(j+1,i))
      END IF
    END DO
  END DO
  
!     3.2 ADVECTIVE FLUXES AT INTERIOR POINTS
!     ---------------------------------------
  
!      ---using VD2-values south of U-node
  DO  i=2,nc
    DO  j=2,nr
      IF (npix(j,i) == 1) THEN
        psi = 0.0
        IF (iadvc == 1) THEN
!         ---fully upwind
          psi = 0.0
        ELSE IF (iadvc == 2) THEN
!         ---Lax-Wendroff
          psi = 1.0
        ELSE IF ((iadvc == 3).OR.(iadvc == 4)) THEN
!         ---TVD-scheme
          udadv = ud2(j,i)
          IF (udadv > 0.0) THEN
            psi = fnlim(advflws(j,i),  advfups(j,i),  &
                advflws(j,i-1),advfups(j,i-1),iadvc)
          ELSE IF (udadv < 0.0) THEN
            psi = fnlim(advflws(j,i),  advfups(j,i),  &
                advflws(j,i+1),advfups(j,i+1),iadvc)
          END IF
        END IF
        advfls(j,i) = (advfups(j,i) +  &
            psi*(advflws(j,i)-advfups(j,i)))/h2atu(j,i)
      END IF
    END DO
  END DO
  
!      ---using VD2-values north of U-node
  DO  i=2,nc
    DO  j=1,nr-1
      IF (npix(j,i) == 1) THEN
        psi = 0.0
        IF (iadvc == 1) THEN
!         ---fully upwind
          psi = 0.0
        ELSE IF (iadvc == 2) THEN
!         ---Lax-Wendroff
          psi = 1.0
        ELSE IF ((iadvc == 3).OR.(iadvc == 4)) THEN
!         ---TVD-scheme
          udadv = ud2(j,i)
          IF (udadv > 0.0) THEN
            psi = fnlim(advflwn(j,i),  advfupn(j,i),  &
                advflwn(j,i-1),advfupn(j,i-1),iadvc)
          ELSE IF (udadv < 0.0) THEN
            psi = fnlim(advflwn(j,i),  advfupn(j,i),  &
                advflwn(j,i+1),advfupn(j,i+1),iadvc)
          END IF
        END IF
        advfln(j,i) = (advfupn(j,i) +  &
            psi*(advflwn(j,i)-advfupn(j,i)))/h2atu(j,i)
      END IF
    END DO
  END DO

!     3.3 ADVECTIVE FLUXES AT U-OPEN BOUNDARIES
!     -----------------------------------------
  
  DO  ii=1,nobu
    i = iobu(ii)
    j = jobu(ii)

    IF (westob(ii)) THEN
      IF (ud2(j,i) <= 0.0) THEN
        advfls(j,i) = ud2(j,i)*vd2(j,i)/h2atc(j,i)
        advfln(j,i) = ud2(j,i)*vd2(j+1,i)/h2atc(j,i)
      ELSE
        advfls(j,i) = advfls(j,i+1)
        advfln(j,i) = advfln(j,i+1)
      END IF
    ELSE
      IF (ud2(j,i) >= 0.0) THEN
        advfls(j,i) = ud2(j,i)*vd2(j,i-1)/h2atc(j,i-1)
        advfln(j,i) = ud2(j,i)*vd2(j+1,i-1)/h2atc(j,i-1)
      ELSE
        advfls(j,i) = advfls(j,i-1)
        advfln(j,i) = advfln(j,i-1)
      END IF
    END IF
  END DO
  
END IF


!     4. DIFFUSIVE FLUXES
!     -------------------

IF (iodif > 0) THEN
  
  
!     4.1 INTERIOR POINTS
!     -------------------
  
!      ---(2,2)-component at centres
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        IF (npiy(j,i) == 1) THEN
          ydiflv(j,i) = ydiflv(j,i) - dheddyvc(j,i)*  &
              vd2(j,i)/(h2atv(j,i)*gy2(j))
        ELSE IF (npiy(j,i) > 1) THEN
          ydiflv(j,i) = ydiflv(j,i) - dheddyvc(j,i)*  &
              vd2(j,i)/(h2atc(j,i)*gy2(j))
        END IF
        IF (npiy(j+1,i) == 1) THEN
          ydiflv(j,i) = ydiflv(j,i) + dheddyvc(j,i)*  &
              vd2(j+1,i)/(h2atv(j+1,i)*gy2(j))
        ELSE IF (npiy(j+1,i) > 1) THEN
          ydiflv(j,i) = ydiflv(j,i) + dheddyvc(j,i)*  &
              vd2(j+1,i)/(h2atc(j,i)*gy2(j))
        END IF
        ydiflv(j,i) = 2.0*ydiflv(j,i)*cosphi(j)
      END IF
    END DO
  END DO
  
!      ---(1,2)-component at U-nodes
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) THEN
        xdiflv(j,i) = dheddyvu(j,i)*((vd2atc(j,i)/h2atc(j,i)  &
            -vd2atc(j,i-1)/h2atc(j,i-1))/gx2u(j,i) +sphcur(j)*ud2(j,i)/h2atu(j,i))
      END IF
    END DO
  END DO
  
!      ---(1,2)-component at V-nodes
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        ydiflu(j,i) = dheddyvv(j,i)*(ud2atc(j,i)/h2atc(j,i)  &
            -ud2atc(j-1,i)/h2atc(j-1,i))/gy2v(j)
      END IF
    END DO
  END DO
  
!      ---(1,1)-component at centres
  IF (igtrh == 1) THEN
    DO  i=1,nc
      DO  j=1,nr
        IF (nwd(j,i) == 1) THEN
          IF (npix(j,i) == 1) THEN
            xdiflu(j,i) = xdiflu(j,i) - ud2(j,i)/h2atu(j,i)
          ELSE IF (npix(j,i) > 1) THEN
            xdiflu(j,i) = xdiflu(j,i) - ud2(j,i)/h2atc(j,i)
          END IF
          IF (npix(j,i+1) == 1) THEN
            xdiflu(j,i) = xdiflu(j,i) + ud2(j,i+1)/h2atu(j,i+1)
          ELSE IF (npix(j,i+1) > 1) THEN
            xdiflu(j,i) = xdiflu(j,i) + ud2(j,i+1)/h2atc(j,i)
          END IF
          xdiflu(j,i) = 2.0*xdiflu(j,i)*dheddyvc(j,i)/gx2(j,i)
        END IF
      END DO
    END DO
  END IF
  
  
!     4.2 DIFFUSIVE FLUXES AT OPEN BOUNDARIES
!     ---------------------------------------
  
!      ---(1,2)-component at U-nodes
  DO  ii=1,nobu
    i = iobu(ii)
    j = jobu(ii)
    IF (westob(ii)) THEN
      xdiflv(j,i) = xdiflv(j,i+1)
    ELSE
      xdiflv(j,i) = xdiflv(j,i-1)
    END IF
  END DO
  
!      ---(1,2)-component at V-nodes
  DO  jj=1,nobv
    i = iobv(jj)
    j = jobv(jj)
    IF (soutob(jj)) THEN
      ydiflu(j,i) = ydiflu(j+1,i)
    ELSE
      ydiflu(j,i) = ydiflu(j-1,i)
    END IF
  END DO
  
END IF

!     5. ADVECTIVE AND DIFFUSIVE TERMS
!     --------------------------------

!     5.1 ADVECTION
!     -------------

IF (iadvc > 0) THEN
  
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        yadvv = (advflc(j,i)-advflc(j-1,i))/(gy2v(j)*cosphiv(j))
        xadvv = 0.5*(advfls(j,i+1) + advfln(j-1,i+1) -  &
            advfls(j,i)   - advfln(j-1,i))/gx2v(j,i)
        vah2d(j,i) = yadvv + xadvv + sphcurv(j)*ud2atv(j,i)**2 /h2atv(j,i)
      END IF
    END DO
  END DO
  
END IF


!     5.2 DIFFUSION
!     -------------

IF (iodif > 0) THEN
  
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) THEN
        ydifv = (ydiflv(j,i)-ydiflv(j-1,i))/(gy2v(j)*cosphiv(j))
        xdifv = 0.5*(xdiflv(j,i+1) + xdiflv(j-1,i+1) -  &
            xdiflv(j,i)   - xdiflv(j-1,i))/gx2v(j,i)
        IF (i == 1) THEN
          ydifu = (ydiflu(j,i+1) - ydiflu(j,i))  &
              /(0.5*gx2v(j,i+1)+1.5*gx2v(j,i))
        ELSE IF (i == nc) THEN
          ydifu = (ydiflu(j,i) - ydiflu(j,i-1))  &
              /(0.5*gx2v(j,i-1)+1.5*gx2v(j,i))
        ELSE
          ydifu = (ydiflu(j,i+1) - ydiflu(j,i-1))  &
              /(0.5*(gx2v(j,i-1)+gx2v(j,i+1))+gx2v(j,i))
        END IF
        vdh2d(j,i) = ydifv + xdifv + ydifu
        vdh2d(j,i) = vdh2d(j,i) + sphcurv(j)*  &
            (0.5*(xdiflu(j-1,i)+xdiflu(j,i))  &
            -2.0*sphcurv(j)*dheddyvv(j,i)*vd2(j,i)/h2atv(j,i))
      END IF
    END DO
  END DO
  
END IF


RETURN

END SUBROUTINE had2dv

!

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:39:42

SUBROUTINE heat
!************************************************************************

!    *HEAT*        SOLVE THE TEMPERATURE EQUATION

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 3 Sep 1999 @(COHERENS)heat.f 8.4

!       DESCRIPTION - SOLAR HEAT ABSORPTION SELECTED BY IOPTHE
!                     = 1 -> ABSORBED AT SEA SURFACE
!                     = 2 -> ABSORBED WITHIN A SURFACE LAYER
!                   - UPDATE TEMPERATURE BY CALLING TRANSPC

!       REFERENCE - Sections III-1.1, III-4.4 of the User Documentation

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS - IRRAD, TRANSPC, ZERO

!************************************************************************
!*    LOCAL VARIABLES

LOGICAL :: call1, flag
INTEGER :: i, j, k
REAL :: sfluxt(nr,nc), sink(nz,nr,nc), source(nz,nr,nc)

SAVE call1, sink
DATA call1 /.true./

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL
!    *FLAG*     LOGICAL   .TRUE. IF TRANSPC IS CALLED FOR THE FIRST TIME
!    *SFLUXT*   REAL      SURFACE FLUX OF TEMPERATURE [deg C m/s]
!    *SINK*     REAL      SINK TERMS IN TEMPERATURE EQUATION
!    *SOURCE*   REAL      SOURCE TERMS IN TEMPERATURE EQUATION

!------------------------------------------------------------------------

!     1. INITIALISATION
!     -----------------


!WRITE (*,*) 'TRACE : heat'

IF (call1) THEN
  
!     --- no sink terms
  CALL zero33(sink,nz,nr,nc)
  
END IF
flag = call1
call1 = .false.



!     2. BOUNDARY CONDITIONS (SURFACE)
!     --------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      sfluxt(j,i) = 0.
    END IF
  END DO
END DO


!     3. SOURCE TERMS
!     ---------------

!     3.1 EVALUATE SOLAR IRRADIANCE
!     -----------------------------

! IF (iopthe == 2) CALL irrad


!     3.2 HEATING SOURCE TERMS
!     ------------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        IF (iopthe == 2) THEN
          source(k,j,i) = 0.0
        ELSE
          source(k,j,i) = 0.0
        END IF
      END DO
    END IF
  END DO
END DO


!     4. UPDATE TEMPERATURE
!     ---------------------

CALL transpc(t,veddyd,source,sink, flag,0,0,iadvs,iadvs,iodif,  &
    ivpobu,ivpobv,tvp,w2, sfluxt,zeros2,zeros3,  &
    zeros4,zeros5,zeros6)

RETURN

END SUBROUTINE heat



SUBROUTINE euler
!************************************************************************

!    *EULER*        SOLVE ADV-DIF EQUATION FOR NON-BUOYANT EULERIAN TRACER

!       AUTHOR - Jochen Kaempf

!       LAST UPDATE - 16 Sep 2006
!************************************************************************
!*    LOCAL VARIABLES

LOGICAL :: call1, flag
INTEGER :: i, j, k
REAL :: sfluxt(nr,nc), sink(nz,nr,nc), source(nz,nr,nc)

SAVE call1, sink
DATA call1 /.true./

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL
!    *FLAG*     LOGICAL   .TRUE. IF TRANSPC IS CALLED FOR THE FIRST TIME
!    *SFLUXT*   REAL      SURFACE FLUX
!    *SINK*     REAL      SINK TERMS IN EQUATION
!    *SOURCE*   REAL      SOURCE TERMS IN EQUATION

!------------------------------------------------------------------------

!     1. INITIALISATION
!     -----------------

IF (call1) THEN
  
!     --- no sink terms
  CALL zero33(sink,nz,nr,nc)
  
END IF
flag = call1
call1 = .false.

!     2. BOUNDARY CONDITIONS (SURFACE)
!     --------------------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      sfluxt(j,i) = 0.
    END IF
  END DO
END DO

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
          source(k,j,i) = 0.0
      END DO
    END IF
  END DO
END DO
!
!     4a. UPDATE TRACER 1
!     ---------------------

trigger = 1.0
CALL transpc(conc1,veddyd,source,sink, flag,0,0,iadvs,iadvs,iodif,  &
    ivpobu,ivpobv,cvp,w2, sfluxt,zeros2,zeros3,  &
    zeros4,zeros5,zeros6)
trigger  = 0.0

!     4a. UPDATE TRACER 2
!     ---------------------

trigger = 2.0
CALL transpc(conc2,veddyd,source,sink, flag,0,0,iadvs,iadvs,iodif,  &
    ivpobu,ivpobv,cvp,w2, sfluxt,zeros2,zeros3,  &
    zeros4,zeros5,zeros6)
trigger  = 0.0

!     4a. UPDATE TRACER 3
!     ---------------------

CALL transpc(conc3,veddyd,source,sink, flag,0,0,iadvs,iadvs,iodif,  &
    ivpobu,ivpobv,cvp,w2, sfluxt,zeros2,zeros3,  &
    zeros4,zeros5,zeros6)

!     4a. UPDATE TRACER 4
!     ---------------------

CALL transpc(conc4,veddyd,source,sink, flag,0,0,iadvs,iadvs,iodif,  &
    ivpobu,ivpobv,cvp,w2, sfluxt,zeros2,zeros3,  &
    zeros4,zeros5,zeros6)

RETURN

END SUBROUTINE euler

!=======================================================================

SUBROUTINE irrad
!***********************************************************************

!    *IRRAD*        HEAT ABSORPTION TERM IN THE TEMPERATURE EQUATION

!       AUTHOR - PATRICK LUYTEN AND PAUL TETT

!       LAST UPDATE - 26 Mar 1999 @(COHERENS)heat.f 8.4

!       DESCRIPTION - UPDATE DIFFUSE ATTENUATION COEFFICIENT FOR THE PHYSICS
!                   - EVALUATE ABSORPTION SOURCE TERM IN TEMPERATURE EQUATION

!       REFERENCE - Sections III-1.4 and III-4.4.10 of the User Documentation

!       CALLING PROGRAM - HEAT

!       EXTERNALS -

!***********************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, j, k
REAL :: radlo1, radlo2, radup1, radup2, r2exp, zbot, ztop
!REAL :: h2atc

!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *RADLO1*   REAL      (INFRARED PART OF) IRRADIANCE AT VERTICAL
!                         NODE (BELOW GRID POINT)          [W/m2]
!    *RADLO2*   REAL      (NON-INFRARED PART OF) IRRADIANCE AT VERTICAL
!                         NODE (BELOW GRID POINT)          [W/m2]
!    *RADUP1*   REAL      (INFRARED PART OF) IRRADIANCE AT VERTICAL
!                         NODE (ABOVE GRID POINT)          [W/m2]
!    *RADUP2*   REAL      (NON-INFRARED PART OF) IRRADIANCE AT VERTICAL
!                         NODE (ABOVE GRID POINT)          [W/m2]
!    *R2EXP*    REAL      EXTRA ATTENUATION FACTOR

!-----------------------------------------------------------------------

!     1. ATTENUATION COEFFICIENT
!     --------------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        atcfcor2(k,j,i) = atcf2(j,i) - epssal*s(k,j,i)
      END DO
    END IF
  END DO
END DO


!     2. HEAT ABSORPTION TERM
!     -----------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      radup1 = r1opt(j,i)*qsol(j,i)
      radup2 = (1.0-r1opt(j,i))*qsol(j,i)
      DO  k=nz,1,-1
        radlo1 = EXP(-atcf1(j,i)*gz2(k,j,i))*radup1
        ztop = h2atc(j,i)*(1.0-gz0(k+1,j,i))
        zbot = h2atc(j,i)*(1.0-gz0(k,j,i))
        IF (zbot <= hexp(j,i)) THEN
          r2exp = r2opt(j,i)**(gz2(k,j,i)/hexp(j,i))
        ELSE IF (ztop >= hexp(j,i)) THEN
          r2exp = 1.0
        ELSE
          r2exp = r2opt(j,i)**(1.0-ztop/hexp(j,i))
        END IF
        radlo2 = r2exp*EXP(-atcfcor2(k,j,i)*gz2(k,j,i))*radup2
        IF (k > 1) THEN
          qheat(k,j,i) = (radup1+radup2-radlo1-radlo2) /gz2(k,j,i)
        ELSE
          qheat(k,j,i) = (radup1+radup2)/gz2(k,j,i)
        END IF
        radup1 = radlo1
        radup2 = radlo2
      END DO
    END IF
  END DO
END DO


RETURN

END SUBROUTINE irrad

!

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:39:48

SUBROUTINE heddy
!************************************************************************

!    *HEDDY*       EVALUATE HORIZONTAL DIFFUSION COEFFICIENTS

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 26 Mar 1999       @(COHERENS)heddy.f 8.4

!       DESCRIPTION - EVALUATES 2-D (IODIF=1,2) AND 3-D (IODIF=2) HORIZONTAL
!                     DIFFUSION COEFFICIENTS AT CELL CENTRES, U-NODES, V-NODES

!       REFERENCE -  Section III-1.3, III-4.3.4c, III-4.4.3c, III-4.6.2 of the
!                    User Documentation

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS -

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, ii, j, jj, k
REAL :: dudxc(nz,nr,nc), dudyv(nz,nr+1,nc)
REAL :: dvdyc(nz,nr,nc), dvdxu(nz,nr,nc+1)
REAL :: deform, fac, ucor, vcor
!REAL :: gx2u, gx2v, gy2v, gz2u, gz2v, h2atc, h2atu, h2atv
!REAL :: u2atc, u2atv, v2atc, v2atu

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *DUDXC*    REAL      X-DERIVATIVE OF U2 AT CELL CENTRE   [1/s]
!    *DUDYV*    REAL      Y-DERIVATIVE OF U2 AT V-NODE        [1/s]
!    *DVDYC*    REAL      Y-DERIVATIVE OF V2 AT CELL CENTRE   [1/s]
!    *DVDXU*    REAL      X-DERIVATIVE OF V2 AT U-NODE        [1/s]
!    *DEFORM*   REAL      TOTAL STRAIN RATE                   [1/s]
!    *UCOR*     REAL      SPHERICAL CORRECTION TERM           [1/s]
!    *VCOR*     REAL      SPHERICAL CORRECTION TERM           [1/s]

!------------------------------------------------------------------------

!     1. DEPTH-INTEGRATED VALUES IF IODIF=1
!     -------------------------------------


IF (iodif == 1) THEN
  
  
!     1.1 CENTRES
!     -----------
  
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) dheddyvc(j,i) = h2atc(j,i)
    END DO
  END DO
  
  
!     1.2 U-NODES
!     ----------
  
!        ---interior points
  DO  i=2,nc
    DO  j=1,nr
      IF (npix(j,i) == 1) dheddyvu(j,i) = h2atu(j,i)
    END DO
  END DO
  
!        ---U-open boundaries
  DO  ii=1,nobu
    i = iobu(ii)
    j = jobu(ii)
    IF (westob(ii)) THEN
      dheddyvu(j,i) = h2atc(j,i)
    ELSE
      dheddyvu(j,i) = h2atc(j,i-1)
    END IF
  END DO
  
  
!     1.3 V-NODES
!     ----------
  
!        ---interior points
  DO  i=1,nc
    DO  j=2,nr
      IF (npiy(j,i) == 1) dheddyvv(j,i) = h2atv(j,i)
    END DO
  END DO
  
!        ---V-open boundaries
  DO  jj=1,nobv
    i = iobv(jj)
    j = jobv(jj)
    IF (soutob(jj)) THEN
      dheddyvv(j,i) = h2atc(j,i)
    ELSE
      dheddyvv(j,i) = h2atc(j-1,i)
    END IF
  END DO
  
  GO TO 1000
  
END IF


!     2. INITIALISE ARRAYS
!     --------------------

CALL zero33(dudxc,nz,nr,nc)
CALL zero33(dvdyc,nz,nr,nc)
CALL zero33(dudyv,nz,nr+1,nc)
CALL zero33(dvdxu,nz,nr,nc+1)


!     3. DERIVATIVES AT INTERIOR POINTS
!     ---------------------------------

!     3.1 CENTRES
!     -----------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        dudxc(k,j,i) = (u2(k,j,i+1)-u2(k,j,i))/gx2(j,i)
        dvdyc(k,j,i) = (v2(k,j+1,i)-v2(k,j,i))/gy2(j)
      END DO
    END IF
  END DO
END DO
311   CONTINUE


!     3.2 V-NODES
!     -----------

DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      DO  k=1,nz
        dudyv(k,j,i) = (u2atc(k,j,i)-u2atc(k,j-1,i))/gy2v(j)
      END DO
    END IF
  END DO
END DO
321   CONTINUE


!     3.3 U-NODES
!     -----------

DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      DO  k=1,nz
        dvdxu(k,j,i) = (v2atc(k,j,i)-v2atc(k,j,i-1))/gx2u(j,i)
      END DO
    END IF
  END DO
END DO
331   CONTINUE


!     4. DERIVATIVES AT OPEN BOUNDARIES
!     ---------------------------------

!     4.1 V-NODES
!     -----------

DO  jj=1,nobv
  i = iobv(jj)
  j = jobv(jj)
  DO  k=1,nz
    IF (soutob(jj)) THEN
      dudyv(k,j,i) = dudyv(k,j+1,i)
    ELSE
      dudyv(k,j,i) = dudyv(k,j-1,i)
    END IF
  END DO
END DO
411   CONTINUE


!     4.2 U-NODES
!     -----------

DO  ii=1,nobu
  i = iobu(ii)
  j = jobu(ii)
  DO  k=1,nz
    IF (westob(ii)) THEN
      dvdxu(k,j,i) = dvdxu(k,j,i+1)
    ELSE
      dvdxu(k,j,i) = dvdxu(k,j,i-1)
    END IF
  END DO
END DO
421   CONTINUE


!     5. EVALUATE HORIZONTAL DIFFUSION COEFFICIENTS
!     ---------------------------------------------

!     5.1 CENTRES
!     -----------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      fac = gx2(j,i)*gy2(j)
      dheddyvc(j,i) = 0.0
      DO  k=1,nz
        ucor = sphcur(j)*u2atc(k,j,i)
        vcor =-sphcur(j)*v2atc(k,j,i)
        deform = (dudxc(k,j,i)+vcor)**2 + dvdyc(k,j,i)**2  &
            + 0.5*(0.5*(dudyv(k,j,i)+dudyv(k,j+1,i)  &
            +dvdxu(k,j,i)+dvdxu(k,j,i+1)) + ucor)**2
        heddyvc(k,j,i) = cm0*fac*SQRT(deform)
        heddydc(k,j,i) = cs0*fac*SQRT(deform)
        dheddyvc(j,i) = dheddyvc(j,i)+heddyvc(k,j,i)*gz2(k,j,i)
      END DO
    END IF
  END DO
END DO
511   CONTINUE


!     5.2 U-NODES
!     -----------

!     ---interior points
DO  i=2,nc
  DO  j=1,nr
    IF (npix(j,i) == 1) THEN
      fac = gx2u(j,i)*gy2(j)
      dheddyvu(j,i) = 0.0
      DO  k=1,nz
        ucor = sphcur(j)*u2(k,j,i)
        vcor =-sphcur(j)*v2atu(k,j,i)
        deform = (0.5*(dudxc(k,j,i)+dudxc(k,j,i-1))+vcor)**2  &
            +0.25*(dvdyc(k,j,i)+dvdyc(k,j,i-1))**2  &
            +0.5*(0.25*(dudyv(k,j,i)+dudyv(k,j,i-1)  &
            +dudyv(k,j+1,i)+dudyv(k,j+1,i-1)) + dvdxu(k,j,i) + ucor)**2
        heddyvu(k,j,i) = cm0*fac*SQRT(deform)
        heddydu(k,j,i) = cs0*fac*SQRT(deform)
        dheddyvu(j,i) = dheddyvu(j,i) + heddyvu(k,j,i)*gz2u(k,j,i)
      END DO
    END IF
  END DO
END DO
521   CONTINUE

!     ---U-open boundaries
DO  ii=1,nobu
  i = iobu(ii)
  j = jobu(ii)
  IF (westob(ii)) THEN
    dheddyvu(j,i) = dheddyvc(j,i)
    DO  k=1,nz
      heddyvu(k,j,i) = heddyvc(k,j,i)
      heddydu(k,j,i) = heddydc(k,j,i)
    END DO
  ELSE
    dheddyvu(j,i) = dheddyvc(j,i-1)
    DO  k=1,nz
      heddyvu(k,j,i) = heddyvc(k,j,i-1)
      heddydu(k,j,i) = heddydc(k,j,i-1)
    END DO
  END IF
END DO
522   CONTINUE


!     5.3 V-NODES
!     -----------

!     ---interior points
DO  i=1,nc
  DO  j=2,nr
    IF (npiy(j,i) == 1) THEN
      fac = gx2v(j,i)*gy2v(j)
      dheddyvv(j,i) = 0.0
      DO  k=1,nz
        ucor = sphcurv(j)*u2atv(k,j,i)
        vcor =-sphcurv(j)*v2(k,j,i)
        deform = (0.5*(dudxc(k,j,i)+dudxc(k,j-1,i))+vcor)**2  &
            +0.25*(dvdyc(k,j,i)+dvdyc(k,j-1,i))**2 +0.5*(dudyv(k,j,i) + ucor +  &
            0.25*(dvdxu(k,j,i)+dvdxu(k,j-1,i)  &
            +dvdxu(k,j,i+1)+dvdxu(k,j-1,i+1)))**2
        heddyvv(k,j,i) = cm0*fac*SQRT(deform)
        heddydv(k,j,i) = cs0*fac*SQRT(deform)
        dheddyvv(j,i) = dheddyvv(j,i) + heddyvv(k,j,i)*gz2v(k,j,i)
      END DO
    END IF
  END DO
END DO
531   CONTINUE

!     ---V-open boundaries
DO  jj=1,nobv
  i = iobv(jj)
  j = jobv(jj)
  IF (soutob(jj)) THEN
    dheddyvv(j,i) = dheddyvc(j,i)
    DO  k=1,nz
      heddyvv(k,j,i) = heddyvc(k,j,i)
      heddydv(k,j,i) = heddydc(k,j,i)
    END DO
  ELSE
    dheddyvv(j,i) = dheddyvc(j-1,i)
    DO  k=1,nz
      heddyvv(k,j,i) = heddyvc(k,j-1,i)
      heddydv(k,j,i) = heddydc(k,j-1,i)
    END DO
  END IF
END DO
532   CONTINUE


1000 CONTINUE

RETURN

END SUBROUTINE heddy

!
!************************************************************************

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:41:40

!    *metin.f*       ENSEMBLE OF ROUTINES TO READ MET FORCING AND TO EVALUATE
!                    SURFACE FLUXES AND SOLAR IRRADIANCE

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 25 May 1999 @(COHERENS)metin.f 8.4

!       DESCRIPTION -

!       REFERENCE -

!       SUBROUTINES - FLUXCO, HUMID, METIN, SOLRAD, SURFLX

!       FUNCTIONS - CD, CHARNO, CE

!************************************************************************

!========================================================================



SUBROUTINE salt
!************************************************************************

!    *SALT*       SOLVE THE SALINITY EQUATION

!       AUTHOR - KEVIN RUDDICK AND PATRICK LUYTEN

!       LAST UPDATE - 19 Apr 1999 @(COHERENS)salt.f 8.4

!       DESCRIPTION - SET ALL SOURCE/SINK TERMS TO ZERO ON FIRST CALL
!                   - UPDATE SALINITY BY CALLING TRANSPC

!       REFERENCES -  Sections III-1.1, III-4.4 of the User Documentation

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS - TRANSPC, ZERO

!************************************************************************
!*    LOCAL VARIABLES

LOGICAL :: call1, flag
REAL :: sink(nz,nr,nc), source(nz,nr,nc)
INTEGER :: j

SAVE call1, source, sink
DATA call1 /.true./

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL
!    *FLAG*     LOGICAL   .TRUE. IF TRANSPC IS CALLED FOR THE FIRST TIME
!    *SOURCE*   REAL      SOURCE TERMS IN SALINITY EQUATION (=0)
!    *SINK*     REAL      SINK TERMS IN SALINITY EQUATION (=0)

!------------------------------------------------------------------------

!     1. INITIALISATION
!     -----------------

!WRITE (*,*) 'TRACE : salt'

IF (call1) THEN
  
!        ---no source/sink terms
  CALL zero33(source,nz,nr,nc)
  CALL zero33(sink,nz,nr,nc)
  
END IF
flag = call1
call1 =.false.


!     2. UPDATE SALINITY
!     ------------------

CALL transpc(s,veddyd,source,sink, flag,0,0,iadvs,iadvs,iodif,  &
    ivpobu,ivpobv,svp,w2, ssalfl,zeros2,zeros3,  &
    zeros4,zeros5,zeros6)

RETURN

END SUBROUTINE salt

!

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:43:05

SUBROUTINE thomv(a,b,c,d,f,nk,nj,ni)
!************************************************************************

!    *THOMV*      SOLVE TRIDIAGONAL LINEAR EQUATION SYSTEM
!                 (COUPLED ON VERTICAL)

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 7 Dec 1994 @(COHERENS)thomv.f 6.1

!       DESCRIPTION - USES THOMAS ALGORITHM TO SOLVE

!                        B(I,J,1) *F(I,J,1) +C(I,J,1)*F(I,J,2)  =D(I,J,1)
!  A(I,J,K) *F(I,J,K-1) +B(I,J,K) *F(I,J,K) +C(I,J,K)*F(I,J,K+1)=D(I,J,K)
!  A(I,J,NK)*F(I,J,NK-1)+B(I,J,NK)*F(I,J,NK)                    =D(I,J,NK)

!        REFERENCE - Press W.H., Flannery B.P., Teukolsky S.A. and Vetterling
!                    W.T., 1989. Numerical Recipes. The art of scientific
!                    computing. Cambridge University Press, Cambridge, 702 pp.

!       CALLING PROGRAM - TRANSP, TRANSPW, UCALC, VCALC

!       EXTERNALS -

!************************************************************************
!*    ARGUMENTS

INTEGER, INTENT(IN)                      :: nk
INTEGER, INTENT(IN)                      :: nj
INTEGER, INTENT(IN)                      :: ni
REAL, INTENT(IN)                         :: a(nk,nj,ni)
REAL, INTENT(IN)                         :: b(nk,nj,ni)
REAL, INTENT(IN)                         :: c(nk,nj,ni)
REAL, INTENT(IN)                         :: d(nk,nj,ni)
REAL, INTENT(OUT)                        :: f(nk,nj,ni)

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *A*        REAL      THOMAS COEFF FOR LOWER  POINT (K-1 LEVEL)  (IN)
!    *B*        REAL      THOMAS COEFF FOR CENTRE POINT (K LEVEL)    (IN)
!    *C*        REAL      THOMAS COEFF FOR UPPER POINT  (K+1 LEVEL)  (IN)
!    *D*        REAL      THOMAS COEFF FOR RIGHT HAND SIDE (EXPLICIT TERMS)
!                                                                    (IN)
!    *F*        REAL      SOLUTION OF TRIDIAGONAL SYSTEM         (IN/OUT)
!    *NI*       INTEGER   DIMENSION OF F IN X DIRECTION              (IN)
!    *NJ*       INTEGER   DIMENSION OF F IN Y DIRECTION              (IN)
!    *NK*       INTEGER   DIMENSION OF F IN Z DIRECTION              (IN)

!------------------------------------------------------------------------

!*    LOCAL VARIABLES

INTEGER :: i, j, k
REAL :: gamma(nz+1)
REAL :: beta

!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *BETA*     REAL      TEMPORARY VARIABLE
!    *GAMMA*    REAL      TEMPORARY VARIABLE

!-----------------------------------------------------------------------

DO  i=1,ni
  DO  j=1,nj
    beta = b(1,j,i)
    f(1,j,i) = d(1,j,i)/ beta
    DO  k=2,nk
      gamma(k) = c(k-1,j,i)/ beta
      beta = b(k,j,i) -a(k,j,i)*gamma(k)
      f(k,j,i) = (d(k,j,i)-a(k,j,i)*f(k-1,j,i))/ beta
    END DO
    DO  k=nk-1,1,-1
      f(k,j,i) = f(k,j,i) - gamma(k+1)*f(k+1,j,i)
    END DO
  END DO
END DO


RETURN

END SUBROUTINE thomv

!

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:43:27

SUBROUTINE tkleng
!***********************************************************************

!      *TKLENG*      SOLVE TRANSPORT EQUATION FOR K*L

!       AUTHOR - PATIRCK LUYTEN

!       LAST UPDATE - 7 May 1998        @(COHERENS)tkleng.f 8.4

!       DESCRIPTION - STORE OLD VALUES OF K*L
!                   - APPLY SURFACE AND BOTTOM BOUNDARY CONDITIONS
!                   - EVALUATE SOURCE/SINK TERMS
!                   - UPDATE K*L BY CALLING TRANSPW
!                   - UPDATE L BY DIVIDING K*L BY K

!       REFERENCE - Sections III-1.2.2b and III-4.4.5.3 of the
!                   User Documentation
!                 - Mellor G.L. and Yamada T., 1982. Development of a
!                   turbulence closure model for geophysical fluid problems.
!                    Reviews of Geophys and Space Physics, 20, 851-875.

!       CALLING PROGRAM - VEDDY2

!       EXTERNALS - TRANSPW

!***********************************************************************
!*    LOCAL VARIABLES

LOGICAL :: call1, flag
INTEGER :: i, iadvt, ibot, idift, isur, j, k

REAL :: btkezl(nr,nc)
REAL :: sink(nz+1,nr,nc), source(nz+1,nr,nc)
REAL :: stkezl(nr,nc), tkezl(nz+1,nr,nc)
REAL :: wpf, zl1, zl2
!REAL :: h2atc

SAVE call1
DATA call1 /.true./

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL
!    *FLAG*     LOGICAL   .TRUE. IF TRANSPW IS CALLED FOR THE FIRST TIME
!    *IADVT*    INTEGER   TYPE OF SCHEME FOR ADVECTION (0/1/2/3/4)
!    *IBOT*     INTEGER   =2 (DIRICHLET NEAR BOTTOM B.C.)
!    *IDIFT*    INTEGER   TYPE OF SCHEME FOR HORIZONTAL DIFFUSION (0/1/2)
!    *ISUR*     INTEGER   =2 (DIRICHLET NEAR SURFACE B.C.)
!    *BTKEZL*   REAL      NEAR BOTTOM VALUE OF KL
!    *SINK*     REAL      SINK TERMS IN KL-EQUATION
!    *SOURCE*   REAL      SOURCE TERMS IN KL-EQUATION
!    *STKEZL*   REAL      NEAR SURFACE VALUE OF KL
!    *TKEZL*    REAL      =T.K.E. TIMES LENGTH SCALE
!    *WPF*      REAL      WALL PROXIMITY FUNCTION

!------------------------------------------------------------------------

!     1. INITIALISE PARAMETERS AND ARRAYS
!     -----------------------------------

!WRITE (*,*) 'TRACE: tkleng'

isur = 2
ibot = 2
IF (iahdht == 0) THEN
  iadvt = 0
  idift = 0
ELSE
  iadvt = iadvs
  idift = iodif
END IF

flag = call1
call1 = .false.


!     2. STORE OLD VALUES
!     -------------------


DO  i=1,nc
  DO  j=1,nr
    DO  k=1,nz+1
      tkezl(k,j,i) = tkeold(k,j,i)*zlw(k,j,i)
    END DO
  END DO
END DO


!     3. BOUNDARY CONDITIONS (SURFACE, BOTTOM)
!     ----------------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      stkezl(j,i) = tkew(nz,j,i) *ckar*(h2atc(j,i)*(1.0-gz0(nz,j,i))+z0sur)
      btkezl(j,i) = tkew(2,j,i) *ckar*(h2atc(j,i)*gz0(2,j,i)+z0bot)
    END IF
  END DO
END DO


!     4. SOURCE/SINK TERMS
!     --------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=3,nz-1
        IF (tkeold(k,j,i) > 0.0) THEN
          source(k,j,i) = 0.5*e1*zlw(k,j,i)*shprod(k,j,i)
          IF (buprod(k,j,i) > 0.0) THEN
            source(k,j,i) = source(k,j,i)+0.5*e1*zlw(k,j,i) *buprod(k,j,i)
            sink(k,j,i) = 0.0
          ELSE
            sink(k,j,i) = -0.5*e1*buprod(k,j,i)/tkeold(k,j,i)
          END IF
          IF (zlw(k,j,i) > 0.0) THEN
            zl1 = gz0(k,j,i)*h2atc(j,i) + z0bot
            zl2 = h2atc(j,i)*(1.0-gz0(k,j,i)) + z0sur
            wpf = 1.0+e2*(zlw(k,j,i)*(zl1+zl2) /(ckar*zl1*zl2))**2
            sink(k,j,i) = sink(k,j,i)+0.5*eps0*wpf  &
                *SQRT(tkeold(k,j,i))/zlw(k,j,i)
          END IF
        ELSE
          source(k,j,i) = 0.0
          sink(k,j,i)   = 0.0
        END IF
      END DO
    END IF
  END DO
END DO


!     5. UPDATE TKEZL
!     ---------------


CALL transpw(tkezl,veddyk,source,sink, flag,isur,ibot,iadvt,iadvt,idift,  &
    ivpobu,ivpobv,turvp, zeros1,zeros2,stkezl,  &
    zeros4,zeros5,btkezl)


!     6. UPDATE ZLW
!     -------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        IF (tkew(k,j,i) > 0.0) THEN
          zlw(k,j,i) = tkezl(k,j,i)/tkew(k,j,i)
        ELSE
          zlw(k,j,i) = 0.0
        END IF
      END DO
    END IF
  END DO
END DO


RETURN

END SUBROUTINE tkleng

!=======================================================================

SUBROUTINE tlendis
!************************************************************************

!    *TLENDIS*  EVALUATE MIXING LENGTH (K-EPS THEORY)

!     AUTHOR - PATRICK LUYTEN

!     LAST UPDATE  - 7 May 1998         @(COHERENS)turben.f 8.4

!     DESCRIPTION - EVALUATES MIXING LENGTH [m] AS A FUNCTION OF
!                   DISSIPATION RATE

!     REFERENCE - Section III-1.2.2 of the User Documentation

!     CALLING PROGRAM - VEDDY2

!     EXTERNALS -

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, j, k

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        IF (dissw(k,j,i) > 0.0) THEN
          zlw(k,j,i) = eps0*tkew(k,j,i)**1.5/dissw(k,j,i)
        ELSE
          zlw(k,j,i) = 0.0
        END IF
      END DO
    END IF
  END DO
END DO


RETURN

END SUBROUTINE tlendis

!


SUBROUTINE transpc(phi,vedphi,source,sink,call1,  &
        isur,ibot,jadvhs,jadvws,jdifs,  &
        ivpu,ivpv,phivp,w2phi,  &
        sflux,stran,sphi,  &
        bflux,btran,bphi)

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:43:51

!***********************************************************************

!    *TRANSPC*  SOLVE TRANSPORT EQUATION FOR A CELL CENTERED SCALAR QUANTITY

!     AUTHOR - PATRICK LUYTEN AND KEVIN RUDDICK


!     LAST UPDATE - 26 Mar 1999         @(COHERENS)transp.f 8.4

!     DESCRIPTION - SOLVES ADVECTION-DIFFUSION EQUATION FOR PHI
!                   (SALINITY, TEMPERATURE, BIOLOGICAL QUANTITITES,
!                   CONTAMINANTS, SEDIMENT CONCENTRATIONS)
!                 - SOURCE TERMS ARE TAKEN EXPLICITLY
!                 - HORIZONTAL ADVECTION/DIFFUSION IS TAKEN EXPLICITLY
!                 - VERTICAL ADVECTION IS SEMI-IMPLICIT VERTICAL DIFFUSION
!                   FULLY IMPLICIT
!                 - ADVECTION SCHEME SELECTED BY IADVS
!                      = 1 => FIRST ORDER UPWIND
!                      = 2 => LAX-WENDROFF/CENTRAL
!                      = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                      = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!                 - HORIZONTAL DIFFUSION SELECTED BY IODIF
!                      = 0 => OFF
!                      = 1 => UNIFORM
!                      = 2 => SMAGORINSKY
!                 - USES CORRECTOR TERMS TO ELIMINATE THE JACOBIAN FROM
!                   THE TIME DERIVATIVE
!                 - USES THE "FILTERED" ADVECTIVE VELOCITES U2F, V2F
!                 - TIME INTEGRATION IS PERFORMED WITH OR WITHOUT OPERATOR
!                   SPLITTING DEPENDING ON THE TYPE OF ADVECTION SCHEME

!     REFERENCE - Sections III-4.4 and V-1.8 of the User Documentation

!     CALLING PROGRAM - BIOLGY, HEAT, MASS, SALT, SEDEUL

!     EXTERNALS - THOMV, XADVDIS, YADVDIS, ZADVDIS, ZERO

!***********************************************************************
!*    ARGUMENTS


REAL, INTENT(IN OUT)                        :: phi(nz,nr,nc)
REAL, INTENT(IN)                     :: vedphi(nz+1,nr,nc)
REAL, INTENT(IN)                         :: source(nz,nr,nc)
REAL, INTENT(IN)                         :: sink(nz,nr,nc)
LOGICAL, INTENT(IN)                  :: call1
INTEGER, INTENT(IN)                  :: isur
INTEGER, INTENT(IN)                  :: ibot
INTEGER, INTENT(IN)                  :: jadvhs
INTEGER, INTENT(IN)                  :: jadvws
INTEGER, INTENT(IN)                  :: jdifs
INTEGER, INTENT(IN)                  :: ivpu(0:nobu)
INTEGER, INTENT(IN)                  :: ivpv(0:nobv)
REAL, INTENT(IN)                     :: phivp(0:nz,0:nvprof)
REAL, INTENT(IN)                     :: w2phi(nz+1,nr,nc)
REAL, INTENT(IN)                     :: sflux(nr,nc)
REAL, INTENT(IN)                     :: stran(nr,nc)
REAL, INTENT(IN)                     :: sphi(nr,nc)
REAL, INTENT(IN)                     :: bflux(nr,nc)
REAL, INTENT(IN)                     :: btran(nr,nc)
REAL, INTENT(IN)                     :: bphi(nr,nc)

!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL (IN)
!    *IBOT*     INTEGER   BOTTOM B.C.  (SEE ZADVDIS)                  (IN)
!    *ISUR*     INTEGER   SURFACE B.C. (SEE ZADVDIS)                  (IN)
!    *JADVHS*   INTEGER   ADVECTION SCHEME FOR HORIZONTAL ADVECTION   (IN)
!                         = 0 => HORIZONTAL ADVECTION DISABLED
!                         = 1 => FIRST ORDER UPWIND
!                         = 2 => LAX-WENDROFF
!                         = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                         = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!    *JADVWS*   INTEGER   ADVECTION SCHEME FOR VERTICAL ADVECTION     (IN)
!                         = 0 => VERTICAL ADVECTION DISABLED
!                         = 1 => FIRST ORDER UPWIND
!                         = 2 => CENTRAL
!                         = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                         = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!    *JDIFS*    INTEGER   SWITCH TO SELECT SCHEME FOR HORIZONTAL DIFFUSION
!                         = 0 => HORIZONTAL DIFFUSION DISABLED        (IN)
!                         = 1 => UNIFORM DIFFSUSION COEFFICIENTS
!                         = 2 => "SMAGORINSKY" FORMULATION
!    *IVPU*     INTEGER   PROFILE NUMBERS AT U-O.B.                   (IN)
!    *IVPV*     INTEGER   PROFILE NUMBERS AT V-O.B.                   (IN)
!    *BFLUX*    REAL      (SEE ZADVDIS)                  [PHI*m/s]    (IN)
!    *BPHI*     REAL      BOTTOM VALUE OF PHI (SEE ZADVDIS)  [PHI]    (IN)
!    *BTRAN*    REAL      BOTTOM TRANSFER COEFF (SEE ZADVDIS)  [m/s]  (IN)
!    *PHI*      REAL      QUANTITY TO BE ADVECTED/DIFFUSED [PHI]  (IN/OUT)
!    *PHIVP*    REAL      OPEN BOUNDARY CONDITION [PHI]               (IN)
!                       PHIVP(K,*)>PHIVP(0,*)  : IMPOSED INFLOW
!                       PHIVP(K,*)<= PHIVP(0,*): NORMAL GRADIENT CONDITION
!    *SFLUX*    REAL      (SEE ZADVDIS)                   [PHI*m/s]   (IN)
!    *SINK*     REAL      SINK TERMS IN TRANSPORT EQUATION (=0)       (IN)
!    *SOURCE*   REAL      SOURCE TERMS IN TRANSPORT EQUATION          (IN)
!    *SPHI*     REAL      SURFACE VALUE OF PHI (SEE ZADVDIS) [PHI]    (IN)
!    *STRAN*    REAL      SURFACE TRANSFER COEFF (SEE ZADVDIS) [m/s]  (IN)
!    *VEDPHI*   REAL      VERTICAL DIFFUSION COEFFICIENT AT W NODES [m2/s]
!                                                                     (IN)
!    *W2PHI*    REAL      VERTICAL (PHYSICAL PLUS SINKING) VELOCITIES [m/s]
!                                                                     (IN)
!-----------------------------------------------------------------------

!*    LOCAL VARIABLES

INTEGER :: i, j, k
REAL :: phia(nz,nr,nc), phib(nz,nr,nc)
REAL :: plhsa(nz,nr,nc), plhsb(nz,nr,nc)
REAL :: plhsc(nz,nr,nc), prhs(nz,nr,nc)
REAL :: ucorr(nz,nr,nc), vcorr(nz,nr,nc), wcorr(nz,nr,nc)
!REAL :: gz2u, gz2v

!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *PLHSA*    REAL      LEFT HAND SIDE OF PHI EQUATION /PHI(I,J,K-1)
!    *PLHSB*    REAL      LEFT HAND SIDE OF PHI EQUATION /PHI(I,J,K)
!    *PLHSC*    REAL      LEFT HAND SIDE OF PHI EQUATION /PHI(I,J,K+1)
!    *PRHS*     REAL      RIGHT HAND SIDE OF PHI EQUATION
!    *PHIA*     REAL      FIRST FRACTIONAL STEP NEW PHI
!    *PHIB*     REAL      SECOND FRACTIONAL STEP NEW PHI
!    *UCORR*    REAL      CORRECTOR TERMS IN X-DIRECTION
!    *VCORR*    REAL      CORRECTOR TERMS IN Y-DIRECTION
!    *WCORR*    REAL      CORRECTOR TERMS IN Z-DIRECTION

!-----------------------------------------------------------------------

!     1. INITIALISE AND CHECKING
!     --------------------------

!     1.1 INITIALISE PHI AT BLOCKED POINTS ON FIRST CALL
!     --------------------------------------------------

IF (call1) THEN
  DO  i=1,nc
    DO  j=1,nr
      DO  k=1,nz
        IF (nwd(j,i) == 0) THEN
          phi(k,j,i) = 0.0
        END IF
      END DO
    END DO
  END DO
END IF

!     1.2 INITIALISE TEMPORARY ARRAYS
!     -------------------------------

CALL zero33(phia,nz,nr,nc)
CALL zero33(phib,nz,nr,nc)

!     2. CORRECTOR TERMS
!     ------------------

IF (jadvhs > 0) THEN
  
!     2.1 X-DIRECTION
!     ---------------
  
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
          ucorr(k,j,i) = 0.0
          IF (npix(j,i) == 1) THEN
            ucorr(k,j,i) = ucorr(k,j,i) - del3*phi(k,j,i)*  &
                u2f(k,j,i)*gz2u(k,j,i)/(gx2(j,i)*gz2(k,j,i))
          ELSE IF (npix(j,i) > 1) THEN
            ucorr(k,j,i) = ucorr(k,j,i) - del3*phi(k,j,i)* u2f(k,j,i)/gx2(j,i)
          END IF
          IF (npix(j,i+1) == 1) THEN
            ucorr(k,j,i) = ucorr(k,j,i) + del3*phi(k,j,i)*  &
                u2f(k,j,i+1)*gz2u(k,j,i+1)/(gx2(j,i)*gz2(k,j,i))
          ELSE IF (npix(j,i+1) > 1) THEN
             ucorr(k,j,i) = ucorr(k,j,i) + del3*phi(k,j,i)*  &
                u2f(k,j,i+1)/gx2(j,i)
          END IF
        END DO
      END IF
    END DO
  END DO
  211   CONTINUE
!     2.2 Y-DIRECTION
!     ---------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
          vcorr(k,j,i) = 0.0
          IF (npiy(j,i) == 1) THEN
            vcorr(k,j,i) = vcorr(k,j,i) - del3*phi(k,j,i)*v2f(k,j,i)*  &
                gz2v(k,j,i)*cosphiv(j)/(gy2(j)*gz2(k,j,i)*cosphi(j))
          ELSE IF (npiy(j,i) > 1) THEN
            vcorr(k,j,i) = vcorr(k,j,i)- del3*phi(k,j,i)*v2f(k,j,i)*  &
                cosphiv(j)/(gy2(j)*cosphi(j))
          END IF
          IF (npiy(j+1,i) == 1) THEN
            vcorr(k,j,i) = vcorr(k,j,i) + del3*phi(k,j,i)*v2f(k,j+1,i)*  &
                gz2v(k,j+1,i)*cosphiv(j+1)/(gy2(j)*gz2(k,j,i)*cosphi(j))
          ELSE IF (npiy(j+1,i) > 1) THEN
            vcorr(k,j,i) = vcorr(k,j,i) + del3*phi(k,j,i)*v2f(k,j+1,i)*  &
                cosphiv(j+1)/(gy2(j)*cosphi(j))
          END IF
        END DO
      END IF
    END DO
  END DO
  221   CONTINUE
  
  
!     2.3 Z-DIRECTION
!     ---------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
          wcorr(k,j,i) = del3*phi(k,j,i)* (w2(k+1,j,i)-w2(k,j,i))/gz2(k,j,i)
        END DO
      END IF
    END DO
  END DO
  231   CONTINUE
  
END IF


!     3. UPDATE PHI WITHOUT OPERATOR SPLITTING
!     ----------------------------------------

IF (jadvhs <= 1) THEN
  
  
!     3.1 ASSEMBLE INERTIA TERMS
!     --------------------------
  
  DO  i=1,nc
    DO  j=1,nr
      DO  k=1,nz
        plhsa(k,j,i) = 0.0
        plhsb(k,j,i) = 1.0
        plhsc(k,j,i) = 0.0
        IF (nwd(j,i) == 1) THEN
          IF (jadvhs == 1) THEN
            prhs(k,j,i) = phi(k,j,i)
          ELSE
            prhs(k,j,i) = phi(k,j,i)*gz1(k,j,i)/gz2(k,j,i)
          END IF
        ELSE
          prhs(k,j,i) = 0.0
        END IF
      END DO
    END DO
  END DO


  
!     3.2 ADD CORRECTOR TERMS
!     -----------------------
  
  IF (jadvhs == 1) THEN
    DO  i=1,nc
      DO  j=1,nr
        IF (nwd(j,i) == 1) THEN
          DO  k=1,nz
            prhs(k,j,i) = prhs(k,j,i) + ucorr(k,j,i)  &
                + vcorr(k,j,i) + wcorr(k,j,i)
          END DO
        END IF
      END DO
    END DO
    321   CONTINUE
  END IF


!     3.3 ADD HORIZONTAL ADVECTIVE/DIFFUSIVE TERMS
!     --------------------------------------------
!
  
  CALL xadvdis(phi,heddydu,u2f,gz2,nz,ivpu, phivp,jadvhs,jdifs,plhsb,prhs)
  CALL yadvdis(phi,heddydv,v2f,gz2,nz,ivpv, phivp,jadvhs,jdifs,plhsb,prhs)
  
  
!     3.4 ADD SOURCE/SINK TERMS
!     -------------------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz
          prhs(k,j,i)  = prhs(k,j,i)  + del3*source(k,j,i)
          plhsb(k,j,i) = plhsb(k,j,i) + del3*sink(k,j,i)
        END DO
      END IF
    END DO
  END DO
  341   CONTINUE

  
!     3.5 ADD VERTICAL ADVECTION/DIFFUSION TERMS
!     ------------------------------------------
  
  CALL zadvdis(phi, vedphi, w2phi, gz2, nwd,isur,ibot,  &
      nz,nr,nc,sflux,stran,sphi,bflux,btran,bphi,  &
      jadvws, plhsa, plhsb, plhsc, prhs)
  
!     3.6 INVERT TRIDIAGONAL SYSTEM TO OBTAIN PHI
!     -------------------------------------------
  
  CALL thomv(plhsa,plhsb,plhsc,prhs,phi,nz,nr,nc)
  
  GO TO 3000
  
END IF

!     4. EXPLICIT X-DIRECTION ADVECTION AND DIFFUSION - STEP A
!     --------------------------------------------------------

!     4.1 ASSEMBLE INERTIA TERMS - STEP A
!     -----------------------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        plhsb(k,j,i) = 1.0
        prhs(k,j,i) = phi(k,j,i) + ucorr(k,j,i)
      END DO
    END IF
  END DO
END DO
411   CONTINUE


!     4.2 ASSEMBLE X-ADVECTION/DIFFUSION TERMS - STEP A
!     -------------------------------------------------


CALL xadvdis(phi,heddydu,u2f,gz2,nz,ivpu, phivp,jadvhs,jdifs,plhsb,prhs)


!     4.3 UPDATE TO GET PHIA-X
!     ------------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        phia(k,j,i) = prhs(k,j,i)/plhsb(k,j,i)
      END DO
    END IF
  END DO
END DO
431   CONTINUE

!     5. EXPLICIT Y-DIRECTION ADVECTION AND DIFFUSION - STEP A
!     --------------------------------------------------------

!     5.1 ASSEMBLE INERTIA TERMS - STEP A
!     -----------------------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        plhsb(k,j,i) = 1.0
        prhs(k,j,i) = phia(k,j,i) + vcorr(k,j,i)
      END DO
    END IF
  END DO
END DO
511   CONTINUE


!     5.2 ASSEMBLE Y-ADVECTION/DIFFUSION TERMS - STEP A
!     -------------------------------------------------


CALL yadvdis(phia,heddydv,v2f,gz2,nz,ivpv, phivp,jadvhs,jdifs,plhsb,prhs)


!     5.3 UPDATE TO GET PHIA-XY
!     -------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        phia(k,j,i) = prhs(k,j,i)/plhsb(k,j,i)
      END DO
    END IF
  END DO
END DO
531   CONTINUE

!     6. IMPLICIT Z-DIRECTION ADVECTION AND DIFFUSION  - STEP A
!     ---------------------------------------------------------

!     6.1 PREPARE ARGUMENTS FOR CALL TO ZADVDIS - STEP A
!     --------------------------------------------------


DO  i=1,nc
  DO  j=1,nr
    DO  k=1,nz
      plhsa(k,j,i)= 0.0
      plhsb(k,j,i)= 1.0
      plhsc(k,j,i)= 0.0
      IF (nwd(j,i) == 1) THEN
        prhs(k,j,i) = phia(k,j,i) + wcorr(k,j,i)
      ELSE
        prhs(k,j,i) = 0.0
      END IF
    END DO
  END DO
END DO


!     6.2 ADD SOURCE/SINK TERMS - STEP A
!     ----------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        prhs(k,j,i)  =  prhs(k,j,i) + del3*source(k,j,i)
        plhsb(k,j,i) = plhsb(k,j,i) + del3*sink(k,j,i)
      END DO
    END IF
  END DO
END DO
621   CONTINUE


!     6.3 EVALUATE VERTICAL ADVECTION AND DIFFUSION - STEP A
!     ------------------------------------------------------

!  (must be performed AFTER SOURCE/SINK since Dirichlet bcs
!   overwrite PRHS,etc.)

CALL zadvdis(phia, vedphi, w2phi, gz2, nwd,isur,ibot,  &
    nz,nr,nc,sflux,stran,sphi,bflux,btran,bphi, jadvws, plhsa, plhsb, plhsc, prhs)


!     6.4 INVERT TRIDIAGONAL SYSTEM TO OBTAIN PHIA
!     --------------------------------------------


CALL thomv(plhsa,plhsb,plhsc,prhs,phia,nz,nr,nc)

!     7. IMPLICIT Z-DIRECTION ADVECTION AND DIFFUSION - STEP B
!     --------------------------------------------------------

!     7.1 ASSEMBLE INERTIA TERMS - STEP B
!     -----------------------------------

DO  i=1,nc
  DO  j=1,nr
    DO  k=1,nz
      plhsa(k,j,i)= 0.0
      plhsb(k,j,i)= 1.0
      plhsc(k,j,i)= 0.0
      IF (nwd(j,i) == 1) THEN
        prhs(k,j,i) = phi(k,j,i) + wcorr(k,j,i)
      ELSE
        prhs(k,j,i) = 0.0
      END IF
    END DO
  END DO
END DO


!     7.2 ADD SOURCE/SINK TERMS - STEP B
!     ----------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        prhs(k,j,i)  = prhs(k,j,i)  + del3*source(k,j,i)
        plhsb(k,j,i) = plhsb(k,j,i) + del3*sink(k,j,i)
      END DO
    END IF
  END DO
END DO
721   CONTINUE


!     7.3 EVALUATE VERTICAL ADVECTION AND DIFFUSION - STEP B
!     ------------------------------------------------------

!  (must be performed AFTER SOURCE/SINK since Dirichlet bcs
!   overwrite PRHS,etc.)

CALL zadvdis(phi, vedphi, w2phi, gz2, nwd,isur,ibot,  &
    nz,nr,nc,sflux,stran,sphi,bflux,btran,bphi, jadvws, plhsa, plhsb, plhsc, prhs)


!     7.4 INVERT TRIDIAGONAL SYSTEM TO OBTAIN PHIB
!     --------------------------------------------

CALL thomv(plhsa,plhsb,plhsc,prhs,phib,nz,nr,nc)

!     8. EXPLICIT Y-DIRECTION ADVECTION AND DIFFUSION - STEP B
!     --------------------------------------------------------

!     8.1 ASSEMBLE INERTIA TERMS - STEP B
!     -----------------------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        plhsb(k,j,i) = 1.0
        prhs(k,j,i) = phib(k,j,i) + vcorr(k,j,i)
      END DO
    END IF
  END DO
END DO
811   CONTINUE


!     8.2 ASSEMBLE ADVECTION/DIFFUSION TERMS - STEP B
!     -----------------------------------------------


CALL yadvdis(phib,heddydv,v2f,gz2,nz,ivpv, phivp,jadvhs,jdifs,plhsb,prhs)


!     8.3 UPDATE TO GET PHIB-XY
!     -------------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        phib(k,j,i) = prhs(k,j,i)/plhsb(k,j,i)
      END DO
    END IF
  END DO
END DO
831   CONTINUE

!     9. EXPLICIT X-DIRECTION ADVECTION AND DIFFUSION - STEP B
!     --------------------------------------------------------

!     9.1 ASSEMBLE INERTIA TERMS - STEP B
!     -----------------------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        plhsb(k,j,i) = 1.0
        prhs(k,j,i) = phib(k,j,i) + ucorr(k,j,i)
      END DO
    END IF
  END DO
END DO
911   CONTINUE


!     9.2 ASSEMBLE X-ADVECTION/DIFFUSION TERMS - STEP B
!     -------------------------------------------------


CALL xadvdis(phib,heddydu,u2f,gz2,nz,ivpu, phivp,jadvhs,jdifs,plhsb,prhs)


!     9.3 UPDATE TO GET PHIB
!     ----------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        phib(k,j,i) = prhs(k,j,i)/plhsb(k,j,i)
      END DO
    END IF
  END DO
END DO
931   CONTINUE


!     10. TAKE AVERAGE OF PHIA AND PHIB FOR SYMMETRICALLY SPLIT OPERATOR
!     ------------------------------------------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz
        phi(k,j,i) = 0.5*(phia(k,j,i)+phib(k,j,i))
      END DO
    END IF
  END DO
END DO

3000 CONTINUE

RETURN

END SUBROUTINE transpc

!===========================================================================

SUBROUTINE transpw(phi,vedw,source,sink,  &
    call1,isur,ibot,jadvhs,jadvws,jdifs, ivpu,ivpv,phivp,  &
    sflux,stran,sphi, bflux,btran,bphi)
!***********************************************************************

!    *TRANSPW*  SOLVE TRANSPORT EQUATION FOR A SCALAR QUANTITY AT W-NODES

!     AUTHOR - PATRICK LUYTEN AND KEVIN RUDDICK

!     LAST UPDATE - 26 Mar 1999         @(COHERENS)transp.f 8.4

!     DESCRIPTION - SOLVES ADVECTION-DIFFUSION EQUATION FOR
!                   QUANTITIES STORED AT W NODES (TURBULENT VARIABLES)
!                 - SOURCE TERMS ARE TAKEN EXPLICITLY, SINK TERMS
!                   QUASI-IMPLCITLY FOLLOWING PATANKAR (1980)
!                 - HORIZONTAL ADVECTION/DIFFUSION IS TAKEN EXPLICITLY
!                 - VERTICAL ADVECTION IS SEMI-IMPLICIT VERTICAL DIFFUSION
!                   FULLY IMPLICIT
!                 - ADVECTION SCHEME SELECTED BY IADVS
!                      = 1 => FIRST ORDER UPWIND
!                      = 2 => LAX-WENDROFF/CENTRAL
!                      = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                      = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!                 - HORIZONTAL DIFFUSION SELECTED BY IODIF (IAHDHT=1)
!                      = 0 => OFF
!                      = 1 => UNIFORM
!                      = 2 => SMAGORINSKY
!                 - USES CORRECTOR TERMS TO ELIMINATE THE JACOBIAN FROM
!                   THE TIME DERIVATIVE
!                 - USES THE "FILTERED" ADVECTIVE VELOCITES U2F, V2F
!                 - TIME INTEGRATION IS PERFORMED WITH OR WITHOUT OPERATOR
!                   SPLITTING DEPENDING ON THE TYPE OF ADVECTION SCHEME
!                 - A NUMBER OF ARRAYS NEED TO BE INTERPOLATED
!                 - THE INDEX OF CELL-CENTERED QUNATITIES IS INCREASED BY 1
!                   FOR CONSISTENCY

!     REFERENCE - Sections III-4.5 and V-1.8 of the User Documentation

!     CALLING PROGRAM - DISSIP, TKLENG, TURBEN

!     EXTERNALS - ERROR, THOMV, XADVDIS, YADVDIS, ZADVDIS, ZERO

!***********************************************************************
!*    ARGUMENTS


REAL, INTENT(IN OUT)                        :: phi(nz+1,nr,nc)
REAL, INTENT(IN)                         :: vedw(nz+1,nr,nc)
REAL, INTENT(IN)                         :: source(nz+1,nr,nc)
REAL, INTENT(IN)                         :: sink(nz+1,nr,nc)
LOGICAL, INTENT(IN)                      :: call1
INTEGER, INTENT(IN)                  :: isur
INTEGER, INTENT(IN)                  :: ibot
INTEGER, INTENT(IN)                  :: jadvhs
INTEGER, INTENT(IN)                  :: jadvws
INTEGER, INTENT(IN)                  :: jdifs
INTEGER, INTENT(IN)                  :: ivpu(0:nobu)
INTEGER, INTENT(IN)                  :: ivpv(0:nobv)
REAL, INTENT(IN)                     :: phivp(0:nz+1,0:nvprof)
REAL, INTENT(IN)                     :: sflux(nr,nc)
REAL, INTENT(IN)                     :: stran(nr,nc)
REAL, INTENT(IN)                     :: sphi(nr,nc)
REAL, INTENT(IN)                     :: bflux(nr,nc)
REAL, INTENT(IN)                     :: btran(nr,nc)
REAL, INTENT(IN)                     :: bphi(nr,nc)

!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL (IN)
!    *IBOT*     INTEGER   BOTTOM B.C.  (SEE ZADVDIS)                  (IN)
!    *ISUR*     INTEGER   SURFACE B.C. (SEE ZADVDIS)                  (IN)
!    *JADVHS*   INTEGER   ADVECTION SCHEME FOR HORIZONTAL/VERTICAL
!                         ADVECTION                                   (IN)
!                         = 0 => ADVECTION DISABLED
!                         = 1 => FIRST ORDER UPWIND
!                         = 2 => LAX-WENDROFF
!                         = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                         = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!    *JADVWS*   INTEGER   ADVECTION SCHEME FOR VERTICAL ADVECTION     (IN)
!                         = 0 => VERTICAL ADVECTION DISABLED
!                         = 1 => FIRST ORDER UPWIND
!                         = 2 => CENTRAL
!                         = 3 => TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                         = 4 => TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!    *JDIFS*    INTEGER   SWITCH TO SELECT SCHEME FOR HORIZONTAL DIFFUSION
!                         = 0 => HORIZONTAL DIFFUSION DISABLED        (IN)
!                         = 1 => UNIFORM DIFFSUSION COEFFICIENTS
!                         = 2 => "SMAGORINSKY" FORMULATION
!    *IVPU*     INTEGER   PROFILE NUMBERS AT U-O.B.                   (IN)
!    *IVPV*     INTEGER   PROFILE NUMBERS AT V-O.B.                   (IN)
!    *BFLUX*    REAL      (SEE ZADVDIS)                  [PHI*m/s]    (IN)
!    *BPHI*     REAL      BOTTOM VALUE OF PHI (SEE ZADVDIS)  [PHI]    (IN)
!    *BTRAN*    REAL      BOTTOM TRANSFER COEFF (SEE ZADVDIS)  [m/s]  (IN)
!    *PHI*      REAL      QUANTITY TO BE ADVECTED/DIFFUSED [PHI]  (IN/OUT)
!    *PHIVP*    REAL      OPEN BOUNDARY CONDITION [PHI]               (IN)
!                       PHIVP(K,*)>PHIVP(0,*)  : IMPOSED INFLOW
!                       PHIVP(K,*)<= PHIVP(0,*): NORMAL GRADIENT CONDITION
!    *SFLUX*    REAL      (SEE ZADVDIS)                   [PHI*m/s]   (IN)
!    *SINK*     REAL      SINK TERMS IN TRANSPORT EQUATION            (IN)
!    *SOURCE*   REAL      SOURCE TERMS IN TRANSPORT EQUATION          (IN)
!    *SPHI*     REAL      SURFACE VALUE OF PHI (SEE ZADVDIS) [PHI]    (IN)
!    *STRAN*    REAL      SURFACE TRANSFER COEFF (SEE ZADVDIS) [m/s]  (IN)
!    *VEDW*     REAL      VERTICAL DIFFUSION COEFFICIENT AT W NODES [m2/s]
!                                                                     (IN)

!-----------------------------------------------------------------------

!*    LOCAL VARIABLES

INTEGER :: i, j, k
REAL :: gz2phi(nz+1,nr,nc)
REAL :: hedphiu(nz+1,nr,nc+1), hedphiv(nz+1,nr+1,nc)
REAL :: phia(nz+1,nr,nc), phib(nz+1,nr,nc)
REAL :: plhsa(nz+1,nr,nc), plhsb(nz+1,nr,nc)
REAL :: plhsc(nz+1,nr,nc), prhs(nz+1,nr,nc)
REAL :: ucorr(nz+1,nr,nc), u2phi(nz+1,nr,nc+1)
REAL :: vcorr(nz+1,nr,nc), v2phi(nz+1,nr+1,nc)
REAL :: vedphi(nz+1,nr,nc)
REAL :: wcorr(nz+1,nr,nc), w2phi(nz+1,nr,nc)
!REAL :: gz1w, gz2ux, gz2vx, gz2w, w2atc

!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *GZ2PHI*   REAL      VERTICAL GRID SPACING ABOUT W NODE [m]
!    *HEDPHIU*  REAL      HORIZONTAL DIFFUSION COEFFICIENT AT X-NODE [m2/s]
!    *HEDPHIV*  REAL      HORIZONTAL DIFFUSION COEFFICIENT AT Y-NODE [m2/s]
!    *PHIA*     REAL      FIRST FRACTIONAL STEP NEW PHI
!    *PHIB*     REAL      SECOND FRACTIONAL STEP NEW PHI
!    *PLHSA*    REAL      LEFT HAND SIDE OF PHI EQUATION /PHI(I,J,K-1)
!    *PLHSB*    REAL      LEFT HAND SIDE OF PHI EQUATION /PHI(I,J,K)
!    *PLHSC*    REAL      LEFT HAND SIDE OF PHI EQUATION /PHI(I,J,K+1)
!    *PRHS*     REAL      RIGHT HAND SIDE OF PHI EQUATION
!    *UCORR*    REAL      CORRECTOR TERMS IN X-DIRECTION
!    *U2PHI*    REAL      FILTERED X-ADVECTION VELOCITY AT X-NODE [m/s]
!    *VCORR*    REAL      CORRECTOR TERMS IN Y-DIRECTION
!    *V2PHI*    REAL      FILTERED Y-ADVECTION VELOCITY SOUTH OF Y-NODE [m/s]
!    *VEDPHI*   REAL      VERTICAL DIFFUSION COEFFICIENT AT THE CENTRE (BELOW PHI)
!    *WCORR*    REAL      CORRECTOR TERMS IN Z-DIRECTION
!    *W2PHI*    REAL      VERTICAL ADVECTIVE VELOCITY AT THE CENTRE (BELOW PHI)

!-----------------------------------------------------------------------

!     1. INITIALISE AND CHECKING
!     --------------------------

!     1.1 ERROR CHECKING
!     ------------------

IF (call1) THEN
  IF ((isur == 0).OR.(ibot == 0)) THEN
    nerrs = 1
    WRITE (6,'(A)') 'TRANSPW not yet implemented with '//  &
        'Neumann boundary conditions'
  END IF
  
  
!     1.2 INITIALISE PHI AT BLOCKED POINTS ON FIRST CALL
!     --------------------------------------------------
  
  
  DO  i=1,nc
    DO  j=1,nr
      DO  k=1,nz+1
        IF (nwd(j,i) == 0) THEN
          phi(k,j,i) = 0.0
        END IF
      END DO
    END DO
  END DO
  
END IF


!     1.3 INITIALISE TEMPORARY ARRAYS
!     -------------------------------


CALL zero33(phia,nz+1,nr,nc)
CALL zero33(phib,nz+1,nr,nc)


!     2. ADVECTIVE VELOCITIES, DIFFUSION COEFFICIENT AND GRID SPACING
!     ---------------------------------------------------------------


!     2.1 X-DIRECTION
!     ---------------

IF ((jadvhs > 0).OR.(jdifs > 0)) THEN
  DO  i=1,nc+1
    DO  j=1,nr
      u2phi(1,j,i) = u2f(1,j,i)
      hedphiu(1,j,i) = heddydu(1,j,i)
      DO  k=2,nz
        u2phi(k,j,i) = 0.5*(u2f(k,j,i)+u2f(k-1,j,i))
        hedphiu(k,j,i) = 0.5*(heddydu(k-1,j,i)+heddydu(k,j,i))
      END DO
      u2phi(nz+1,j,i) = u2f(nz,j,i)
      hedphiu(nz+1,j,i) = heddydu(nz,j,i)
    END DO
  END DO
END IF


!     2.2 Y-DIRECTION
!     ---------------


IF ((jadvhs > 0).OR.(jdifs > 0)) THEN
  DO  i=1,nc
    DO  j=1,nr+1
      v2phi(1,j,i) = v2f(1,j,i)
      hedphiv(1,j,i) = heddydv(1,j,i)
      DO  k=2,nz
        v2phi(k,j,i) = 0.5*(v2f(k,j,i)+v2f(k-1,j,i))
        hedphiv(k,j,i) = 0.5*(heddydv(k-1,j,i)+heddydv(k,j,i))
      END DO
      v2phi(nz+1,j,i) = v2f(nz,j,i)
      hedphiv(nz+1,j,i) = heddydv(nz,j,i)
    END DO
  END DO
END IF


!     2.3 Z-DIRECTION
!     ---------------

!     ---vertical velocity
IF (jadvws > 0) THEN
  DO  i=1,nc
    DO  j=1,nr
      w2phi(1,j,i)  = 0.0
      DO  k=2,nz+1
        w2phi(k,j,i) = w2atc(k-1,j,i)
      END DO
    END DO
  END DO
  231   CONTINUE
ELSE
  CALL zero33(w2phi,nz+1,nr,nc)
END IF

!     ---vertical diffusion coefficient
DO  i=1,nc
  DO  j=1,nr
    vedphi(1,j,i) = vedw(1,j,i)
    DO  k=2,nz+1
      vedphi(k,j,i) = 0.5*(vedw(k-1,j,i)+vedw(k,j,i))
    END DO
  END DO
END DO
232   CONTINUE

!     ---vertical grid spacing
DO  i=1,nc
  DO  j=1,nr
    gz2phi(1,j,i) = gz2(1,j,i)
    DO  k=2,nz
      gz2phi(k,j,i) = gz2w(k,j,i)
    END DO
    gz2phi(nz+1,j,i) = gz2(nz,j,i)
  END DO
END DO
233   CONTINUE


!     3. CORRECTOR TERMS
!     ------------------

IF (jadvhs > 0) THEN
  
  
!     3.1 X-DIRECTION
!     ---------------
  
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=2,nz
          ucorr(k,j,i) = 0.0
          IF (npix(j,i) == 1) THEN
            ucorr(k,j,i) = ucorr(k,j,i) - del3*phi(k,j,i)*  &
                u2phi(k,j,i)*gz2ux(k,j,i)/(gx2(j,i)*gz2w(k,j,i))
          ELSE IF (npix(j,i) > 1) THEN
            ucorr(k,j,i) = ucorr(k,j,i) - del3*phi(k,j,i)*  &
                u2phi(k,j,i)/gx2(j,i)
          END IF
          IF (npix(j,i+1) == 1) THEN
            ucorr(k,j,i) = ucorr(k,j,i) + del3*phi(k,j,i)*  &
                u2phi(k,j,i+1)*gz2ux(k,j,i+1)/(gx2(j,i)*gz2w(k,j,i))
          ELSE IF (npix(j,i+1) > 1) THEN
            ucorr(k,j,i) = ucorr(k,j,i) + del3*phi(k,j,i)*  &
                u2phi(k,j,i+1)/gx2(j,i)
          END IF
        END DO
      END IF
    END DO
  END DO
  311   CONTINUE
  
  
!     3.2 Y-DIRECTION
!     ---------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=2,nz
          vcorr(k,j,i) = 0.0
          IF (npiy(j,i) == 1) THEN
            vcorr(k,j,i) = vcorr(k,j,i) - del3*phi(k,j,i)*v2phi(k,j,i)*  &
                gz2vx(k,j,i)*cosphiv(j)/(gy2(j)*gz2w(k,j,i)*cosphi(j))
          ELSE IF (npiy(j,i) > 1) THEN
            vcorr(k,j,i) = vcorr(k,j,i) - del3*phi(k,j,i)*v2phi(k,j,i)*  &
                cosphiv(j)/(gy2(j)*cosphi(j))
          END IF
          IF (npiy(j+1,i) == 1) THEN
            vcorr(k,j,i) = vcorr(k,j,i) + del3*phi(k,j,i)*v2phi(k,j+1,i)*  &
                gz2vx(k,j+1,i)*cosphiv(j+1) /(gy2(j)*gz2w(k,j,i)*cosphi(j))
          ELSE IF (npiy(j+1,i) > 1) THEN
            vcorr(k,j,i) = vcorr(k,j,i) + del3*phi(k,j,i)*v2phi(k,j+1,i)*  &
                cosphiv(j+1)/(gy2(j)*cosphi(j))
          END IF
        END DO
      END IF
    END DO
  END DO
  321   CONTINUE
  
  
!     3.3 Z-DIRECTION
!     ---------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=2,nz
          wcorr(k,j,i) = del3*phi(k,j,i)*  &
              (w2phi(k+1,j,i)-w2phi(k,j,i))/gz2w(k,j,i)
        END DO
      END IF
    END DO
  END DO
  331   CONTINUE
  
END IF


!     4. UPDATE PHI WITHOUT OPERATOR SPLITTING
!     ----------------------------------------

IF (jadvhs <= 1) THEN
  
  
!     4.1 ASSEMBLE INERTIA TERMS
!     --------------------------
  
  
  DO  i=1,nc
    DO  j=1,nr
      plhsa(1,j,i) = 0.0
      plhsb(1,j,i) = 1.0
      plhsc(1,j,i) = 0.0
      prhs(1,j,i) = 0.0
      DO  k=2,nz
        plhsa(k,j,i) = 0.0
        plhsb(k,j,i) = 1.0
        plhsc(k,j,i) = 0.0
        IF (nwd(j,i) == 1) THEN
          IF (jadvhs == 1) THEN
            prhs(k,j,i) = phi(k,j,i)
          ELSE
            prhs(k,j,i) = phi(k,j,i)*gz1w(k,j,i)/gz2w(k,j,i)
          END IF
        ELSE
          prhs(k,j,i) = 0.0
        END IF
      END DO
      plhsa(nz+1,j,i) = 0.0
      plhsb(nz+1,j,i) = 1.0
      plhsc(nz+1,j,i) = 0.0
      prhs(nz+1,j,i) = 0.0
    END DO
  END DO
  411   CONTINUE
  
  
!     4.2 ADD CORRECTOR TERMS
!     -----------------------
  
  IF (jadvhs == 1) THEN
    DO  i=1,nc
      DO  j=1,nr
        IF (nwd(j,i) == 1) THEN
          DO  k=2,nz
            prhs(k,j,i) = prhs(k,j,i) + ucorr(k,j,i)  &
                + vcorr(k,j,i) + wcorr(k,j,i)
          END DO
        END IF
      END DO
    END DO
    421   CONTINUE
  END IF
  
  
!     4.3 ADD HORIZONTAL ADVECTION/DIFFUSION TERMS
!     --------------------------------------------
  
  
  CALL xadvdis(phi,hedphiu,u2phi,gz2phi,nz+1,ivpu,  &
      phivp,jadvhs,jdifs,plhsb,prhs)
  CALL yadvdis(phia,hedphiv,v2phi,gz2phi,nz+1,ivpv,  &
      phivp,jadvhs,jdifs,plhsb,prhs)
  
  
!     4.4 ADD SOURCE/SINK TERMS
!     -------------------------
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=2,nz
          prhs(k,j,i)  = prhs(k,j,i)  + del3*source(k,j,i)
          plhsb(k,j,i) = plhsb(k,j,i) + del3*sink(k,j,i)
        END DO
      END IF
    END DO
  END DO
  441   CONTINUE
  
  
!     4.5 ADD VERTICAL ADVECTION/DIFFUSION TERMS
!     ------------------------------------------
  
  
  CALL zadvdis(phia, vedphi, w2phi, gz2phi, nwd,isur,ibot,  &
      nz+1,nr,nc,sflux,stran,sphi,bflux,btran,bphi,  &
      jadvws, plhsa, plhsb, plhsc, prhs)
  
  
!     4.6 INVERT TRIDIAGONAL SYSTEM TO OBTAIN PHI
!     -------------------------------------------
  
  
  CALL thomv(plhsa,plhsb,plhsc,prhs,phi,nz+1,nr,nc)
  
  GO TO 3000
  
END IF


!     5. EXPLICIT X-DIRECTION ADVECTION AND DIFFUSION - STEP A
!     --------------------------------------------------------

!     5.1 ASSEMBLE INERTIA TERMS - STEP A
!     -----------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        plhsb(k,j,i) = 1.0
        prhs(k,j,i) = phi(k,j,i) + ucorr(k,j,i)
      END DO
    END IF
  END DO
END DO
511   CONTINUE


!     5.2 ASSEMBLE X-ADVECTION/DIFFUSION TERMS - STEP A
!     -------------------------------------------------


CALL xadvdis(phi,hedphiu,u2phi,gz2phi,nz+1,ivpu,  &
    phivp,jadvhs,jdifs,plhsb,prhs)


!     5.3 UPDATE TO GET PHIA-X
!     ------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        phia(k,j,i) = prhs(k,j,i)/plhsb(k,j,i)
      END DO
    END IF
  END DO
END DO
531   CONTINUE


!     6. EXPLICIT Y-DIRECTION ADVECTION AND DIFFUSION - STEP A
!     --------------------------------------------------------

!     6.1 ASSEMBLE INERTIA TERMS - STEP A
!     -----------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        plhsb(k,j,i) = 1.0
        prhs(k,j,i) = phia(k,j,i) + vcorr(k,j,i)
      END DO
    END IF
  END DO
END DO
611   CONTINUE


!     6.2 ASSEMBLE Y-ADVECTION/DIFFUSION TERMS - STEP A
!     -------------------------------------------------


CALL yadvdis(phia,hedphiv,v2phi,gz2phi,nz+1,ivpv,  &
    phivp,jadvhs,jdifs,plhsb,prhs)


!     6.3 UPDATE TO GET PHIA-XY
!     -------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        phia(k,j,i) = prhs(k,j,i)/plhsb(k,j,i)
      END DO
    END IF
  END DO
END DO
631   CONTINUE


!     7. IMPLICIT Z-DIRECTION ADVECTION AND DIFFUSION  - STEP A
!     ---------------------------------------------------------

!     7.1 PREPARE ARGUMENTS FOR CALL TO ZADVDIS - STEP A
!     --------------------------------------------------


DO  i=1,nc
  DO  j=1,nr
    plhsa(1,j,i)= 0.0
    plhsb(1,j,i)= 1.0
    plhsc(1,j,i)= 0.0
    prhs(1,j,i) = 0.0
    DO  k=2,nz
      plhsa(k,j,i)= 0.0
      plhsb(k,j,i)= 1.0
      plhsc(k,j,i)= 0.0
      IF (nwd(j,i) == 1) THEN
        prhs(k,j,i) = phia(k,j,i) + wcorr(k,j,i)
      ELSE
        prhs(k,j,i) = 0.0
      END IF
    END DO
    plhsa(nz+1,j,i)= 0.0
    plhsb(nz+1,j,i)= 1.0
    plhsc(nz+1,j,i)= 0.0
    prhs(nz+1,j,i) = 0.0
  END DO
END DO
711   CONTINUE


!     7.2 ADD SOURCE/SINK TERMS - STEP A
!     ----------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        prhs(k,j,i)  =  prhs(k,j,i) + del3*source(k,j,i)
        plhsb(k,j,i) = plhsb(k,j,i) + del3*sink(k,j,i)
      END DO
    END IF
  END DO
END DO
721   CONTINUE


!     7.3 EVALUATE VERTICAL ADVECTION AND DIFFUSION - STEP A
!     ------------------------------------------------------

!  (must be performed AFTER SOURCE/SINK since Dirichlet bcs
!   overwrite PRHS,etc.)

CALL zadvdis(phia, vedphi, w2phi, gz2phi, nwd,isur,ibot,  &
    nz+1,nr,nc,sflux,stran,sphi,bflux,btran,bphi,  &
    jadvws, plhsa, plhsb, plhsc, prhs)


!     7.4 INVERT TRIDIAGONAL SYSTEM TO OBTAIN PHIA
!     --------------------------------------------


CALL thomv(plhsa,plhsb,plhsc,prhs,phia,nz+1,nr,nc)


!     8. IMPLICIT Z-DIRECTION ADVECTION AND DIFFUSION - STEP B
!     --------------------------------------------------------

!     8.1 ASSEMBLE INERTIA TERMS - STEP B
!     -----------------------------------


DO  i=1,nc
  DO  j=1,nr
    plhsa(1,j,i)= 0.0
    plhsb(1,j,i)= 1.0
    plhsc(1,j,i)= 0.0
    prhs(1,j,i) = 0.0
    DO  k=2,nz
      plhsa(k,j,i)= 0.0
      plhsb(k,j,i)= 1.0
      plhsc(k,j,i)= 0.0
      IF (nwd(j,i) == 1) THEN
        prhs(k,j,i) = phi(k,j,i) + wcorr(k,j,i)
      ELSE
        prhs(k,j,i) = 0.0
      END IF
    END DO
    plhsa(nz+1,j,i)= 0.0
    plhsb(nz+1,j,i)= 1.0
    plhsc(nz+1,j,i)= 0.0
    prhs(nz+1,j,i) = 0.0
  END DO
END DO
811   CONTINUE


!     8.2 ADD SOURCE/SINK TERMS - STEP B
!     ----------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        prhs(k,j,i)  =  prhs(k,j,i) + del3*source(k,j,i)
        plhsb(k,j,i) = plhsb(k,j,i) + del3*sink(k,j,i)
      END DO
    END IF
  END DO
END DO
821   CONTINUE


!     8.3 EVALUATE VERTICAL ADVECTION AND DIFFUSION - STEP B
!     ------------------------------------------------------

!  (must be performed AFTER SOURCE/SINK since Dirichlet bcs
!   overwrite PRHS,etc.)

CALL zadvdis(phi, vedphi, w2phi, gz2phi, nwd,isur,ibot,  &
    nz+1,nr,nc,sflux,stran,sphi,bflux,btran,bphi,  &
    jadvws, plhsa, plhsb, plhsc, prhs)


!     8.4 INVERT TRIDIAGONAL SYSTEM TO OBTAIN PHIB
!     --------------------------------------------


CALL thomv(plhsa,plhsb,plhsc,prhs,phib,nz+1,nr,nc)


!     9. EXPLICIT Y-DIRECTION ADVECTION AND DIFFUSION - STEP B
!     --------------------------------------------------------

!     9.1 ASSEMBLE INERTIA TERMS - STEP B
!     -----------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        plhsb(k,j,i) = 1.0
        prhs(k,j,i) = phib(k,j,i) + vcorr(k,j,i)
      END DO
    END IF
  END DO
END DO
911   CONTINUE


!     9.2 ASSEMBLE ADVECTION/DIFFUSION TERMS - STEP B
!     -----------------------------------------------


CALL yadvdis(phib,hedphiv,v2phi,gz2phi,nz+1,ivpv,  &
    phivp,jadvhs,jdifs,plhsb,prhs)


!     9.3 UPDATE TO GET PHIB-XY
!     -------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        phib(k,j,i) = prhs(k,j,i)/plhsb(k,j,i)
      END DO
    END IF
  END DO
END DO
931   CONTINUE


!     10. EXPLICIT X-DIRECTION ADVECTION AND DIFFUSION - STEP B
!     --------------------------------------------------------

!     10.1 ASSEMBLE INERTIA TERMS - STEP B
!     -----------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        plhsb(k,j,i) = 1.0
        prhs(k,j,i) = phib(k,j,i) + ucorr(k,j,i)
      END DO
    END IF
  END DO
END DO


!     10.2 ASSEMBLE X-ADVECTION/DIFFUSION TERMS - STEP B
!     -------------------------------------------------


CALL xadvdis(phib,hedphiu,u2phi,gz2phi,nz+1,ivpu,  &
    phivp,jadvhs,jdifs,plhsb,prhs)


!     10.3 UPDATE TO GET PHIB
!     ----------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        phib(k,j,i) = prhs(k,j,i)/plhsb(k,j,i)
      END DO
    END IF
  END DO
END DO


!     11. TAKE AVERAGE OF PHIA AND PHIB FOR SYMMETRICALLY SPLIT OPERATOR
!     ------------------------------------------------------------------

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz+1
        phi(k,j,i) = 0.5*(phia(k,j,i)+phib(k,j,i))
      END DO
    END IF
  END DO
END DO


3000 CONTINUE

RETURN

END SUBROUTINE transpw


! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:44:02

SUBROUTINE turben
!************************************************************************

!    *TURBEN*     SOLVE TRANSPORT EQUATION FOR TURBULENCE ENERGY

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 7 May 1998        @(COHERENS)turben.f 8.4

!       DESCRIPTION - STORE OLD VALUES OF K
!                   - APPLY SURFACE AND BOTTOM BOUNDARY CONDITIONS
!                   - EVALUATE SOURCE/SINK TERMS
!                   - UPDATE TURBULENGE ENERGY BY CALLING TRANSPW
!                   - APPLY (NUMERICAL) LOWER LIMIT TO DIVISION BY ZERO

!       REFERENCE - Section III-1.2.2 and III-4.5.3 of the User Documentation

!       CALLING PROGRAM - VEDDY2

!       EXTERNALS - PRODUC, TRANSPW

!************************************************************************
!*    LOCAL VARIABLES

LOGICAL :: call1, flag
INTEGER :: i, iadvt, ibot, idift, isur, j, k
REAL :: btke(nr,nc), sink(nz+1,nr,nc)
REAL :: source(nz+1,nr,nc), stke(nr,nc)

SAVE call1
DATA call1 /.true./

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS THE FIRST CALL
!    *FLAG*     LOGICAL   .TRUE. IF TRANSPW IS CALLED FOR THE FIRST TIME
!    *IADVT*    INTEGER   TYPE OF SCHEME FOR ADVECTION (0/1/2/3/4)
!    *IBOT*     INTEGER   =1 (DIRICHLET BOTTOM B.C.)
!    *IDIFT*    INTEGER   TYPE OF SCHEME FOR HORIZONTAL DIFFUSION (0/1/2)
!    *ISUR*     INTEGER   =1 (DIRICHLET SURFACE B.C.)
!    *BTKE*     REAL      BOTTOM VALUE OF T.K.E.
!    *SINK*     REAL      SINK TERMS IN T.K.E. EQUATION
!    *SOURCE*   REAL      SOURCE TERMS IN T.K.E. EQUATION
!    *STKE*     REAL      SURFACE VALUE OF T.K.E.

!------------------------------------------------------------------------

!     1. INITIALISE PARAMETERS AND ARRAYS
!     -----------------------------------

!WRITE (*,*) 'TRACE: turben'

isur = 1
ibot = 1
IF (iahdht == 0) THEN
  iadvt = 0
  idift = 0
ELSE
  iadvt = iadvs
  idift = iodif
END IF

flag = call1
call1 = .false.


!     2. STORE OLD VALUES
!     -------------------


DO  i=1,nc
  DO  j=1,nr
    DO  k=1,nz+1
      tkeold(k,j,i) = tkew(k,j,i)
    END DO
  END DO
END DO


!     3. BOUNDARY CONDITIONS (SURFACE, BOTTOM)
!     ----------------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      stke(j,i) = sstot(j,i)/SQRT(cmu)
      btke(j,i) = bstot(j,i)/SQRT(cmu)
    END IF
  END DO
END DO


!     4. SOURCE/SINK TERMS
!     --------------------


CALL produc

DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        IF (tkew(k,j,i) > 0.0) THEN
          source(k,j,i) = shprod(k,j,i)
          sink(k,j,i) = dissw(k,j,i)/tkew(k,j,i)
          IF (buprod(k,j,i) > 0.0) THEN
            source(k,j,i) = source(k,j,i)+buprod(k,j,i)
          ELSE
            sink(k,j,i) = sink(k,j,i)-buprod(k,j,i)/tkew(k,j,i)
          END IF
        ELSE
          source(k,j,i) = 0.0
          sink(k,j,i) = 0.0
        END IF
      END DO
    END IF
  END DO
END DO


!     5. UPDATE T.K.E.
!     ---------------


CALL transpw(tkew,veddyk,source,sink, flag,isur,ibot,iadvt,iadvt,idift,  &
    ivpobu,ivpobv,turvp, zeros1,zeros2,stke,  &
    zeros4,zeros5,btke)


!     6. LIMITING CONDITIONS (NUMERICAL)
!     ----------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=1,nz+1
        IF (tkew(k,j,i) < 1.0E-20) tkew(k,j,i) = 0.0
      END DO
    END IF
  END DO
END DO


RETURN

END SUBROUTINE turben

!=======================================================================

SUBROUTINE produc
!************************************************************************

!    *PRODUC*   EVALUATE THE SHEAR AND BUOYANCY SOURCE/SINK TERMS IN THE
!               T.K.E EQUATION

!     AUTHOR - PATRICK LUYTEN

!     LAST UPDATE  - 7 Dec 1994         @(COHERENS)turben.f 6.1

!     DESCRIPTION - THE SHEAR AND BUOYANCY SOURCE/SINK TERMS ARE STORED
!                   IN SHPROD AND BUPROD

!     REFERENCE - Section III-1.2.2 of the User Documentation

!     CALLING PROGRAM - TURBEN

!     EXTERNALS - BUOFR2, SHFR2

!************************************************************************
!*    LOCAL VARIABLES

INTEGER :: i, j, k
REAL :: z1
!REAL :: buofr2, shfr2



!     1. EVALUATE SHEAR AND BUOYANCY SOURCE/SINK TERMS
!     ------------------------------------------------


DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      DO  k=2,nz
        IF (tkew(k,j,i) > 0.0) THEN
          IF (itcpar == 1) THEN
            z1 = SQRT(tkew(k,j,i))*zlw(k,j,i)
          ELSE
            z1 = (tkew(k,j,i)/dissw(k,j,i))*tkew(k,j,i)
          END IF
          shprod(k,j,i) = smu(k,j,i)*shfr2(k,j,i)*z1
          buprod(k,j,i) = -shb(k,j,i)*buofr2(k,j,i)*z1
        ELSE
          shprod(k,j,i) = 0.0
          buprod(k,j,i) = 0.0
        END IF
      END DO
    END IF
  END DO
END DO

RETURN

END SUBROUTINE produc

!
!=======================================================================

SUBROUTINE czero33(cvar,nx,lenx)
!************************************************************************

!    *CZERO*      INITIALISE A CHARACTER ARRAY

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 12 May 1998 @(COHERENS)Util.f 8.4

!       DESCRIPTION - INSERTS BLANKS INTO THE CHARACTER ARRAY CVAR

!       REFERENCE -

!       CALLING PROGRAM - VARIOUS

!       EXTERNALS - NONE

!************************************************************************

!*    ARGUMENTS


INTEGER, INTENT(IN)                      :: nx
CHARACTER (LEN=*), INTENT(IN OUT)           :: cvar(nx)
INTEGER, INTENT(IN)                      :: lenx

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CVAR*     CHAR      CHARACTER ARRAY TO BE INITIALISED      (IN/OUT)
!    *NX*       INTEGER   DIMENSION OF THE ARRAY                     (IN)
!    *LENX*     INTEGER   LENGTH OF THE ARRAY ELEMENTS               (IN)

!------------------------------------------------------------------------

!*    LOCAL VARIABLES

INTEGER :: l, n


DO  n=1,nx
  DO  l=1,lenx
    cvar(n)(l:l) = ' '
  END DO
END DO


RETURN

END SUBROUTINE czero33

!=======================================================================

SUBROUTINE complx(cr,ci,camp,cphas)
!***********************************************************************

!    *COMPLX*    RETURNS AMPLITUDE AND PHASE OF THE COMPLEX NUMBER (CR,CI)

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 6 Nov 1997 @(COHERENS)Util.f 8.1

!       REFERENCE -

!       DESCRIPTION -

!       CALLING PROGRAM - ANALYS2

!       EXTERNALS -

!***********************************************************************
!*    ARGUMENTS



REAL, INTENT(IN)                        :: cr
REAL, INTENT(IN)                         :: ci
REAL, INTENT(OUT)                        :: camp
REAL, INTENT(OUT)                        :: cphas


!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CR*       REAL      REAL PART OF COMPLEX NUMBER                   (IN)
!    *CI*       REAL      IMAGINARY PART OF COMPLEX NUMBER              (IN)
!    *CAMP*     REAL      AMPLITUDE  OF COMPLEX NUMBER                 (OUT)
!    *CPHAS*    REAL      PHASE OF COMPLEX NUMBER (BETWEEN 0 AND 2*PI) (OUT)

!------------------------------------------------------------------------

camp = SQRT(cr**2+ci**2)
IF (camp > 0.0) THEN
  cphas = ATAN2(ci,cr)
  IF (cphas < 0.0) cphas = cphas + 2.0*pi
ELSE
  cphas = 0.0
END IF


RETURN

END SUBROUTINE complx

!=======================================================================


!========================================================================

SUBROUTINE lsqfit(x,y,ndatmax,ndat,a,b,corr)
!************************************************************************

!    *LSQFIT*     LINEAR REGRESSION ANALYSIS

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 16 Jun 1998 @(COHERENS)Util.f 8.4

!       DESCRIPTION - FITS NDAT DATA POINTS (X,Y) TO THE
!                     STRAIGHT LINE Y=A+B*X
!                   - DETERMINES LINEAR CORRELATION COEFFICIENT

!       REFERENCE - Press W.H., Flannery B.P., Teukolsky S.A. and Vetterling
!                   W.T., 1989. Numerical Recipes. The art of scientific
!                   computing. Cambridge University Press, Cambridge, 702 pp.

!       CALLING PROGRAM - USED IN THE PRINT ROUTINE OF TEST CASE PYCNO

!       EXTERNALS -

!************************************************************************

!*    ARGUMENTS

INTEGER, INTENT(IN)                  :: ndatmax
REAL, INTENT(IN)                         :: x(ndatmax)
REAL, INTENT(IN)                         :: y(ndatmax)
INTEGER, INTENT(IN)                      :: ndat
REAL, INTENT(OUT)                        :: a
REAL, INTENT(OUT)                        :: b
REAL, INTENT(OUT)                        :: corr




!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *X*        REAL      X-COORDINATES OF DATA POINTS                (IN)
!    *Y*        REAL      Y-COORDINATES OF DATA POINTS                (IN)
!    *NDATMAX*  INTEGER   DIMENSION OF ARRAYS XDAT AND YDAT           (IN)
!    *NDAT*     INTEGER   NUMBER OF DATA POINTS                       (IN)
!    *A*        REAL      PARAMETER "A" OF STRAIGHT LINE             (OUT)
!    *B*        REAL      (SLOPE) PARAMETER "B" OF STRAIGHT LINE     (OUT)
!    *CORR*     REAL      (SQUARED) LINEAR CORRELATION COEFFICIENT   (OUT)

!------------------------------------------------------------------------

!*    LOCAL VARIABLES

INTEGER :: n
REAL :: chi2, rdat, sigvar, sxoss, sumt2, sumx, sumy, t


!     ---initialise values
sumx = 0.0
sumy = 0.0
sumt2 = 0.0
b = 0.0

!     ---sum of data values
DO  n=1,ndat
  sumx = sumx + x(n)
  sumy = sumy + y(n)
END DO

!     ---evaluate coefficients B and A
rdat = REAL(ndat)
sxoss = sumx/rdat
DO  n=1,ndat
  t = x(n) - sxoss
  sumt2 = sumt2 + t*t
  b = b + t*y(n)
END DO
b = b/sumt2
a = (sumy-sumx*b)/rdat

!     ---correlation coefficient
chi2 = 0.0
sigvar = 0.0
DO  n=1,ndat
  chi2 = chi2 +(y(n)-a-b*x(n))**2
  sigvar = sigvar + (y(n)-sumy/rdat)**2
END DO
corr = 1.0 - chi2/sigvar


RETURN

END SUBROUTINE lsqfit

!************************************************************************

! Code converted using TO_F90 by Alan Miller
! Date: 2004-05-31  Time: 10:45:17

!    *veddy.f*    EVALUATE EDDY COEFFICIENTS

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 25 May 1999 @(COHERENS)veddy.f 8.4

!       DESCRIPTION - TYPE OF SCHEME SELECTED BY IOPTK
!              = 1 => ALGEBRAIC RELATIONS (FURTHER SELECTIONS BY ITFORM)
!              = 2 => TURBULENCE ENERGY MODEL (FURTHER SELECTIONS BY
!                     NTRANS, ITCPAR, ISTPAR, ILENG, ILIM, IAHDHT)

!       REFERENCE - Section III-1.2 of the User Documentation

!       ROUTINES -  VEDDY1, VEDDY2

!************************************************************************

!========================================================================

SUBROUTINE veddy1
!************************************************************************

!    *VEDDY1*      EVALUATE VERTICAL EDDY VISCOSITY AND DIFFUSIVITY
!                  IF IOPTK=1 (ALGEBRAIC RELATIONS)

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 25 May 1999 @(COHERENS)veddy.f 8.4

!       DESCRIPTION - TURBULENCE SCHEME SELECTED BY ITFORM
!             = 1 => PACANOWSKI-PHILANDER RELATIONS
!             = 2 => MUNK-ANDERSON RELATIONS
!             = 3 => EDDY COEFFICIENTS PROPORTIONAL TO THE MAGNITUDE OF THE
!                    DEPTH-AVERAGED CURRENT TIMES THE DEPTH AND A DAMPING
!                    FUNCTION FOR STRATIFICATION
!             = 4 => EDDY COEFFICIENTS PROPORTIONAL TO THE SQUARED MAGNITUDE
!                    OF THE DEPTH-AVERAGED CURRENT TIMES A DAMPING FUNCTION
!                    FOR STRATIFICATION
!             = 5 => EDDY COEFFICIENTS PROPORTIONAL TO THE MAGNITUDE OF THE
!                    DEPTH-AVERAGED CURRENT TIMES THE DEPTH OF THE
!                    BOTTOM BOUNDARY LAYER AND A DAMPING FUNCTION FOR
!                    STRATIFICATION

!       REFERENCE - Section III-1.2.1 of the User Documentation
!        1) Pacanowski R.C. and Philander S.G.H., 1981. Parameterization of
!           vertical mixing in numerical models of tropical oceans. J. of Phys.
!           Oceanogr., 11, 1443-1451.
!        2) Munk W.H. and Anderson E.R., 1948. Notes on a theory of the
!           thermocline. J. Marine Res., VII(3), 276-295.
!        3) Davies A.M., 1993. A bottom boundary layer-resolving
!           three-dimensional tidal model : A sensitivity study of eddy
!           viscosity formulation. J. of Phys. Oceanogr., 23, 1437-1453.

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS - FVEDMA, GVEDMA, RICH

!************************************************************************
!*    LOCAL VARIABLES

LOGICAL :: call1
INTEGER :: i, j, k
REAL :: phiz(nz+1)
REAL :: curmag, delta, denom, epsmin, vedw, ved0
!REAL :: fvedma, gvedma, h2atc, rich, ud2atc, vd2atc

SAVE call1, phiz
DATA call1 /.true./

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL
!    *PHIZ*     REAL      NORMALISED VERTICAL PROFILE FOR EDDY VISCOSITY
!    *CURMAG*   REAL      SQUARED DEPTH-AVERAGED CURRENT MAGNITUDE [m2/s2]
!    *EPSMIN*   REAL      FACTOR TO PREVENT DIVISION BY ZERO

!------------------------------------------------------------------------

!     1. OUTPUT MESSAGES AND INTIAL ARRAYS ON FIRST CALL
!     --------------------------------------------------


IF (call1) THEN
  
  WRITE (*,*) 'Turbulence model: algebraic'
  IF (itform == 1) THEN
    WRITE (*,*) '  using Pacanowski-Philander relations'
  ELSE IF (itform == 2) THEN
    WRITE (*,*) '  using Munk-Anderson relations'
  ELSE IF (itform == 3) THEN
    WRITE (*,*) '  using U*H law'
  ELSE IF (itform == 4) THEN
    WRITE (*,*) '  using U*U law'
  ELSE IF (itform == 5) THEN
    WRITE (*,*) '  using U*DELTA law'
  END IF
  
!      ---initialise array PHIZ
  IF (itform > 2) THEN
    denom = 1.0+0.5*add1*(adr1-1.0)+0.5*add2*(adr2-1.0)
    DO  k=1,nz+1
      IF (gz0(k,j,i) < add1) THEN
        phiz(k) = (1.0-adr1)*gz0(k,j,i)/add1+adr1
      ELSE IF (gz0(k,j,i) > (1.0-add2)) THEN
        phiz(k) = adr2-(adr2-1.0)*(1.0-gz0(k,j,i))/add2
      ELSE
        phiz(k) = 1.0
      END IF
      phiz(k) = phiz(k)/denom
    END DO
  END IF
  
  call1 =.false.
  
END IF


!     2.EVALUATE EDDY COEFFICIENTS
!     ----------------------------

!     2.1 PACANOWSKI-PHILANDER
!     ------------------------

IF (itform == 1) THEN
  epsmin = rmaxnut**(-1.0/ppn)
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz+1
          denom = MAX(1.0+ppa*rich(k,j,i),epsmin)
          veddyv(k,j,i) = ppved0/(denom**ppn) + ppvisbg
          veddyd(k,j,i) = veddyv(k,j,i)/denom + ppdifbg
        END DO
      END IF
    END DO
  END DO
  211   CONTINUE
  
  
!     2.2 MUNK-ANDERSON RELATIONS
!     ---------------------------
  
ELSE IF (itform == 2) THEN
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=1,nz+1
          veddyv(k,j,i) = amved0*fvedma(k,j,i) + vismol
          veddyd(k,j,i) = amved0*gvedma(k,j,i) + difmol
        END DO
      END IF
    END DO
  END DO
  221   CONTINUE
  
  
!     2.3 U*H-LAW
!     -----------
  
ELSE IF (itform == 3) THEN
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        curmag = ud2atc(j,i)**2 + vd2atc(j,i)**2
        vedw = adlam*SQRT(sstot(j,i))
        DO  k=1,nz+1
          ved0 = adk1*SQRT(curmag)*phiz(k)
          veddyv(k,j,i) = (ved0+vedw)*fvedma(k,j,i) + vismol
          veddyd(k,j,i) = (ved0+vedw)*gvedma(k,j,i) + difmol
        END DO
      END IF
    END DO
  END DO
  231   CONTINUE
  
  
!     2.4 U*U-LAW
!     -----------
  
ELSE IF (itform == 4) THEN
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        curmag = (ud2atc(j,i)/h2atc(j,i))**2 + (vd2atc(j,i)/h2atc(j,i))**2
        vedw = adlam*SQRT(sstot(j,i))
        DO  k=1,nz+1
          ved0 = (adk2*1.0E+04)*curmag*phiz(k)
          veddyv(k,j,i) = (ved0+vedw)*fvedma(k,j,i) + vismol
          veddyd(k,j,i) = (ved0+vedw)*gvedma(k,j,i) + difmol
        END DO
      END IF
    END DO
  END DO
  241   CONTINUE
  
  
!     2.5 U*DELTA-LAW
!     ---------------
  
ELSE IF (itform == 5) THEN
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        curmag = (ud2atc(j,i)/h2atc(j,i))**2 + (vd2atc(j,i)/h2atc(j,i))**2
        delta = MIN(adcnu*(1.0E+04*SQRT(bstot(j,i))),h2atc(j,i))
        vedw = adlam*SQRT(sstot(j,i))
        DO  k=1,nz+1
          ved0 = adk1*SQRT(curmag)*delta*phiz(k)
          veddyv(k,j,i) = (ved0+vedw)*fvedma(k,j,i) + vismol
          veddyd(k,j,i) = (ved0+vedw)*gvedma(k,j,i) + difmol
        END DO
      END IF
    END DO
  END DO
  251   CONTINUE
  
END IF


RETURN

END SUBROUTINE veddy1


!======================================================================

SUBROUTINE veddy2
!************************************************************************

!    *VEDDY2*       EVALUATE EDDY VISCOSITY AND DIFFUSIVITY IF IOPTK=2
!                   (TURBULENCE ENERGY MODELS)

!       AUTHOR - PATRICK LUYTEN

!       LAST UPDATE - 7 May 1998 @(COHERENS)veddy.f 8.4

!       DESCRIPTION - USE LEVEL 2 OR LEVEL 2.5 CLOSURE SCHEME
!                     SELECTED BY THE FOLLOWING SWITCHES :
!         1) NTRANS : SELECTS NUMBER OF TRANSPORT EQUATIONS
!              = 0 => LEVEL 2 SCHEME
!              = 1 => LEVEL 2.5 SCHEME WITH TRANSPORT EQN. FOR T.K.E.
!                    (SOLVED BY TURBEN)
!              = 2 => LEVEL 2.5 SCHEME WITH TRANSPORT EQN. FOR T.K.E.
!                     AND EITHER FOR KL (TKLENG) OR EPS (DISSIP)
!         2) ITCPAR : SELECTS TYPE OF PARAMETERISATION
!              = 1 => K-L OR K-KL
!              = 2 => K-EPS
!         3) ISTPAR : SELECTS FORM OF THE STABILITY FUNCTIONS
!              = 1 => AS FUNCTION OF THE STABILITY PARAMETER (G_h OR A_N)
!              = 2 => AS FUNCTION OF RICHARDSON NUMBER IN ANALOGY WITH
!                     MUNK-ANDERSON RELATIONS
!         4) ILENG : SELECTS LENGTH SCALE PRESCRIPTION WHEN NTRANS=0,1
!              = 1 => PARABOLIC LAW
!              = 2 => "MODIFIED" PARABOLIC LAW
!              = 3 => XING FORMULATION
!              = 4 => BLACKADAR FORMULATION
!         5) ILIM : SELECTS LIMITING CONDITION
!              = 0 => LIMITING CONDITION DISABLED
!              = 1 => LIMITING CONDITION ENABLED
!         6) IAHDHT : SWITCH TO DISABLE/ENABLE ADVECTION AND HORIZONTAL
!                     DIFFUSION OF TURBULENCE
!              = 0 => ADVECTION AND HORIZONTAL DIFFUSION DISABLED
!              = 1 => ADVECTION AND HORIZONTAL DIFFUSION ENABLED

!       REFERENCE - Section III-1.2.2 of the User Documentation
!        1) Mellor G.L. and Yamada T., 1982. Development of a turbulence
!           closure model for geophysical fluid problems. Reviews of Geophys.
!           and Space Physics, 20, 851-875.
!        2) Galperin B., Kantha L.H., Hassid S. and Rosati A., 1988. A
!           quasi-equilibrium turbulent energy model for geophysical flows.
!           J. Atmos. Sc., 45, 55-62.
!        3) Luyten P.J., Deleersnijder E., Ozer J. and Ruddick K.G., 1996.
!           Presentation of a family of turbulence closure models for
!           stratified shallow water flows and preliminary application to the
!           Rhine outflow region. Cont. Shelf Res., 16, 101-130.

!       CALLING PROGRAM - MAINPROG

!       EXTERNALS - BUOFR2, DISLEN, DISMIN, DISSIP, FNSHB, FNSMU,
!                   FNTURF, TKLENG, TLENDIS, TLENG, TURBEN, ZLMAX

!************************************************************************
!*    LOCAL VARIABLES

LOGICAL :: call1
INTEGER :: i, j, k
REAL :: ghan, z1
!REAL :: buofr2, dismin, fnshb, fnsmu, fnturf, zlmax

SAVE call1
DATA call1 /.true./

!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CALL1*    LOGICAL   FLAG TO DETERMINE IF THIS IS THE FIRST CALL
!    *GHAN*     REAL      = G_h (K-L), = ALPHA_N (K-EPS)
!                         GHAN > 0 STABLE ; < 0 UNSTABLE STRATIFICATION

!------------------------------------------------------------------------

!     1. INITIALISATION ON FIRST CALL
!     -------------------------------


IF (call1) THEN
  
  
!     1.1 OUTPUT MESSAGES
!     -------------------
  
!  WRITE (*,*) 'Turbulence model: second order closure'
  IF (ntrans == 0) THEN
    WRITE (*,*) '  level2'
  ELSE IF (ntrans == 1) THEN
    WRITE (*,*) '  one transport equation'
  ELSE IF (ntrans == 2) THEN
    WRITE (*,*) '  two transport equations'
  END IF
  IF (itcpar == 1) THEN
    WRITE (*,*) '  using k-l theory'
  ELSE IF (itcpar == 2) THEN
    WRITE (*,*) '  using K-EPS theory'
  END IF
  IF (istpar == 1) THEN
    WRITE (*,*) '  stability functions depending on stability'// ' parameter'
  ELSE IF (istpar == 2) THEN
    WRITE (*,*) '  stability functions depending on Richardson'// '  number'
  END IF
  IF (ilim == 0) THEN
    WRITE (*,*) '  without limiting conditions'
  ELSE IF (ilim == 1) THEN
    WRITE (*,*) '  with limiting conditions'
  END IF
  IF (iahdht == 0) THEN
    WRITE (*,*) ' advection and horizontal diffusion disabled'
  ELSE IF (iahdht == 1) THEN
    WRITE (*,*) ' advection and horizontal diffusion enabled'
  END IF
  
  call1 = .false.
  
END IF


!     2. UPDATE TURBULENCE VARIABLES
!     ------------------------------

!     2.1 ZERO-EQUATION MODELS
!     -------------------------

IF (ntrans == 0)  THEN
  
!      ---turbulence energy
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        tkew(1,j,i) = bstot(j,i)/SQRT(cmu)
        DO  k=2,nz
          tkew(k,j,i) = fnturf(k,j,i)*zlw(k,j,i)**2
          IF (tkew(k,j,i) < 1.0E-20) tkew(k,j,i) = 0.0
        END DO
        tkew(nz+1,j,i) = sstot(j,i)/SQRT(cmu)
      END IF
    END DO
  END DO
!      ---mixing length
  CALL tleng
!      ---dissipation rate
  CALL dislen
  
  
!     2.2 ONE- AND TWO-EQUATION MODELS
!     --------------------------------
  
  
ELSE IF (ntrans > 0) THEN
  
!      ---turbulence energy
  CALL turben
  
  
!     2.2.1 K-L CASE
!     --------------
  
  IF (itcpar == 1) THEN
!       ---mixing length
    IF (ntrans == 1) THEN
      CALL tleng
    ELSE IF (ntrans == 2) THEN
      CALL tkleng
    END IF
!       ---limiting conditions
    IF (ilim == 1) THEN
      DO  i=1,nc
        DO  j=1,nr
          IF (nwd(j,i) == 1) THEN
            DO  k=2,nz
              zlw(k,j,i) = MIN(zlw(k,j,i),zlmax(k,j,i))
              tkew(k,j,i) = MAX(tkew(k,j,i),tkemin)
            END DO
          END IF
        END DO
      END DO
      221   CONTINUE
    END IF
!       ---dissipation rate
    CALL dislen
    
    
!     2.2.2 K-EPS CASE
!     ----------------
    
  ELSE IF (itcpar == 2) THEN
!       ---dissipation rate
    IF (ntrans == 1) THEN
      CALL tleng
      CALL dislen
    ELSE IF (ntrans == 2) THEN
      CALL dissip
    END IF
!       ---limiting conditions
    IF (ilim == 1) THEN
      DO  i=1,nc
        DO  j=1,nr
          IF (nwd(j,i) == 1) THEN
            DO  k=2,nz
              tkew(k,j,i) = MAX(tkew(k,j,i),tkemin)
              dissw(k,j,i) = MAX(dissw(k,j,i),dismin(k,j,i))
            END DO
          END IF
        END DO
      END DO
      222   CONTINUE
    END IF
!       ---mixing length
    CALL tlendis
    
  END IF
  
END IF


!     2.3 ERROR CHECKING
!     ------------------


nerrs = 0
DO  i=1,nc
  DO  j=1,nr
    IF (nwd(j,i) == 1) THEN
      IF (itcpar == 1) THEN
        DO  k=2,nz
          IF (tkew(k,j,i) > 0.0.AND.zlw(k,j,i) == 0.0) THEN
            nerrs = nerrs + 1
            WRITE (*,*) k, j, i, 'ZLW', zlw(k,j,i)
            WRITE (*,*) '  ZLW must be strictly '//  &
                'positive if TKEW has a non-zero value:',tkew(k,j,i)
          END IF
        END DO
      ELSE IF (itcpar == 2) THEN
        DO  k=2,nz
          IF (tkew(k,j,i) > 0.0.AND.dissw(k,j,i) == 0.0) THEN
            nerrs = nerrs + 1
            WRITE (*,*) k, j, i, 'DISSW', dissw(k,j,i)
            WRITE (*,*) '  DISSW must be strictly '//  &
                'positive if TKEW has a non-zero value:',tkew(k,j,i)
          END IF
        END DO
      END IF
    END IF
  END DO
END DO

!     3. EVALUATE TURBULENT DIFFUSION COEFFICIENTS
!     --------------------------------------------

!     3.1 K-L MODEL
!     -------------

IF (itcpar == 1) THEN
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=2,nz
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
        END DO
      END IF
    END DO
  END DO
  310   CONTINUE
  
  
!     3.2 K-EPS MODEL
!     ---------------
  
ELSE IF (itcpar == 2) THEN
  
  DO  i=1,nc
    DO  j=1,nr
      IF (nwd(j,i) == 1) THEN
        DO  k=2,nz
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
        END DO
      END IF
    END DO
  END DO
  320   CONTINUE
  
END IF


RETURN

8001 FORMAT('Invalid value for element (',i3,',',i3,',',i3,  &
    ') of array ',a,' : ',1PE14.7)

END SUBROUTINE veddy2

SUBROUTINE MC3D
!***********************************************************************
!
!    *MC3D*       UPDATE THE PARTICLE POSITIONS BY ADVECTION AND 
!                 DIFFUSION INSIDE THE DOMAIN
!
!       AUTHOR - KAJO MIRBACH AND PATRICK LUYTEN
!
!       LAST UPDATE - 19 Jan 1998	@(COHERENS)sedlag.f	8.4
!
!       DESCRIPTION - Section V-1.6 of the User Documentation
!
!       REFERENCE - Sections III-4.8.3 and 4.8.4 of the User Documentation
!
!       CALLING PROGRAM - SEDLAG
!
!       EXTERNALS - RANDOM
!
!***********************************************************************
!
!    LOCAL VARIABLES
!
 REAL, PARAMETER :: CLIM = 1.0E-12
 REAL, PARAMETER :: DCLIM = 1.0E-08
 REAL, PARAMETER :: EPSMIN = 1.0E-05
 REAL, PARAMETER :: EPSMAX = 1.0E+05
!

      INTEGER :: IBOT, IEAST, II, INORT, IOLD, IRND, ISOUT, ITOP, IWEST
      INTEGER :: IZAL, JJ, JOLD, KK, KOLD, M, N
      INTEGER :: I,J,K,idum

      REAL :: DTMIN, DTREST, DTX, DTY, DTZ, DUDX, DVDY, DWDZ
      REAL :: DX2, DY2, DZ2, WBOT, WTOP
      REAL :: XADV, XDIFF, XEAST, XEDGE, XEND, XOLD, XRAN, XTMP, XTOT
      REAL :: XU2, XWEST, YADV, YDIFF, YEDGE, YEND, YNORT, YOLD, YSOUT
      REAL :: YTMP, YTOT, YV2, ZADV, ZBOT, ZEDGE, ZEND, ZOLD
      REAL :: ZTMP, ZTOP, ZW2
 !     REAL :: GX2U, GY2V, GZSC, H2ATC, RANDOM
      REAL :: UVTURB,VEDC,ZPOSS, GZZ, POSNEW
      REAL :: dxx,dyy,dist
      REAL :: udiff, vdiff, wdiff
     
!------------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *IBOT*     INTEGER   COUNTER INDICATING WHETHER THE ADVECTED PARTICLE
!                         CROSSED THE BOTTOM CELL BOUNDARY (0/1)
!    *IEAST*    INTEGER   COUNTER INDICATING WHETHER THE ADVECTED PARTICLE
!                         CROSSED THE EASTERN CELL BOUNDARY (0/1)
!    *II*       INTEGER   CURRENT X-INDEX OF PARTICLE
!    *INORT*    INTEGER   COUNTER INDICATING WHETHER THE ADVECTED PARTICLE
!                         CROSSED THE NORTHERN CELL BOUNDARY (0/1)
!    *IOLD*     INTEGER   OLD X-INDEX OF PARTICLE (AT START OF ADVECTION OR 
!                         DIFFUSION)
!    *IRND*     INTEGER   COUNTER INDICATING WHETHER THE ADVECTED PARTICLE
!                         CROSSED ONE OF THE BOX EDGES (0/1/2/3)
!    *ISOUT*    INTEGER   COUNTER INDICATING WHETHER THE ADVECTED PARTICLE
!                         CROSSED THE SOUTHERN CELL BOUNDARY (0/1)
!    *ITOP*     INTEGER   COUNTER INDICATING WHETHER THE ADVECTED PARTICLE
!                         CROSSED THE TOP CELL BOUNDARY (0/1)
!    *IWEST*    INTEGER   COUNTER INDICATING WHETHER THE ADVECTED PARTICLE
!                         CROSSED THE WESTERN CELL BOUNDARY (0/1)
!    *IZAL*     INTEGER   COUNTS NUMBER OF REMAINING TIME STEPS "DTREST" WITHIN
!                         ONE TIME STEP DEL3
!    *JJ*       INTEGER   CURRENT Y-INDEX OF PARTICLE
!    *JOLD*     INTEGER   OLD Y-INDEX OF PARTICLE (AT START OF ADVECTION OR 
!                         DIFFUSION)
!    *KK*       INTEGER   CURRENT Z-INDEX OF PARTICLE
!    *KOLD*     INTEGER   OLD Z-INDEX OF PARTICLE (AT START OF ADVECTION OR 
!                         DIFFUSION)
!    *CLIM*     REAL      MINIMUM VALUE FOR CURRENT (TO PREVENT DIVISION 
!                         BY ZERO)
!    *DCLIM*    REAL      (NUMERICAL) LOWER LIMIT FOR THE COMPONENTS OF THE
!                         VELOCITY GRADIENT (TO PREVENT DIVISION BY ZERO) [1/s]
!    *DTMIN*    REAL      IF IRND>0, EQUAL TO TO THE MINIMUM TIME TO REACH ONE
!                         OF THE BOX EDGES BY ADVECTION                     [s]
!    *DTREST*   REAL      REMAINING TIME STEP EVALUATED WHEN THE ADVECTED 
!                         PARTICLE CROSSED A BOX EDGE                       [s]
!    *DTX*      REAL      IF IRND>0, TIME NEEDED BY THE ADVECTED PARTICLE TO
!                         REACH WESTERN OR EASTERN CELL BOUNDARY            [s]
!    *DTY*      REAL      IF IRND>0, TIME NEEDED BY THE ADVECTED PARTICLE TO
!                         REACH SOUTHERN OR NORTHERN CELL BOUNDARY          [s]
!    *DTZ*      REAL      IF IRND>0, TIME NEEDED BY THE ADVECTED PARTICLE TO
!                         REACH BOTTOM OR TOP CELL BOUNDARY                 [s]
!    *DUDX*     REAL      X-GRADIENT OF U2                                [1/s]
!    *DVDY*     REAL      Y-GRADIENT OF V2                                [1/s]
!    *DWDZ*     REAL      Z-GRADIENT OF W2                                [1/s]
!    *DX2*      REAL      HORIZONTAL GRID SPACING IN X-DIRECTION AT CURRENT
!                         CELL POSITION                                     [m]
!    *DY2*      REAL      HORIZONTAL GRID SPACING IN Y-DIRECTION AT CURRENT
!                         CELL POSITION                                     [m]
!    *DZ2*      REAL      VERTICAL GRID SPACING AT CURRENT CELL POSITION    [m]
!    *EPSMIN*   REAL      NUMERICAL MINIMUM VALUE FOR "X" TO PREVENT INACCURATE
!                         EVALUATION OF LOG(X+1)
!    *EPSMAX*   REAL      NUMERICAL MAXIMUM VALUE FOR "X" TO PREVENT INACCURATE
!                         EVALUATION OF LOG(X+1)
!    *REST*     REAL      TEMPORARY WORK SPACE
!    *WBOT*     REAL      VERTICAL VELOCITY AT BOTTOM CELL FACE           [m/s]
!    *WTOP*     REAL      VERTICAL VELOCITY AT TOP CELL FACE              [m/s]
!    *XADV*     REAL      ADVECTED DISTANCE IN X-DIRECTION DURING THE TIME
!                         STEPS "DTREST" OR "DTMIN"                         [m]
!    *XDIFF*    REAL      DIFFUSED DISTANCE IN X-DIRECTION DURING DEL3      [m]
!    *XEAST*    REAL      DISTANCE TO THE EASTERN CELL BOUNDARY [m]
!    *XEDGE*    REAL      EQUAL TO XEAST IF IEAST=1 OR TO -XWEST IF IWEST=1 [m]
!    *XEND*     REAL      X-COORDINATE WITH RESPECT TO THE OLD CELL CENTRE,
!                         REACHED BY THE ADVECTED PARTICLE AT THE END OF THE 
!                         TIME STEP "DTREST"                                [m]
!    *XOLD*     REAL      OLD X-COORDINATE OF PARTICLE (AT START OF ADVECTION
!                         OR DIFFUSION)                                     [m]
!    *XRAN*     REAL      RANDOM NUMBER BETWEEN 0 AND 1
!    *XTMP*     REAL      TEMPORARY WORK SPACE
!    *XTOT*     REAL      TEMPORARY WORK SPACE
!    *XU2*      REAL      VALUE OF U2 AT OLD PARTICLE POSITION            [m/s]
!    *XWEST*    REAL      DISTANCE TO WESTERN CELL BOUNDARY [m]
!    *YADV*     REAL      ADVECTED DISTANCE IN Y-DIRECTION DURING THE TIME
!                         STEPS "DTREST" OR "DTMIN"                         [m]
!    *YDIFF*    REAL      DIFFUSED DISTANCE IN Y-DIRECTION DURING DEL3      [m]
!    *YEDGE*    REAL      EQUAL TO YNORT IF INORT=1 OR TO -YSOUT IF ISOUT=1 [m]
!    *YEND*     REAL      Y-COORDINATE WITH RESPECT TO THE OLD CELL CENTRE,
!                         REACHED BY THE ADVECTED PARTICLE AT THE END OF THE 
!                         TIME STEP "DTREST"                                [m]
!    *YNORT*    REAL      DISTANCE TO THE NORTHERN CELL BOUNDARY            [m]
!    *YOLD*     REAL      OLD Y-COORDINATE OF PARTICLE (AT START OF ADVECTION
!                         OR DIFFUSION)                                     [m]
!    *YSOUT*    REAL      DISTANCE TO THE SOUTHERN CELL BOUNDARY            [m]
!    *YTMP*     REAL      TEMPORARY WORK SPACE
!    *YTOT*     REAL      TEMPORARY WORK SPACE
!    *YV2*      REAL      VALUE OF V2 AT OLD PARTICLE POSITION            [m/s]
!    *ZADV*     REAL      ADVECTED DISTANCE IN Z-DIRECTION DURING THE TIME
!                         STEPS "DTREST" OR "DTMIN"                         [-]
!    *ZBOT*     REAL      DISTANCE TO THE BOTTOM CELL BOUNDARY              [-]
!    *ZDIFF*    REAL      DIFFUSED DISTANCE IN Z-DIRECTION DURING DEL3      [-]
!    *ZEDGE*    REAL      EQUAL TO ZTOP IF ITOP=1 OR TO -ZBOT IF IBOT=1     [-]
!    *ZEND*     REAL      Z-COORDINATE WITH RESPECT TO THE OLD CELL CENTRE,
!                         REACHED BY THE ADVECTED PARTICLE AT THE END OF THE 
!                         TIME STEP "DTREST"                                [-]
!    *ZOLD*     REAL      OLD SIGMA-COORDINATE OF PARTICLE (AT START OF 
!                         ADVECTION OR DIFFUSION)                           [-]
!    *ZTMP*     REAL      TEMPORARY WORK SPACE
!    *ZTOP*     REAL      DISTANCE TO THE TOP CELL BOUNDARY                 [-]
!    *ZW2*      REAL      VALUE OF W2 AT OLD PARTICLE POSITION DIVIDED
!                         BY TOTAL DEPTH                                  [1/s]
!
!-----------------------------------------------------------------------
!

  !    NERRS = 0

DO I=1,NC
DO J=1,NR
   IF (NWD(J,I)==1) THEN
     DO K=1,NZ
        UVTURB = SQRT(6.0*HEDDYDC(K,J,I)/DEL3)
        VEDC = 0.5*(VEDDYD(K,J,I)+VEDDYD(K+1,J,I))
        UTURB(K,J,I) = UVTURB
        VTURB(K,J,I) = UVTURB
        WTURB(K,J,I) = SQRT(6.0*VEDC/DEL3)
     END DO
   ENDIF
 END DO
 END DO

!
! Old positions
!
DO n = 1,ntrac
  xtraold(n) = REAL(IT(n)-0.5)*gx2(1,1)+XT(n)
  ytraold(n) = REAL(JT(n)-0.5)*gy2(1)+YT(n)
END DO

!     1. START LOOP
!     -------------
!
! loop 1
DO M= 1,ntrac !   1,ntrac

!
!     2. ADVECTIVE TRANSPORT
!     ----------------------
!
!write(6,*)"**********************", trac_start(m)

  N = M
!***********************
IF(trac_start(n))THEN
!**********************

  DTREST = DEL3
  IZAL = 0

!
!     2.1 INITIALISE INDICES AND PARAMETERS
!     -------------------------------------
! loop 2
  DO idum = 1,5
! if 1
    IF((DTREST>0).and.(IT(N)>0))THEN
      IZAL = IZAL + 1

!        ---old grid indices
       IOLD = IT(N)
       JOLD = JT(N)
       KOLD = KT(N)

!        ---current grid indices
       II = IT(N)
       JJ = JT(N)
       KK = KT(N)

!        ---counters
       IWEST = 0
       IEAST = 0
       ISOUT = 0
       INORT = 0
       IBOT = 0
       ITOP = 0
       IRND = 0

!        ---old box coordinates
        XOLD = XT(N)
        YOLD = YT(N)
        ZOLD = ZT(N)

!WRITE(6,*)"**1***",II,JJ,KK,XOLD,YOLD,ZOLD

!        ---box dimensions (divided by 2)
         DX2 = 0.5*GX2(JJ,II)
         DY2 = 0.5*GY2(JJ)
         DZ2 = 0.5*GZSC(KK,JJ,II)

!        ---distance to cell boundaries
         XWEST = DX2 + XOLD
         XEAST = DX2 - XOLD
         YSOUT = DY2 + YOLD
         YNORT = DY2 - YOLD
         ZBOT = DZ2 + ZOLD
         ZTOP = DZ2 - ZOLD

!
!     2.2 ADVECTED DISTANCE WITHIN DTREST
!     -----------------------------------
!
         XRAN = RANDOM(0)
         UDIFF = UTURB(KK,JJ,II)*(XRAN-0.5)*2.0

!        ---X-direction
         DUDX = (U2(KK,JJ,II+1)-U2(KK,JJ,II))/GX2(JJ,II)
         XU2 = UGEO+U2(KK,JJ,II) + XWEST*DUDX+UDIFF

         IF (ABS(XU2)<CLIM) XU2 = 0.0
         IF ((ABS(DUDX)<DCLIM).OR.(ABS(DUDX*DTREST)<DCLIM)) THEN
            XADV = XU2*DTREST*(1.0+0.5*DUDX*DTREST)
         ELSE
            XADV = XU2*(EXP(DUDX*DTREST)-1.0)/DUDX
         ENDIF
         XEND = XOLD + XADV

!        ---Y-direction
         XRAN = RANDOM(0)
         VDIFF = VTURB(KK,JJ,II)*(XRAN-0.5)*2.0

         DVDY = (V2(KK,JJ+1,II)-V2(KK,JJ,II))/GY2(JJ)
         YV2 = V2(KK,JJ,II) + YSOUT*DVDY + VDIFF
         IF (ABS(YV2)<CLIM) YV2 = 0.0
         IF ((ABS(DVDY)<DCLIM).OR.(ABS(DVDY*DTREST)<DCLIM)) THEN
            YADV = YV2*DTREST*(1.0+0.5*DVDY*DTREST)
         ELSE
            YADV = YV2*(EXP(DVDY*DTREST)-1.0)/DVDY
         ENDIF
         YEND = YOLD + YADV

!        ---Z-direction

         XRAN = RANDOM(0)
         WDIFF = WTURB(KK,JJ,II)/H2ATC(JJ,II)*(XRAN-0.5)*2.0

         ZPOSS = REAL(KK-0.5)*gzsc(KK,JJ,II)+ZOLD
         ZPOSS = 0.0

         WBOT = W2(KK,JJ,II)/H2ATC(JJ,II) 
         WTOP = W2(KK+1,JJ,II)/H2ATC(JJ,II)

         DWDZ = (WTOP-WBOT)/GZSC(KK,JJ,II)
         ZW2 = WBOT + ZBOT*DWDZ + WDIFF
         IF (ABS(ZW2)<CLIM) ZW2 = 0.0
         IF ((ABS(DWDZ)<DCLIM).OR.(ABS(DWDZ*DTREST)<DCLIM)) THEN
            ZADV = ZW2*DTREST*(1.0+0.5*DWDZ*DTREST)
         ELSE
            ZADV = ZW2*(EXP(DWDZ*DTREST)-1.0)/DWDZ
         ENDIF
         ZEND = ZOLD + ZADV

!IF(N==10) WRITE(6,*)"idum",idum
!IF(N==10) WRITE(6,*)"ZOLD,ZADV,ZEND,**2",ZOLD,ZADV,ZEND

!        ---particle cannot cross surface/bottom
         IF (KK==NZ) ZEND = MIN(ZEND,0.9*DZ2)
         IF (KK==1) ZEND = MAX(ZEND,-0.9*DZ2)

!IF(N==10) WRITE(6,*)"ZEND,DZ2,**3",ZEND,DZ2
!
!     2.3 ADVECTION TERMINATED IF NO CELL BOUNDARY IS CROSSED
!     -------------------------------------------------------
!
! if 2

   IF ( (ABS(XEND)<=DX2).AND.(ABS(YEND)<=DY2).AND.(ABS(ZEND)<=DZ2) ) THEN
      XT(N) = XEND
      YT(N) = YEND
      ZT(N) = ZEND
      DTREST = 0.0

!WRITE(6,*)"**2***",II,JJ,KK,XEND,YEND,ZEND

   ELSE

!
!     2.4 LOCATION AND NUMBER OF CROSSED BOUNDARIES
!     ---------------------------------------------
!
!        ---east
         IF (XEND>DX2) THEN
            IEAST = 1
            II = II + 1
            IRND = IRND + 1
!        ---west
         ELSEIF (XEND<-DX2) THEN
            IWEST = 1
            II = II - 1
            IRND = IRND + 1
         ENDIF
!        ---north
         IF (YEND>DY2) THEN
            INORT = 1
            JJ = JJ + 1
            IRND = IRND + 1
!        ---south
         ELSEIF (YEND<-DY2) THEN
            ISOUT = 1
            JJ = JJ - 1
            IRND = IRND + 1
         ENDIF
!        ---top
         IF (ZEND>DZ2) THEN
            ITOP = 1
            KK = KK + 1
            IRND = IRND + 1
!        ---bottom
         ELSEIF (ZEND<=-DZ2) THEN
            IBOT = 1
            KK = KK - 1
            IRND = IRND + 1
         ENDIF

!IF(N==10) WRITE(6,*)"KOLD,KK,**4b",KOLD,KK
!
!     2.5 EVALUATE TIME NEEDED TO REACH X-, Y- OR Z- CELL FACE
!     --------------------------------------------------------
!
!        ---X-boundary
         XEDGE = XEAST*IEAST - XWEST*IWEST
         XTMP = XEDGE*DUDX
         IF (XU2==0.0) THEN
            DTX = 0.0
         ELSEIF (ABS(XTMP)<EPSMIN*ABS(XU2)) THEN
            DTX = XEDGE/XU2*(1.0-0.5*XTMP/XU2)
         ELSE
            DTX = (LOG(XTMP/XU2+1.0))/DUDX
         ENDIF
!        ---Y-boundary
         YEDGE = YNORT*INORT - YSOUT*ISOUT
         YTMP = YEDGE*DVDY
         IF (YV2==0.0) THEN
            DTY = 0.0
         ELSEIF (ABS(YTMP)<EPSMIN*ABS(YV2)) THEN
            DTY = YEDGE/YV2*(1.0-0.5*YTMP/YV2)
         ELSE
            DTY = (LOG(YTMP/YV2+1.0))/DVDY
         ENDIF
!        ---Z-boundary
         ZEDGE = ZTOP*ITOP - ZBOT*IBOT
         ZTMP = ZEDGE*DWDZ
         IF (ZW2==0.0) THEN
            DTZ = 0.0
         ELSEIF (ABS(ZTMP)<EPSMIN*ABS(ZW2)) THEN
            DTZ = ZEDGE/ZW2*(1.0-0.5*ZTMP/ZW2)
         ELSE
            DTZ = (LOG(ZTMP/ZW2+1.0))/DWDZ
         ENDIF

!
!     2.6 DEFINE DTMIN AS THE SMALLEST TIME TO REACH A CELL BOUNDARY
!     --------------------------------------------------------------
!
         IF (IRND>1) THEN
            IF (KK==KOLD) THEN
               DTZ = DTX + DTY
            ELSEIF (JJ==JOLD) THEN
               DTY = DTX + DTZ
            ELSEIF (II==IOLD) THEN
               DTX = DTY + DTZ
            ENDIF
            DTMIN = MIN(DTX,DTY,DTZ)
         ELSE
            DTMIN = (IWEST+IEAST)*DTX+(ISOUT+INORT)*DTY+(IBOT+ITOP)*DTZ
         ENDIF

!
!     2.7 RESET BOX INDICES (IF NECESSARY)
!     ------------------------------------
!

         IF (IRND>1) THEN
!           ---X-boundary
            IF (DTX>DTMIN) THEN
               II = II - IEAST + IWEST
               IEAST = 0
               IWEST = 0
            ENDIF
!           ---Y-boundary
            IF (DTY>DTMIN) THEN
               JJ = JJ - INORT + ISOUT
               INORT = 0
               ISOUT = 0
            ENDIF
!           ---Z-boundary
            IF (DTZ>DTMIN) THEN
               KK = KK - ITOP + IBOT
               ITOP = 0
               IBOT = 0
            ENDIF
         ENDIF
!
!
!     2.10 EVALUATE NEW BOX COORDINATES
!     ---------------------------------
!
!        ---X-coordinate
         IF ((ABS(DUDX)<DCLIM).OR.(ABS(DUDX*DTMIN)<DCLIM)) THEN
            XADV = XU2*DTMIN*(1.0+0.5*DUDX*DTMIN)
         ELSE
            XADV = XU2*(EXP(DUDX*DTMIN)-1.0)/DUDX
         ENDIF
         IF (IWEST==1.OR.IEAST==1) THEN
            XT(N) = 0.5*GX2(JJ,II)*(IWEST-IEAST)
         ELSE
            XT(N) = XOLD + XADV
            XT(N) = MAX(-0.5*GX2(JJ,II),XT(N))
            XT(N) = MIN(0.5*GX2(JJ,II),XT(N))
         ENDIF
         IT(N) = II

!        ---Y-coordinate
         IF ((ABS(DVDY)<DCLIM).OR.(ABS(DVDY*DTMIN)<DCLIM)) THEN
            YADV = YV2*DTMIN*(1.0+0.5*DVDY*DTMIN)
         ELSE
            YADV = YV2*(EXP(DVDY*DTMIN)-1.0)/DVDY
         ENDIF
         IF (ISOUT==1.OR.INORT==1) THEN
            YT(N) = 0.5*GY2(JJ)*(ISOUT-INORT)
         ELSE
            YT(N) = YOLD + YADV
            YT(N) = MAX(-0.5*GY2(JJ),YT(N))
            YT(N) = MIN(0.5*GY2(JJ),YT(N))
         ENDIF
         JT(N) = JJ

!        ---Z-coordinate
         IF ((ABS(DWDZ)<DCLIM).OR.(ABS(DWDZ*DTMIN)<DCLIM)) THEN
            ZADV = ZW2*DTMIN*(1.0+0.5*DWDZ*DTMIN)
         ELSE
            ZADV = ZW2*(EXP(DWDZ*DTMIN)-1.0)/DWDZ
         ENDIF
         IF (IBOT==1.OR.ITOP==1) THEN
            ZT(N) = 0.5*GZSC(KK,JJ,II)*(IBOT-ITOP)
!IF(N==10) WRITE(6,*)"ZT(N),**5a",ZT(N)
         ELSE
            ZT(N) = ZOLD + ZADV
            ZT(N) = MAX(-0.5*GZSC(KK,JJ,II),ZT(N))
            ZT(N) = MIN(0.5*GZSC(KK,JJ,II),ZT(N))
!IF(N==10) WRITE(6,*)"ZT(N),**5a",ZT(N)
         ENDIF

         KT(N) = KK
!

!     2.11 EVALUATE NEW "REST" TIME
!     -----------------------------
!

         DTREST = DTREST - DTMIN
         DTREST = MAX(0.0,DTREST)
! end if 2    
    END IF

 ! end if 1
  END IF

! end do 2
END DO

!WRITE(6,*)"**3***",II,JJ,KK,XT(N),YT(N),ZT(N)


!
!     3. DIFFUSION
!     ------------
!     3.1 REINITIALISE INDICES
!     ------------------------
!        ---old indices
!         IOLD = IT(N)
!         JOLD = JT(N)
!         KOLD = KT(N)
!        ---current indices
!         II = IT(N)
!         JJ = JT(N)
 !        KK = KT(N)
!     3.2 X-DIFFUSION
!     ---------------
!        ---evaluate diffused distance
 !        XRAN = RANDOM(0)
 !        XDIFF = DEL3*UTURB(KOLD,JOLD,IOLD)*(XRAN-0.5)*2.0
 !        XTOT = XT(N) + XDIFF
!        ---new X-coordinate
 !        IF (ABS(XTOT)<=0.5*GX2(JJ,II)) THEN
!           --no boundary crossing
 !           XT(N) = XTOT
 !        ELSEIF (XTOT<=(-0.5*GX2(JJ,II))) THEN
!           --westward diffusion
 !           IF (NPIX(JJ,II)==1) THEN
 !              XT(N) = XTOT + GX2U(JJ,II)
 !              II = II - 1
 !           ELSE
!              --reflection at solid/open boundary
 !              XT(N) = -XDIFF - GX2(JJ,II) - XT(N)
 !           ENDIF
  !       ELSEIF (XTOT>0.5*GX2(JJ,II)) THEN
!           --eastward diffusion
 !          IF (NPIX(JJ,II+1)==1) THEN
 !             XT(N) = XTOT - GX2U(JJ,II+1)
 !              II = II + 1
 !           ELSE
!              --reflection at solid/open boundary
  !             XT(N) = GX2(JJ,II) - XDIFF - XT(N)
  !          ENDIF
  !      ENDIF

!
!     3.3 Y-DIFFUSION
!     ---------------
!
!        ---evaluate diffused distance
!         XRAN = RANDOM(0)
!         YDIFF = DEL3*VTURB(KOLD,JOLD,IOLD)*(XRAN-0.5)*2.0
!         YTOT = YT(N) + YDIFF
!        ---new Y-coordinate
 !        IF (ABS(YTOT)<=0.5*GY2(JJ)) THEN
!           --no boundary crossing
  !          YT(N) = YTOT
   !      ELSEIF (YTOT<=(-0.5*GY2(JJ))) THEN
!           --southward diffusion
 !          IF (NPIY(JJ,II)==1) THEN
 !             YT(N) = YTOT + GY2V(JJ)
  !            JJ = JJ - 1
  !         ELSE
!              --reflection at solid/open boundary
   !            YT(N) = -YDIFF - GY2(JJ) - YT(N)
   !         ENDIF
   !      ELSEIF (YTOT>0.5*GY2(JJ)) THEN
!           --northward diffusion
  !          IF (NPIY(JJ+1,II)==1) THEN
  !             YT(N) = YTOT - GY2V(JJ+1)
  !             JJ = JJ + 1
  !          ELSE
!!              --reflection at solid/open boundary
   !            YT(N) = GY2(JJ) - YDIFF - YT(N)
   !         ENDIF
   !      ENDIF
   !      IT(N) = II
   !      JT(N) = JJ
!
!     3.4 VERTICAL DIFFUSION
!     ----------------------
!
!        ---evaluate diffused distance
!  XRAN = RANDOM(0)
!  ZDIFF = DEL3*(WTURB(KOLD,JOLD,IOLD)/H2ATC(JOLD,IOLD))*(XRAN-0.5)*2.0
!  GZZ = GZSC(KK,JJ,II)  
!  POSNEW = (REAL(KK)-0.5)*GZZ+ZT(N)+ZDIFF
!  POSNEW = MAX(0.0,POSNEW) ! stops at sea surface
!  POSNEW = MIN(1.0,POSNEW) ! stops at seafloor
!  IF(POSNEW==1.0)THEN
!    KK = NZ
!    ZT(N) = 0.5*GZZ
!  ELSE
 !   KK = INT(POSNEW*NZ)+1
 !   ZT(N) = POSNEW-REAL(KK)/REAL(NZ)+0.5*GZZ
!        ---vertical grid index
 ! END IF
 ! KT(N) = KK

!***********************
END IF
!**********************

  IF(IT(N)<=2.or.IT(N)>=nc-1) trac_stop(n) = .true.
  IF(JT(n)<=1.or.JT(n)>=nr) trac_stop(n) = .true.
  IF(KT(n)<1.or.KT(n)>nz) trac_stop(n) = .true.

END DO

! speed of floats

DO m = 1,ntrac
  xpos(m) = REAL(IT(m)-0.5)*gx2(1,1)+XT(m)
  ypos(m) = REAL(JT(m)-0.5)*gy2(1)+YT(m)
  dxx = xpos(m) - xtraold(m)
  dyy = ypos(m) - ytraold(m) 
  dist = sqrt(dxx*dxx+dyy*dyy)
  posU(m) = dist/del3
END DO

RETURN

END SUBROUTINE MC3D


END MODULE cohrun



















