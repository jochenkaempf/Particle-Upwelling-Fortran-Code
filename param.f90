MODULE param

INTEGER(4), PARAMETER :: nc = 50
INTEGER(4), PARAMETER :: nr = 40
INTEGER(4), PARAMETER :: nz = 60
INTEGER(4), PARAMETER :: nobu = 78
INTEGER(4), PARAMETER :: nobv = 50
INTEGER(4), PARAMETER :: nvprof = 128
!****************************************
! particle module
INTEGER, PARAMETER :: ntrac = 4000
!****************************************

INTEGER(4), PARAMETER :: ncon = 0
INTEGER(4), PARAMETER :: nconto = 1
INTEGER(4), PARAMETER :: nconc = 0
INTEGER(4), PARAMETER :: maxnop = 0
INTEGER(4), PARAMETER :: noutmax = 5
INTEGER(4), PARAMETER :: nanalmax = 1
INTEGER(4), PARAMETER :: navrmax = 5
INTEGER(4), PARAMETER :: maxout = 1
INTEGER(4), PARAMETER :: maxavr = 1

REAL :: iVOL, iVOLsum, ttime, etam
REAL :: trigger,ad

!-----------------------------------------------------------------------
!
!* LAGRANGIAN MODULE
!

INTEGER :: QTNIT1(1:4),QTNIT2(1:4),QTNIT3(1:4),QXNIT1(1:4),QXNIT2(1:4),QXNIT3(1:4)

INTEGER :: IT(0:ntrac), JT(0:ntrac), KT(0:ntrac)
INTEGER :: LST(0:ntrac), NUMPIN(0:ntrac)
REAL :: posH(0:ntrac), posU(0:ntrac), posC(0:ntrac), xtraold(0:ntrac),ytraold(0:ntrac)
INTEGER :: NUMP
REAL :: SPMCON(NZ,NR,NC), ST(0:MAXNOP), STIN(NZ,NR,NC)
REAL :: STOUT(NZ,NR,NC)
REAL :: UTURB(NZ,NR,NC),VTURB(NZ,NR,NC),WTURB(NZ,NR,NC)
REAL :: XT(0:ntrac), YT(0:ntrac), ZT(0:ntrac)
REAL :: xpos(ntrac),ypos(ntrac),zpos(ntrac),depos(ntrac)
REAL :: XLEN,YLEN,ZLEN
INTEGER :: ITini(0:ntrac), JTini(0:ntrac), KTini(0:ntrac)
REAL :: XTini(0:ntrac), YTini(0:ntrac), ZTini(0:ntrac)
REAL :: xposini(ntrac),yposini(ntrac),zposini(ntrac),ttmin
INTEGER :: ibatch, jmerk(nc)
LOGICAL :: trac_start(ntrac),FLOAT_ON,trac_stop(ntrac)
!-----------------------------------------------------------------------
!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *NC*       INTEGER   NUMBER OF GRID CELLS IN THE X-DIRECTION
!    *NR*       INTEGER   NUMBER OF GRID CELLS IN THE Y-DIRECTION
!    *NZ*       INTEGER   NUMBER OF GRID CELLS IN THE VERTICAL DIRECTION
!    *NOBU*     INTEGER   NUMBER OF OPEN BOUNDARY POINTS AT U-FACES
!    *NOBV*     INTEGER   NUMBER OF OPEN BOUNDARY POINTS AT V-FACES
!    *NVPROF*   INTEGER   NUMBER OF VERTICAL PROFILES AT OPEN BOUNDARIES
!    *NCON*     INTEGER   NUMBER OF TIDAL CONSTITUENTS
!    *NCONTO*   INTEGER   NUMBER OF FREQUENCIES FOR HARMONIC ANALYSIS
!                         (BETWEEN 1 AND 9)
!    *NCONC*    INTEGER   NUMBER OF CONTAMINANTS DISTRIBUTIONS
!                         FOR EULERIAN CONTAMINANT MODULE
!    *MAXNOP*   INTEGER   MAXIMUM ALLOWED NUMBER OF PARTICLES FOR LAGRANGIAN
!                         TRANSPORT MODULE
!    *NOUTMAX*  INTEGER   MAXIMUM NUMBER OF 2D OR 3D FIELDS FOR TIME SERIES
!                         OUTPUT (MIN. VALUE OF 1)
!    *NANALMAX* INTEGER   MAXIMUM NUMBER OF 2D OR 3D FIELDS FOR HARMONIC
!                         ANALYSIS (MIN. VALUE OF 1)
!    *NAVRMAX*  INTEGER   MAXIMUM NUMBER OF 2D OR 3D FIELDS FOR TIME AVERAGED
!                         OUTPUT (MIN. VALUE OF 1)
!    *MAXOUT*   INTEGER   MAXIMUM NUMBER OF TIME SERIES OUTPUT FILES
!                         (BETWEEN 1 AND 9)
!    *MAXAVR*   INTEGER   MAXIMUM NUMBER OF TIME AVERAGED OUTPUT FILES
!                         (BETWEEN 1 AND 9)

!------------------------------------------------------------------------

!*    COMMON *BOUNDS* - OPEN BOUNDARY DATA

INTEGER :: itypobu(0:nobu), itypobv(0:nobv)
INTEGER :: ivpobu(0:nobu), ivpobv(0:nobv)
INTEGER :: lstobu(0:nobu), lstobv(0:nobv)
INTEGER :: ifo(208),jfo(208)

REAL :: ampobu(0:nobu,0:ncon), ampobv(0:nobv,0:ncon)
REAL :: convp(0:nz,0:nvprof,0:nconc)
REAL :: phaobu(0:nobu,0:ncon), phaobv(0:nobv,0:ncon)
REAL :: p2bvp(0:nz,0:nvprof), p2cvp(0:nz,0:nvprof)
REAL :: p2mvp(0:nz,0:nvprof), p2nvp(0:nz,0:nvprof)
REAL :: p2nhsvp(0:nz,0:nvprof),p2nosvp(0:nz,0:nvprof)
REAL :: p2ovp(0:nz,0:nvprof), p2znvp(0:nz,0:nvprof)
REAL :: qstobu(0:nobu), qstobv(0:nobv)
REAL :: rmassu(0:nobu,nz), rmassv(0:nobv,nz)
REAL :: r1obu(0:nobu), r1obv(0:nobv)
REAL :: r2obu(0:nobu), r2obv(0:nobv)
REAL :: sedc1vp(0:nz,0:nvprof), svp(0:nz,0:nvprof)
REAL :: tvp(0:nz,0:nvprof), turvp(0:nz+1,0:nvprof)
REAL :: uvp(0:nz,0:nvprof),cvp(0:nz,0:nvprof)

INTEGER :: ICOAST(1000),JCOAST(1000), NCOAST
LOGICAL :: YES

!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ITYPOBU*  INTEGER   TYPE OF U BOUNDARY COND (0,1,2,3,4)
!    *ITYPOBV*  INTEGER   TYPE OF V BOUNDARY COND (0,1,2,3,4)
!    *IVPOBU*   INTEGER   PROFILE TYPE AT U OPEN BOUNDARY
!    *IVPOBV*   INTEGER   PROFILE TYPE AT V OPEN BOUNDARY
!    *LSTOBU*   INTEGER   PARTICLE LABEL AT U-OPEN BOUNDARIES
!                         (TO TRACE ITS ORIGIN)
!    *LSTOBV*   INTEGER   PARTICLE LABEL AT V-OPEN BOUNDARIES
!                         (TO TRACE ITS ORIGIN)
!    *AMPOBU*   REAL      AMPLITUDE OF INCOMING VARIABLE AT U-NODES       [m]
!    *AMPOBV*   REAL      AMPLITUDE OF INCOMING VARIABLE AT V-NODES       [m]
!    *CONVP*    REAL      VERTICAL PROFILES FOR CONTAMINANTS
!    *PHAOBU*   REAL      PHASE OF INCOMING VARIABLE AT U-NODES         [rad]
!    *PHAOBV*   REAL      PHASE OF INCOMING VARIABLE AT V-NODES         [rad]
!    *P2BVP*    REAL      VERTICAL PROFILES FOR MICROPLANKTON CARBON
!                                                                 [mmol C/m3]
!    *P2CVP*    REAL      VERTICAL PROFILES FOR DETRITAL CARBON   [mmol C/m3]
!    *P2MVP*    REAL      VERTICAL PROFILES FOR DETRITAL NITROGEN [mmol N/m3]
!    *P2NVP*    REAL      VERTICAL PROFILES FOR MICROPLANKTON NITROGEN
!                                                                 [mmol N/m3]
!    *P2NHSVP*  REAL      VERTICAL PROFILES FOR AMMONIUM          [mmol N/m3]
!    *P2NOSVP*  REAL      VERTICAL PROFILES FOR NITRATE           [mmol N/m3]
!    *P2OVP*    REAL      VERTICAL PROFILES FOR OXYGEN            [mmol O/m3]
!    *P2ZNVP*   REAL      VERTICAL PROFILES FOR ZOOPLANKTON NITROGEN
!                                                                 [mmol N/m3]
!    *QSTOBU*   REAL      SPM INPUT AT U-OPEN BOUNDARIES FOR PARTICLE MODULE
!                         [ton/s] AT RIVER, [ton/m3] AT OPEN SEA BOUNDARIES
!    *QSTOBV*   REAL      SPM INPUT AT V-OPEN BOUNDARIES FOR PARTICLE MODULE
!                         [ton/s] AT RIVER, [ton/m3] AT OPEN SEA BOUNDARIES
!    *RMASSU*   REAL      REMAINING MASS TO BE DISTRIBUTED AT NEXT TIME STEP
!                         AT U-OPEN BOUNDARIES                          [ton]
!    *RMASSV*   REAL      REMAINING MASS TO BE DISTRIBUTED AT NEXT TIME STEP
!                         AT V-OPEN BOUNDARIES                          [ton]
!    *R1OBU*    REAL      INCOMING CHARACTERISTIC AT U OPEN BOUNDARIES [m2/s]
!    *R1OBV*    REAL      INCOMING CHARACTERISTIC AT V OPEN BOUNDARIES [m2/s]
!    *R2OBU*    REAL      OUTGOING CHARACTERISTIC AT U OPEN BOUNDARIES [m2/s]
!    *R2OBV*    REAL      OUTGOING CHARACTERISTIC AT V OPEN BOUNDARIES [m2/s]
!    *SEDC1VP*  REAL      VERTICAL PROFILES FOR SEDIMENT CONC.         [g/m3]
!    *SVP*      REAL      VERTICAL PROFILES FOR SALINITY                [PSU]
!    *TVP*      REAL      VERTICAL PROFILES FOR TEMPERATURE           [deg C]
!    *TURVP*    REAL      VERTICAL PROFILES FOR TURBULENCE
!    *UVP*      REAL      VERTICAL PROFILES FOR CURRENT                 [m/s]

!  THE 'QVP' ARRAYS HAVE THE FOLLWING MEANING :
!     QVP(K,*) >  QVP(0,*) => QVP(K,*) IS AN IMPOSED EXTERNAL VALUE OF Q
!     QVP(K,*) <= QVP(0,*) => ZERO NORMAL GRADIENT CONDITION

!*    COMMON *CONCN* - CONCENTRATIONS OF VARIOUS SUBSTANCES

REAL :: buoy(nz,nr,nc), concn1(nz,nr,nc,0:nconc), ro(nz,nr,nc)
REAL :: s(nz,nr,nc), sbet(nz,nr,nc), t(nz,nr,nc), tbet(nz,nr,nc),tr(nz,nr,nc) 
REAL :: sref, tref, roref(nz,nr,nc),ssref(nz,nr,nc),ttref(nz,nr,nc),uuref(nz,nr,nc)
REAL :: cc1ref(nz,nr,nc),cc2ref(nz,nr,nc),cc3ref(nz,nr,nc),cc4ref(nz,nr,nc)
REAL :: rpress(nz+1,nr,nc)
REAL :: shape(nr)

! eulerian tracer module: 3 separate concentration fields 
REAL :: conc1(nz,nr,nc),conc2(nz,nr,nc),conc3(nz,nr,nc),conc4(nz,nr,nc)

!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *BUOY*     REAL      BUOYANCY [m/s2]
!    *CONCN*    REAL      CONCENTRATIONS OF CONTAMINANTS   [g/m3]
!    *RO*       REAL      DENSITY                         [kg/m3]
!    *S*        REAL      SALINITY                          [PSU]
!    *SBET*     REAL      SALINITY EXPANSION COEFFICIENT  [1/PSU]
!    *T*        REAL      TEMPERATURE                     [deg C]
!    *TBET*     REAL      THERMAL EXPANSION COEFFICIENT [1/deg C]

!------------------------------------------------------------------

!*    COMMON *CONSTS* - UNIVERSAL CONSTANTS USED IN PROGRAM

REAL, PARAMETER :: g = 9.81
REAL, PARAMETER :: pi = 3.1415927410
REAL, PARAMETER :: conv = pi/180.0

!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CONV*     REAL      CONVERSION FACTOR DEGREES TO RADIANS (=PI/180)
!    *G*        REAL      GRAVITATIONAL CONSTANT    [m2/s]
!    *PI*       REAL      NUMBER PI
!    *REARTH*   REAL      MEAN RADIUS OF THE EARTH     [m]

!-------------------------------------------------------------

!*    COMMON *COUNT* - COUNTERS FOR MULTIPLE TIME-STEPS

INTEGER :: icbbc, iccbc, ichbc, icmet, icpbc, icwav, ic2bc, ic3d


!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ICBBC*    INTEGER   COUNTER FOR NEW INPUT OF OPEN BOUNDARY DATA
!                         (BIOLOGY, SEDIMENTS)
!    *ICCBC*    INTEGER   COUNTER FOR NEW INPUT OF OPEN BOUNDARY DATA
!                         (CONTAMINANTS)
!    *ICHBC*    INTEGER   COUNTER FOR NEW INPUT OF OPEN BOUNDARY DATA
!                         (3-D CURRENTS, SALINITY, TEMPERATURE)
!    *ICMET*    INTEGER   COUNTER FOR NEW METEOROLOGICAL INPUT
!    *ICPBC*    INTEGER   COUNTER FOR NEW INPUT OF OPEN BOUNDARY DATA
!                         (LAGRANGIAN SPM MODULE)
!    *ICWAV*    INTEGER   COUNTER FOR NEW INPUT OF WAVE DATA
!    *IC2BC*    INTEGER   COUNTER FOR NEW INPUT OF OPEN BOUNDARY DATA
!                         (2-D MODE)
!    *IC3D*     INTEGER   RATIO OF 3-D TO 2-D TIME STEP

!----------------------------------------------------------

!*    COMMON *CRRNTS* - 2-D (DEPTH-INTEGRATED) AND 3-D CURRENTS


REAL :: uadhdev(nr,nc+1), uah2d(nr,nc+1)
REAL :: udh2d(nr,nc+1), udp(nr,nc+1), udqden(nr,nc+1)
REAL :: ud2(nr,nc+1), ud2f(nr,nc+1), uqden(nz,nr,nc+1), uqden0(nz,nr,nc+1)
REAL :: u1(nz,nr,nc+1), u2(nz,nr,nc+1), u2f(nz,nr,nc+1),vort(nz,nr,nc)
REAL :: vadhdev(nr+1,nc), vah2d(nr+1,nc)
REAL :: vdh2d(nr+1,nc), vdp(nr+1,nc), vdqden(nr+1,nc)
REAL :: vd2(nr+1,nc), vd2f(nr+1,nc), vqden(nz,nr+1,nc),vqden0(nz,nr+1,nc)
REAL :: v1(nz,nr+1,nc), v2(nz,nr+1,nc), v2f(nz,nr+1,nc)
REAL :: ugeo
REAL :: w2(nz+1,nr,nc), w2phys(nz,nr,nc)

REAL :: ufinal(nz,nr,nc),vfinal(nz,nr,nc),wfinal(nz+1,nr,nc),fcount
REAL :: tfinal(nz,nr,nc),sfinal(nz,nr,nc),rfinal(nz,nr,nc),pfinal(nz,nr,nc)
REAL :: zetafinal(nr,nc),zetaf(nr,nc)
REAL ::  hedfinal(nz,nr,nc),vedfinal(nz,nr,nc)

!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *UADHDEV*  REAL      DEPTH-INTEGRATED HOR. ADV AND DIFFUSION OF
!                         U2 DEVIATION                               [m2/s2]
!    *UAH2D*    REAL      2-D HORIZONTAL ADVECTION OF UD2            [m2/s2]
!    *UDH2D*    REAL      2-D HORIZONTAL DIFFUSION OF UD2            [m2/s2]
!    *UDP*      REAL      DEPTH-INTEGRATED U-CURRENT AT PREDICTOR STEP)
!                                                                     [m2/s]
!    *UDQDEN*   REAL      DEPTH-INTEGRATED BAROCLINIC PRESSURE GRADIENT
!                         IN X-DIRECTION                             [m2/s2]
!    *UD2*      REAL      DEPTH-INTEGRATED U-CURRENT                  [m2/s]
!    *UD2F*     REAL      DEPTH-INTEGRATED U-CURRENT AVERAGED OVER 3-D TIME
!                         STEP                                        [m2/s]
!    *UQDEN*    REAL      BAROCLINIC PRESSURE GRADIENT IN X-DIRECTION [m/s2]
!    *U1*       REAL      U-CURRENT AT OLD TIME STEP                   [m/s]
!    *U2*       REAL      U-CURRENT AT NEW TIME STEP                   [m/s]
!    *U2F*      REAL      FILTERED ADVECTIVE U-VELOCITY                [m/s]
!    *VADHDEV*  REAL      DEPTH-INTEGRATED HOR. ADV AND DIFFUSION OF
!                         V2 DEVIATION                               [m2/s2]
!    *VAH2D*    REAL      2-D HORIZONTAL ADVECTION OF VD2            [m2/s2]
!    *VDH2D*    REAL      2-D HORIZONTAL DIFFUSION OF VD2            [m2/s2]
!    *VDP*      REAL      DEPTH-INTEGRATED V-CURRENT AT PREDICTOR STEP
!                                                                     [m2/s]
!    *VDQDEN*   REAL      DEPTH-INTEGRATED BAROCLINIC PRESSURE GRADIENT
!                         IN Y-DIRECTION                             [m2/s2]
!    *VD2*      REAL      DEPTH-INTEGRATED V-CURRENT                  [m2/s]
!    *VD2F*     REAL      DEPTH-INTEGRATED V-CURRENT AVERAGED OVER 3-D TIME
!                         STEP                                        [m2/s]
!    *VQDEN*    REAL      BAROCLINIC PRESSURE GRADIENT IN Y-DIRECTION [m/s2]
!    *V1*       REAL      V-CURRENT AT OLD TIME STEP                   [m/s]
!    *V2*       REAL      V-CURRENT AT NEW TIME STEP                   [m/s]
!    *V2F*      REAL      FILTERED ADVECTIVE U-VELOCITY                [m/s]
!    *W2*       REAL      TRANSFORMED VERTICAL VELOCITY AT W-NODES     [m/s]
!    *W2PHYS*   REAL      "PHYSICAL" VERTICAL VELOCITY AT CELL CENTRE  [m/s]

!------------------------------------------

!*    COMMON *DEPTHS* - TOPOGRAPHIC DATA

REAL :: dep(nr,nc),dep2(nr,nc)
REAL :: depun
!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *DEP*      REAL      2-D ARRAY OF MEAN WATER DEPTHS             [m]
!    *DEPUN*    REAL      UNIFORM WATER DEPTH (1-D APPLICATION ONLY) [m]

!--------------------------------------------------------------

!*    COMMON *ELEV* - SURFACE ELEVATIONS (RELATIVE TO UNDISTURBED M.W.L.)

REAL :: zeta1(nr,nc), zeta2(nr,nc), zAD(nr,nc)

!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ZETA1*    REAL      SURFACE ELEVATION AT LOWER  TIME LEVEL [m]
!    *ZETA2*    REAL      SURFACE ELEVATION AT HIGHER TIME LEVEL [m]

!-----------------------------------------------------------------
!*    COMMON *GRID* - MODEL GRID LAYOUT DATA

LOGICAL :: soutob(0:nobv), westob(0:nobu)
INTEGER :: iobu(0:nobu), iobv(0:nobv), jobu(0:nobu), jobv(0:nobv)
INTEGER :: npix(nr,nc+1), npiy(nr+1,nc), nwd(nr,nc)
REAL :: coriol(nr), coriolv(nr+1), cosphi(nr), cosphiv(nr+1)
REAL :: dlat(nr), dlon(nc), gx0(nc+1), gx2(nr,nc)
REAL :: gy0(nr+1), gy2(nr), gz0(nz+1,nr,nc), gz1(nz,nr,nc)
REAL :: gz2(nz,nr,nc), sphcur(nr), sphcurv(nr+1)
REAL :: dlaref, dloref

!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *SOUTOB*   LOGICAL   ORIENTATION OF OPEN BOUNDARY POINTS AT V-NODES
!                         (.TRUE./.FALSE. FOR SOUTH/NORTH)
!    *WESTOB*   LOGICAL   ORIENTATION OF OPEN BOUNDARY POINTS AT U-NODES
!                         (.TRUE./.FALSE. FOR WEST/EAST)
!    *IOBU*     INTEGER   X-INDICES OF OPEN BOUNDARY POINTS AT U-NODES
!    *IOBV*     INTEGER   X-INDICES OF OPEN BOUNDARY POINTS AT V-NODES
!    *JOBU*     INTEGER   Y-INDICES OF OPEN BOUNDARY POINTS AT U-NODES
!    *JOBV*     INTEGER   Y-INDICES OF OPEN BOUNDARY POINTS AT V-NODES
!    *NPIX*     INTEGER   X FACE POINTERS
!                         = 0 => DRY/CLOSED BOUNDARY
!                         = 1 => WET
!                         = 2 => OPEN SEA BOUNDARY
!                         = 3 => RIVER BOUNDARY
!    *NPIY*     INTEGER   Y FACE POINTERS
!                         = 0 => DRY/CLOSED BOUNDARY
!                         = 1 => WET
!                         = 2 => OPEN SEA BOUNDARY
!                         = 3 => RIVER OPEN BOUNDARY
!    *NWD*      INTEGER   CELL POINTERS
!                         = 0 =>  DRY
!                         = 1 =>  WET
!    *CORIOL*   REAL      CORIOLIS FREQUENCY (LATITUDE DEPENDENT IF IGTRH=1)
!                         AT GRID CENTRES                                 [1/s]
!    *CORIOLV*  REAL      CORIOLIS FREQUENCY (LATITUDE DEPENDENT IF IGTRH=1)
!                         AT V-NODES                                      [1/s]
!    *COSPHI*   REAL      COSINE OF LATITUDE AT GRID CENTRES AND U-NODES
!                         (=1 IF IGTRH=0)
!    *COSPHIV*  REAL      COSINE OF LATITUDE AT V-NODES (=1 IF IGTRH=0)
!    *DLAREF*   REAL      REFERENCE LATITUDE IN DECIMAL DEGREES (POSITIVE N)
!    *DLAT*     REAL      LATITUDE AT GRID CENTRES IN DECIMAL DEGREES
!                         (POSITIVE N)
!                         (=DLAREF IF IGTRH=0)
!    *DLON*     REAL      LONGITUDE AT GRID CENTRES IN DECIMAL DEGREES
!                         (POSTIVE E)
!                         (=DLOREF IF IGTRH=0)
!    *DLOREF*   REAL      REFERENCE LONGITUDE IN DECIMAL DEGREES (POSITIVE E)
!    *GX0*      REAL      X-COORDINATES OF CELL CORNERS         [m or degrees]
!    *GX2*      REAL      HORIZONTAL GRID SPACING IN X-DIRECTION           [m]
!    *GY0*      REAL      Y-COORDINATES OF CELL CORNERS         [m or degrees]
!    *GY2*      REAL      HORIZONTAL GRID SPACING IN Y-DIRECTION           [m]
!    *GZ0*      REAL      COORDINATES OF SIGMA GRID (BETWEEN 0 TO 1)
!    *GZ1*      REAL      VERTICAL GRID SPACINGS AT OLD TIME STEP          [m]
!    *GZ2*      REAL      VERTICAL GRID SPACINGS AT NEW TIME STEP          [m]
!    *SPHCUR*   REAL      SPHERICAL CORRECTION FACTOR FOR CURVATURE AT CENTRES
!                         AND U-NODES (=TAN(LATITUDE)/EARTH RADIUS OR =0
!                         IF IGTRH=0                                     [1/m]
!    *SPHCURV*  REAL      SPHERICAL CORRECTION FACTOR FOR CURVATURE AT V-NODES
!                         (=TAN(LATITUDE)/(EARTH RADIUS OR =0 IF IGTRH=0 [1/m]

!----------------------------------------------------------------
!*    COMMON *MET* - METEOROLOGICAL FORCING DATA

REAL :: cloud2(nr,nc), hum2(nr,nc), p2(nr,nc), rain2(nr,nc)
REAL :: sat2(nr,nc), sst2(nr,nc), windu2(nr,nc), windv2(nr,nc)
!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CLOUD2*   REAL      FRACTIONAL CLOUD COVER (BETWEEN 0 AND 1)
!    *HUM2*     REAL      RELATIVE HUMIDITY (BETWEEN 0 AND 1)
!    *P2*       REAL      ATMOSPHERIC PRESSURE            [N/m2]
!    *RAIN2*    REAL      PRECIPITATION RATE         [kg/(m2*s)]
!    *SAT2*     REAL      SURFACE AIR TEMPERATURE        [deg C]
!    *SST2*     REAL      SEA SURFACE TEMPERATURE        [deg C]
!    *WINDU2*   REAL      X-COMPONENT OF 10m WIND VELOCITY [m/s]
!    *WINDV2*   REAL      Y-COMPONENT OF 10m WIND VELOCITY [m/s]

!--------------------------------------------------------------
!     COMMON *OPTICS* - OPTICAL PARAMETERS AND VARIABLES

REAL :: atcfcord(nz,nr,nc), atcfcor2(nz,nr,nc)
REAL :: atcf1(nr,nc), atcf2(nr,nc), hexp(nr,nc)
REAL :: par(nz,nr,nc), qheat(nz,nr,nc), r1opt(nr,nc), r2opt(nr,nc)
!-----------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ATCFCORD* REAL      DIFFUSE ATTENUATION COEFFICIENT FOR PAR      [1/m]
!    *ATCFCOR2* REAL      DIFFUSE ATTENUATION COEFFICIENT FOR SHORT-WAVE
!                         RADIATION                                    [1/m]
!    *ATCF1*    REAL      ATTENUATION COEFFICIENT FOR INFRARED PART
!                         OF SOLAR RADIATION                           [1/m]
!    *ATCF2*    REAL      ATTENUATION COEFFICIENT FOR NON-INFRARED PART
!                         OF SOLAR RADIATION AND PURE SEA WATER        [1/m]
!    *HEXP*     REAL      DEPTH OF THE SURFACE LAYER WHERE ATTENUATION IS
!                         HYPER-EXPONENTIAL                              [m]
!    *PAR*      REAL      PHOTOSYNTHETICALLY AVAILABLE RADIATION      [W/m2]
!    *QHEAT*    REAL      HEATING SOURCE TERM IN TEMPERATURE EQUATION [W/m3]
!    *R1OPT*    REAL      INFRARED FRACTION OF SOLAR RADIATION
!                         ABSORBED NEAR THE SEA SURFACE
!    *R2OPT*      REAL    EXTRA NEAR SURFACE ATTENUATION FACTOR

!--------------------------------------------------------------
!     COMMON  *PARTUR* - SWITCHES FOR TURBULENCE CLOSURE

INTEGER :: iahdht, ileng, ilim, istpar, itcpar, itform, ntrans

!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *IAHDHT*   INTEGER   SWITCH TO DISABLE/ENABLE ADVECTION AND HORIZONTAL
!                         DIFFUSION OF TURBULENCE
!             = 0 => ADVECTION AND HORIZONTAL DIFFUSION DISABLED
!             = 1 => ADVECTION AND HORIZONTAL DIFFUSION ENABLED

!    *ILENG*    INTEGER   SELECTS MIXING LENGTH FORMULATION (IOPTK=2)
!             = 1 => PARABOLIC LAW
!             = 2 => "MODIFIED" PARABOLIC LAW
!             = 3 => "XING" FORMULATION
!             = 4 => "BLACKADAR" FORMULATION

!    *ILIM*     INTEGER   SELECTS LIMITING CONDITION FOR STABLE STRATIFICATION
!                         (IOPTK=2)
!             = 0 => LIMITING CONDITION DISABLED
!             = 1 => LIMITING CONDITION ENABLED

!    *ISTPAR*   INTEGER   SELECTS FORM OF THE STABILITY FUNCTIONS (IOPTK=2)
!             = 1 => AS FUNCTION OF STABILITY PARAMETER G_h (ITCPAR=1)
!                    OR ALPHA_N (ITCPAR=2)
!             = 2 => AS FUNCTION OF RICHARDSON NUMBER IN ANALOGY WITH
!                    MUNK-ANDERSON RELATIONS

!    *ITCPAR*   INTEGER   SELECTS TYPE OF SECOND TURBULENT VARIABLE (IOPTK=2)
!             = 1 => L OR KL
!             = 2 => EPSILON

!    *ITFORM*   INTEGER   SELECTS TYPE OF SCHEME IF IOPTK=1
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

!    *NTRANS*   INTEGER   SELECTS NUMBER OF TRANSPORT EQUATIONS (IOPTK=2)
!             = 0 => LEVEL 2
!             = 1 => LEVEL 2.5: TRANSPORT EQUATION FOR K
!             = 2 => LEVEL 2.5: TRANSPORT EQUATIONS FOR K AND FOR KL (ITCPAR=1)
!                               OR EPSILON (ITCPAR=2)

!----------------------------------------------------------------
!*    COMMON *PHYCON* - PHYSICAL MODEL PARAMETERS

REAL :: atcf1un, atcf2un, ckar, cp, epssal, fsun, gsun, hexpun
REAL :: hsun, r0ref, r1optun, r2optun, sbetun, tbetun
REAL :: twun

!-----------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ATCF1UN*  REAL      DIFFUSE ATTENUATION COEFFICIENT FOR INFRARED PART
!                         OF SOLAR RADIATION                            [m-1]
!    *ATCF2UN*  REAL      DIFFUSE ATTENUATION COEFFICIENT FOR NON-INFRARED
!                         PART OF SOLAR RADIATION AND PURE SEA WATER    [m-1]
!    *CKAR*     REAL      VON KARMAN CONSTANT
!    *CP*       REAL      SPECIFIC HEAT OF SEAWATER AT CONSTANT PRESSURE
!                                                                  [J/kg/K)]
!    *EPSSAL*   REAL      SLOPE OF THE LINEAR RELATIONSHIP BETWEEN THE DIFFUSE
!                         ATTENUATION COEFFICIENT AND SALINITY      [1/PSU/m]
!    *FSUN*     REAL      UNIFORM X-COMP. OF SURFACE STRESS           [m2/s2]
!    *GSUN*     REAL      UNIFORM Y-COMP. OF SURFACE STRESS           [m2/s2]
!    *HEXPUN*   REAL      UNIFORM DEPTH OF THE SURFACE LAYER WHERE ATTENUATION
!                         IS HYPER-EXPONENTAL                             [m]
!    *HSUN*     REAL      UNIFORM WAVE HEIGHT                             [m]
!    *R0REF*    REAL      REFERENCE DENSITY                           [kg/m3]
!    *R1OPTUN*  REAL      UNIFORM INFRARED FRACTION OF SOLAR RADIATION
!                         ABSORBED NEAR THE SEA SURFACE
!    *R2OPTUN*  REAL      UNIFORM EXTRA NEAR SURFACE ATTENUATION FACTOR
!    *SBETUN*   REAL      UNIFORM SALINITY EXPANSION COEFFICIENT      [1/PSU]
!    *SREF*     REAL      REFERENCE SALINITY                            [PSU]
!    *TBETUN*   REAL      UNIFORM THERMAL EXPANSION COEFFICIENT     [1/deg C]
!    *TREF*     REAL      REFERENCE TEMPERATURE                       [deg C]
!    *TWUN*     REAL      UNIFORM WAVE PERIOD                             [s]

!-----------------------------------------------------------------
!*    COMMON *RUNP* - MODEL RUN PARAMETERS

CHARACTER (LEN=6) :: title
INTEGER :: iadvc, iadvs, iadvwb,ibstr, idrag, ifluff,igrdim
INTEGER :: igtrh, iodif, ioptb, ioptc, ioptd, iopthe,ioptk
INTEGER :: ioptm, ioptp, iopts, ioptsa,ioptw, iopt2, iopt3
INTEGER :: ioutf, ioutp, iouts, ioutt, itdif

!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *IADVC*    INTEGER   SWITCH TO SELECT ADVECTION SCHEME FOR MOMENTUM
!                         = 0 -> NO ADVECTION
!                         = 1 -> UPWIND
!                         = 2 -> LAX-WENDROFF (HORIZ.) OR CENTRAL (VERTICAL)
!                         = 3 -> TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                         = 4 -> TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!    *IADVS*    INTEGER   SWITCH TO SELECT ADVECTION SCHEME FOR SCALARS
!                         = 0 -> NO ADVECTION
!                         = 1 -> UPWIND
!                         = 2 -> LAX-WENDROFF (HORIZ.) OR CENTRAL (VERTICAL)
!                         = 3 -> TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                         = 4 -> TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!    *IADVWB*   INTEGER   SWITCH TO SELECT VERTICAL ADVECTION SCHEME
!                         FOR BIOLOGY IF IADVS = 0
!                         = 0 -> NO VERTICAL ADVECTION
!                         = 1 -> UPWIND
!                         = 2 -> CENTRAL
!                         = 3 -> TVD SCHEME WITH SUPERBEE LIMITING FUNCTION
!                         = 4 -> TVD SCHEME WITH MONOTONIC LIMITING FUNCTION
!    *IBSTR*    INTEGER   SWITCH TO SELECT BOTTOM STRESS SCHEME
!                         = 0 -> ZERO BOTTOM STRESS
!                         = 1 -> QUADRATIC BOTTOM STRESS USING CONSTANT
!                                BOTTOM ROUGHNESS LENGTH
!                         = 2 -> LINEAR BOTTOM STRESS
!    *IDRAG*    INTEGER   SWITCH TO SELECT WIND DEPENDENCE OF SURFACE DRAG
!                         COEFFICIENT
!                         = 0 -> CONSTANT VALUE
!                         = 1 -> LARGE AND POND (1981)
!                         = 2 -> SMITH AND BANKE (1975)
!                         = 3 -> GEERNAERT ET AL. (1986)
!                         = 4 -> CHARNOCK'S RELATION
!    *IFLUFF*   INTEGER   SWITCH TO DISABLE/ENABLE "FLUFF LAYER" IN BIOLOGICAL
!                         MODEL (0/1)
!    *IGRDIM*   INTEGER   SWITCH TO SELECT THE DIMENSION OF THE GRID
!                         = 1 -> 1-D APPLICATION IN THE VERTICAL
!                         = 3 -> 3-D APPLICATION
!    *IGTRH*    INTEGER   SWITCH TO SELECT HORIZONTAL COORDINATE TRANSFORM
!                         =0 -> CARTESIAN
!                         =1 -> SPHERICAL (LONGITUDE, LATITUDE)
!    *IODIF*    INTEGER   SWITCH TO SELECT HORIZONTAL DIFFUSION COEFFICIENTS
!                         = 0 -> HORIZONTAL DIFFUSION DISABLED
!                         = 1 -> UNIFORM VALUES
!                         = 2 -> AS FUNCTION OF STRAIN RATE
!                                (FOLLOWING SMAGORINSKY)
!    *IOPTB*    INTEGER   SWITCH TO SELECT BIOLOGICAL MODEL OPTION (0/1)
!    *IOPTC*    INTEGER   SWITCH TO SELECT CONTAMINANT DISPERSION OPTION
!               = 0 -> CONTAMINANT TRANSPORT DISABLED
!               = 1 -> CONTAMINANT TRANSPORT ENABLED WITHOUT VERT/HOR DIFFUSION
!               = 2 -> CONTAMINANT TRANSPORT ENABLED WITH VERT/HOR DIFFUSION
!    *IOPTD*    INTEGER   SWITCH TO SELECT EQUATION OF STATE
!                         = 0 -> UNIFORM DENSITY
!                         = 1 -> USING LINEAR EQUATION OF STATE
!                         = 2 -> USING INTERNATIONAL EQUATION OF STATE
!    *IOPTHE*   INTEGER   SWITCH TO SELECT TEMPERATURE EQUATION
!                         = 0 -> NO TEMPERATURE CALCULATIONS
!                         = 1 -> WITH SOLAR ENERGY ABSORBED AT THE SEA SURFACE
!                         = 2 -> WITH SOLAR ENERGY ABSORBED WITHIN THE WATER
!                                COLUMN
!    *IOPTK*    INTEGER   SWITCH TO SELECT EDDY VISCOSITY/DIFFUSIVITY FORM
!                         = 0 -> CONSTANT VALUES
!                         = 1 -> ALGEBRAIC EXPRESSIONS SELECTED BY ITFORM
!                         = 2 -> USING TUBULENCE ENERGY MODEL SELECTED BY
!                                NTRANS, ITCPAR, ISTPAR, ILENG, ILIM, IAHDHT
!    *IOPTM*    INTEGER   SWITCH TO SELECT TYPE OF MET INPUT
!                 = 0 -> NO MET SURFACE FORCING
!                 = 1 -> UNIFORM IN SPACE, CONSTANT IN TIME, ZERO HEAT FLUXES
!                 = 2 -> UNIFORM IN SPACE, NON-CONSTANT IN TIME
!                        MET. DATA READ FROM DATA FILE, SURFACE STRESS,
!                        HEAT AND SALINITY FLUXES CALC. BY THE PROGRAM
!                 = 3 -> NON-UNIFORM IN SPACE, NON-CONSTANT IN TIME
!                        MET. DATA READ FROM DATA FILE, SURFACE STRESS,
!                        HEAT AND SALINITY FLUXES CALC. BY THE PROGRAM
!                 = 4 -> UNIFORM IN SPACE, NON-CONSTANT IN TIME
!                        WIND DATA, HEAT AND SALINITY FLUXES READ FROM
!                        DATA FILE, WIND STRESS CALC. BY THE PROGRAM
!                 = 5 -> NON-UNIFORM IN SPACE, NON-CONSTANT IN TIME
!                        WIND DATA, HEAT AND SALINITY FLUXES READ FROM
!                        DATA FILE, WIND STRESS CALC. BY THE PROGRAM
!    *IOPTP*    INTEGER   SWITCH TO SELECT LAGRANGIAN PARTICLE TRANSPORT
!                 = 0 -> PARTICLE TRANSPORT DISABLED
!                 = 1 -> PARTICLE TRANSPORT ENABLED WITHOUT VERT/HOR DIFFUSION
!                 = 2 -> PARTICLE TRANSPORT ENABLED WITH VERT/HOR DIFFUSION
!    *IOPTS*    INTEGER   SWITCH TO SELECT SEDIMENT TRANSPORT (0/1)
!    *IOPTSA*   INTEGER   SWITCH TO SELECT SALINITY EQUATION
!                         = 0 -> NO SALINITY CALCULATIONS
!                         = 1 -> USING PRESCRIBED SALINITY SURFACE FLUX
!                         = 2 -> SURFACE SALINITY FLUX EVALUATED USING
!                                EVAPORATION/PRECIPITATION BALANCE
!    *IOPTW*    INTEGER   SWITCH TO SELECT WAVE INPUT
!                         = 0 -> NO WAVE EFFECTS
!                         = 1 -> CONSTANT WAVE INPUT (WAVE HEIGHT, PERIOD)
!                         = 2 -> WAVE DATA UNIFORM IN SPACE,
!                                NON-CONSTANT IN TIME
!                         = 3 -> WAVE DATA NON-UNIFORM IN SPACE,
!                                NON-CONSTANT IN TIME
!    *IOPT2*    INTEGER   SWITCH TO SELECT 2-D CURRENTS (0,1)
!    *IOPT3*    INTEGER   SWITCH TO SELECT 3-D CURRENTS (0,1)
!    *IOUTF*    INTEGER   SWITCH TO SELECT TIME-AVERAGED OUTPUT (0,1)
!    *IOUTP*    INTEGER   SWITCH TO SELECT PARTICLE OUTPUT (0,1)
!    *IOUTS*    INTEGER   SWITCH TO SELECT TIME SERIES OUTPUT
!    *IOUTT*    INTEGER   SWITCH TO SELECT HARMONIC OUTPUT
!                         = 0 -> HARMONIC ANALYSIS ENABLED
!                         = 1 -> HARMONIC OUTPUT OF USER-DEFINED VARIABLES
!                         = 2 -> HARMONIC OUTPUT OF USER-DEFINED VARAIBLES
!                                AND TIDAL ELLIPSE PARAMETERS
!    *ITDIF*    INTEGER   SWITCH TO SELECT DEPENDENCE OF DRAG AND THERMAL
!                         EXCHANGE COEFFICIENT ON AIR/SEA TEMPERATURE
!                         DIFFERENCE (0,1)
!    *TITLE*    CHARACTER TITLE OF MODEL RUN

!---------------------------------------------------------------
!*    COMMON *STRESS* - SURFACE AND BOTTOM STRESSES,SURFACE FLUXES

REAL :: bstot(nr,nc), cdb(nr,nc), cdb100(nr,nc), cdz0(nr,nc)
REAL :: evapr(nr,nc), fb(nr,nc+1), fbk(nr,nc+1), fs(nr,nc)
REAL :: gb(nr+1,nc), gbk(nr+1,nc), gs(nr,nc), qnsol(nr,nc)
REAL :: qsol(nr,nc), ssalfl(nr,nc), sstot(nr,nc)
REAL :: cdlin, cdz0un

!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *BSTOT*    REAL      BOTTOM STRESS  MAGNITUDE AT CELL CENTRE      [m2/s2]
!    *CDB*      REAL      QUADRATIC BOTTOM DRAG COEFFICIENT AT CENTRE OF
!                         BOTTOM GRID CELL
!    *CDB100*   REAL      QUADRATIC BOTTOM DRAG COEFFICIENT (CELL CENTRE)
!                         AT A REFERENCE HEIGHT OF 1m ABOVE THE BED
!    *CDLIN*    REAL      (UNIFORM) LINEAR BOTTOM FRICTION COEFFICIENT   [m/s]
!    *CDZ0*     REAL      ARRAY OF BOTTOM ROUGHNESS LENGTHS                [m]
!    *CDZ0UN*   REAL      UNIFORM BOTTOM ROUGHNESS LENGTH                  [m]
!    *EVAPR*    REAL      EVAPORATION MINUS PRECIPITATION RATE       [kg/m2/s]
!    *FB*       REAL      (MINUS) X-COMPONENT OF BOTTOM STRESS AT U-NODES
!                                                                      [m2/s2]
!    *FBK*      REAL      X-COMPONENT OF BOTTOM FRICTION COEFFICIENT AT
!                         U-NODES                                        [m/s]
!    *FS*       REAL      X-COMPONENT OF SURFACE STRESS AT CELL CENTRE [m2/s2]
!    *GB*       REAL      (MINUS) Y-COMPONENT OF BOTTOM STRESS AT V-NODES
!                                                                      [m2/s2]
!    *GBK*      REAL      Y-COMPONENT OF BOTTOM FRICTION COEFFICIENT AT
!                         V-NODES                                        [m/s]
!    *GS*       REAL      Y-COMPONENT OF SURFACE STRESS AT CELL CENTRE [m2/s2]
!    *QNSOL*    REAL      UPWARD NON-SOLAR HEAT FLUX AT SEA SURFACE     [W/m2]
!    *QSOL*     REAL      DOWNWARD SOLAR HEAT FLUX AT SEA SURFACE       [W/m2]
!    *SSALFL*   REAL      SURFACE SALINITY FLUX [PSU kg/(m2/s]
!    *SSTOT*    REAL      SURFACE STRESS MAGNITUDE AT CELL CENTRE      [m2/s2]

!-----------------------------------------------------------------
!*    COMMON *TIDE* - TIDAL PARAMETERS

REAL :: phase0(0:ncon), sigma(0:ncon)

!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *PHASE0*   REAL      PHASE CORRECTION FACTOR (FOR ASTRONOMICAL TIDE) [rad]
!    *SIGMA*    REAL      FREQUENCIES OF TIDAL CONSTITUENTS             [rad/s]

!------------------------------------------------------------------
!*    COMMON *TIME* - TIME VARIABLES

INTEGER :: ibdate, ibyear, idate, iedate, ieyear, iyear, nstep, nt
REAL :: delm, delt, delw, del3, dtmax, hour

!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *IBDATE*   INTEGER   BEGIN DATE (MMDDHHMM - MONTH, DAY, HOUR, MINUTE)
!    *IBYEAR*   INTEGER   START YEAR
!    *IDATE*    INTEGER   CURRENT DATE (MMDDHHMM - MONTH, DAY, HOUR, MINUTE)
!    *IEDATE*   INTEGER   END DATE (MMDDHHMM - MONTH, DAY, HOUR, MINUTE)
!    *IEYEAR*   INTEGER   END YEAR
!    *IYEAR*    INTEGER   CURRENT YEAR
!    *NSTEP*    INTEGER   TOTAL NO. OF 2-D TIME STEPS
!    *NT*       INTEGER   (2-D) TIME STEP COUNTER
!    *DELM*     REAL      MET. INPUT TIME STEP (INTEGER MULTIPLE OF DELT) [s]
!    *DELT*     REAL      2-D TIME STEP                                   [s]
!    *DELW*     REAL      WAVE INPUT TIME STEP (INTEGER MULTIPLE OF DELT) [s]
!    *DEL3*     REAL      TIME-STEP FOR 3-D MODE                          [s]
!    *DTMAX*    REAL      MAXIMUM 2-D TIME STEP ALLOWED BY CFL-STABILITY LIMIT
!                                                                         [s]
!    *HOUR*     REAL      CURRENT TIME IN HOURS SINCE START OF SIMULATION

!--------------------------------------------------------------------
!*    COMMON *TURBKE* - TURBULENCE VARIABLES

REAL :: buprod(nz+1,nr,nc), dissw(nz+1,nr,nc)
REAL :: shb(nz+1,nr,nc), shprod(nz+1,nr,nc)
REAL :: smu(nz+1,nr,nc), tkeold(nz+1,nr,nc)
REAL :: tkew(nz+1,nr,nc), zlw(nz+1,nr,nc)

!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *BUPROD*   REAL      BUOYANCY SOURCE OR SINK TERM IN TUTRBULENCE
!                         KINETIC ENERGY EQUATION AT W NODES       [W/kg]
!    *DISSW*    REAL      TURBULENCE DISSIPATION RATE AT W NODES   [W/kg]
!    *SHB*      REAL      STABILITY COEFFICIENT AT W-NODES
!                         (S_h FOR K-L OR =S_b FOR K-EPS)             [-]
!    *SHPROD*   REAL      SHEAR PRODUCTION OF TURBULENCE KINETIC ENERGY
!                         AT W-NODES                               [W/kg]
!    *SMU*      REAL      STABILITY COEFFICIENT AT W-NODES
!                         (S_m FOR K-L OR =S_u FOR K-EPS)             [-]
!    *TKEOLD*   REAL      TURBULENCE KINETIC ENERGY AT W-NODES AND OLD TIME
!                         STEP                                      [J/kg]
!    *TKEW*     REAL      TURBULENCE KINETIC ENERGY AT W-NODES AND NEW TIME
!                         STEP                                      [J/kg]
!    *ZLW*      REAL      TURBULENCE MIXING LENGTH AT W-NODES          [m]

!------------------------------------------------------------------------
!*    COMMON *TURCON* - PARAMETERS FOR TURBULENCE MODELS

!*    COMMON *TURCON1* - ALGEBRAIC FORMULATIONS (IOPTK=1)

REAL :: adcnu, add1, add2, adk1, adk2, adlam, adr1, adr2
REAL :: ama, amb, amn1, amn2, amved0
REAL :: ppa, ppdifbg, ppn, ppved0, ppvisbg, rmaxlat, rmaxnut

!-----------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ADCNU*    REAL      C_NU COEFFICIENT FOR "DELTA*U"-LAW (ITFORM=5)
!    *ADD1*     REAL      FRACTIONAL DEPTH OF BOTTOM LAYER WHERE
!                         EDDY VISCOSITY IS REDUCED (ITFORM=3,4,5)
!    *ADD2*     REAL      FRACTIONAL DEPTH OF SURFACE LAYER WHERE
!                         EDDY VISCOSITY IS REDUCED (ITFORM=3,4,5)
!    *ADK1*     REAL      K1-COEFFICIENT FOR "U*H"-LAW (ITFORM=3)
!    *ADK2*     REAL      K2-COEFFICIENT FOR "U*U"-LAW (ITFORM=4)
!    *ADLAM*    REAL      SURFACE ROUGHNESS LENGTH (ITFORM=3,4,5)         [m]
!    *ADR1*     REAL      RATIO OF BOTTOM TO INTERIOR EDDY VISCOSITY
!                         (ITFORM=3,4,5)
!    *ADR2*     REAL      RATIO OF SURFACE TO INTERIOR EDDY VISCOSITY
!                         (ITFORM=3,4,5)
!    *AMA*      REAL      PARAMETER ALPHA_m IN MUNK-ANDERSON RELATION
!                         (ITFORM=2)
!    *AMB*      REAL      PARAMETER BETA_m IN MUNK-ANDERSON RELATION (ITFORM=2)
!    *AMN1*     REAL      N1-EXPONENT IN MUNK-ANDERSON RELATION (ITFORM=2)
!    *AMN2*     REAL      N2-EXPONENT IN MUNK-ANDERSON RELATION (ITFORM=2)
!    *AMVED0*   REAL      NEUTRAL VALUE OF EDDY COEFFICIENT IN MUNK-ANDERSON
!                         RELATION (ITFORM=2)                          [m2/s]
!    *PPA*      REAL      PARAMETER ALPHA_p IN PACONOWSKI-PHILANDER
!                         FORMULATION (ITFORM=1)
!    *PPDIFBG*  REAL      BACKGROUND DIFFUSIVITY IN PACANOWSKI-PHILANDER
!                         FORMULATION (ITFORM=1)                       [m2/s]
!    *PPN*      REAL      Np-EXPONENT IN PACANOWSKI-PHILANDER FORMULATION
!                         (ITFORM=1)
!    *PPVED0*   REAL      NEUTRAL VALUE OF EDDY VISCOSITY IN
!                         PACANOWSKI-PHILANDER FORMLULATION (ITFORM=1) [m2/s]
!    *PPVISBG*  REAL      BACKGROUND VISCOSITY IN PACANOWSKI-PHILANDER
!                         FORMULATION (ITFORM=1)                       [m2/s]
!    *RMAXLAT*  REAL      MAX VALUE OF THE RATIO OF THE EDDY DIFFUSIVITY TO ITS
!                         NEUTRAL VALUE IN THE MUNK-ANDERSON RELATIONS
!    *RMAXNUT*  REAL      MAX VALUE OF THE RATIO OF THE EDDY VISCOSITY TO ITS
!                         NEUTRAL VALUE IN THE PACONOWSKI-PHILANDER AND
!                         MUNK-ANDERSON RELATIONS

!-----------------------------------------------------------------------

!*    COMMON *TURCON2* - TURBULENCE ENERGY MODELS (IOPTK=2)

REAL :: clev2(4), clev25(6)
REAL :: ablac, bxing, cmu, c1e, c2e, c3e1, c3e2, eps0, e1, e2
REAL :: gamax, gamin, sige, sigk, tkemin, z0bot, z0sur

COMMON/turcon2/ clev2, clev25, ablac, bxing, cmu, c1e, c2e, c3e1,  &
    c3e2, eps0, e1, e2, gamax, gamin, sige, sigk, tkemin, z0bot, z0sur

!-----------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ABLAC*    REAL      CONSTANT IN BLACKADAR MIXING LENGTH FORMULATION
!    *BXING*    REAL      CONSTANT IN "XING" MIXING LENGTH FORMULATION
!    *CLEV2*    REAL      COEFFICIENTS FOR CALCULATING THE SQUARED
!                         TURBULENCE FREQUENCY (LEVEL 2)
!    *CLEV25*   REAL      COEFFICIENTS FOR CALCULATING THE STABILITY
!                         FUNCTIONS SMU,SHB
!    *CMU*      REAL      VALUE OF SMU FOR NON-STRATIFIED FLOWS (K-EPS)
!    *C1E*      REAL      EPS-EQUATION (SHEAR PRODUCTION TERM)
!    *C2E*      REAL      EPS-EQUATION (DISSIPATION TERM)
!    *C3E1*     REAL      EPS-EQUATION (BUOYANCY TERM, STABLE CASE)
!    *C3E2*     REAL      EPS-EQUATION (BUOYANCY TERM, UNSTABLE CASE)
!    *EPS0*     REAL      =ZL*DISS/TKE**1.5 (CONVERSION FACTOR)
!    *E1*       REAL      KL-EQUATION (SHEAR/BUOYANCY TERMS)
!    *E2*       REAL      KL-EQUATION (WALL PROXIMITY FUNCTION)
!    *GAMAX*    REAL      MAXIMUM VALUE OF G_h OR ALPHA_N (STABLE STRAT)
!    *GAMIN*    REAL      MINIMUM VALUE OF G_h OR ALPHA_N (UNSTABLE STRAT)
!    *SIGE*     REAL      EPS-EQUATION (DIFFUSION TERM)
!    *SIGK*     REAL      K-EQUATION (DIFFUSION TERM)
!    *TKEMIN*   REAL      MINIMUM VALAUE FOR TURBULENCE KINETIC ENERGY
!                         (ILIM=1)                                      [J/kg]
!    *Z0BOT*    REAL      BOTTOM ROUGHNESS LENGTH SCALE                    [m]
!    *Z0SUR*    REAL      SURFACE ROUGHNESS LENGTH SCALE                   [m]

!-----------------------------------------------------------------------

!*    COMMON *TURCON3* - UNIFORM DIFFUSION COEFFICIENTS

REAL :: cm0, cs0, difmol, heddun, hedvun, vismol

!-----------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *CM0*      REAL      COEFFICIENT IN SMAGORINSKY FORMULATION FOR MOMENTUM
!    *CS0*      REAL      COEFFICIENT IN SMAGORINSKY FORMULATION FOR SCALARS
!    *DIFMOL*   REAL      BACKGROUND OR UNIFORM EDDY DIFFUSIVITY        [m2/s]
!    *HEDDUN*   REAL      UNIFORM HORIZONTAL DIFFUSION COEFFICIENT FOR
!                         SCALARS (IODIF=1)                             [m2/s]
!    *HEDVUN*   REAL      UNIFORM HORIZONTAL DIFFUSION COEFFICIENT FOR
!                         MOMENTUM (IODIF=1)                            [m2/s]
!    *VISMOL*   REAL      BACKGROUND OR UNIFORM EDDY VISCOSITY          [m2/s]

!------------------------------------------------------------------
!*    COMMON *VISC* - EDDY VISCOSITIES, EDDY DIFFUSIVITIES

REAL :: dheddyvc(nr,nc), dheddyvu(nr,nc+1), dheddyvv(nr+1,nc)
REAL :: heddydc(nz,nr,nc), heddydu(nz,nr,nc+1), heddydv(nz,nr+1,nc)
REAL :: heddyp(nr,nc)
REAL :: heddyvc(nz,nr,nc), heddyvu(nz,nr,nc+1), heddyvv(nz,nr+1,nc)
REAL :: veddyd(nz+1,nr,nc), veddye(nz+1,nr,nc)
REAL :: veddyk(nz+1,nr,nc), veddyv(nz+1,nr,nc)



!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *DHEDDYVC* REAL      DEPTH-INTEGRATED HORIZONTAL DIFFUSION COEFFICIENT
!                         FOR MOMENTUM AT CELL CENTRE                   [m2/s]
!    *DHEDDYVU* REAL      DEPTH-INTEGRATED HORIZONTAL DIFFUSION COEFFICIENT
!                         FOR MOMENTUM AT U-NODE                        [m2/s]
!    *DHEDDYVV* REAL      DEPTH-INTEGRATED HORIZONTAL DIFFUSION COEFFICIENT
!                         FOR MOMENTUM AT V-NODE                        [m2/s]
!    *HEDDYDC*  REAL      HORIZONTAL DIFFUSION COEFFICIENT FOR SCALARS
!                         AT CELL CENTRE                                [m2/s]
!    *HEDDYDU*  REAL      HORIZONTAL DIFFUSION COEFFICIENT FOR SCALARS
!                         AT U-NODE                                     [m2/s]
!    *HEDDYDV*  REAL      HORIZONTAL DIFFUSION COEFFICIENT FOR SCALARS
!                         AT V-NODE                                     [m2/s]
!    *HEDDYP*   REAL      HORIZONTAL DIFFUSION COEFFICIENT FOR PARTICLE MODULE
!                         AT CELL CENTRE                                [m2/s]
!    *HEDDYVC*  REAL      HORIZONTAL DIFFUSION COEFFICIENT FOR MOMENTUM
!                         AT CELL CENTRE                                [m2/s]
!    *HEDDYVU*  REAL      HORIZONTAL DIFFUSION COEFFICIENT FOR MOMENTUM
!                         AT U-NODE                                     [m2/s]
!    *HEDDYVV*  REAL      HORIZONTAL DIFFUSION COEFFICIENT FOR MOMENTUM
!                         AT V-NODE                                     [m2/s]
!    *VEDDYD*   REAL      VERTICAL EDDY DIFFUSIVITY (W-NODE)           [m2/s]
!    *VEDDYE*   REAL      VERTICAL DIFFUSION COEFFIFIENT FOR DISSW (W-NODE)
!                                                                       [m2/s]
!    *VEDDYK*   REAL      VERTICAL DIFFUSION COEFFICIENT FOR T.K.E. (W-NODE)
!                                                                       [m2/s]
!    *VEDDYV*   REAL      VERTICAL EDDY VISCOSITY (W-NODE)             [m2/s]

!----------------------------------------------------------------
!* COMMON *WORKSP* - EMPTY ARRAYS USED AS DUMMY ARGUMENTS
!                         FOR SUBROUTINE CALLS

REAL :: zeros1(nr,nc), zeros2(nr,nc), zeros3(nr,nc)
REAL :: zeros4(nr,nc), zeros5(nr,nc), zeros6(nr,nc)

!------------------------------------------------------------------------

!*    NAME      TYPE      PURPOSE
!     ----      ----      -------
!    *ZEROS1*   REAL      ZERO VALUE FOR "SFLUX" (TRANSP ROUTINES)
!    *ZEROS2*   REAL      ZERO VALUE FOR "STRAN" (TRANSP ROUTINES)
!    *ZEROS3*   REAL      ZERO VALUE FOR "SPHI"  (TRANSP ROUTINES)
!    *ZEROS4*   REAL      ZERO VALUE FOR "BFLUX" (TRANSP ROUTINES)
!    *ZEROS5*   REAL      ZERO VALUE FOR "BTRAN" (TRANSP ROUTINES)
!    *ZEROS6*   REAL      ZERO VALUE FOR "BPHI"  (TRANSP ROUTINES)
!-------------------------------------------------------------------

END MODULE param
