!-------------------------------------------------------------------------------
! General Sigma Coordinate Model  - GSCM Version 1.0
!-------------------------------------------------------------------------------
! By:
!   Jochen Kaempf
!     School of Chemistry, Physics and Earth Sciences, Flinders University Adelaide
! 
!-------------------------------------------------------------------------------
PROGRAM RUNA
  !-------------------------------------
  USE param
  USE functions
  USE cohini
  USE cohrun
  INTEGER(4) :: NTTT, NTTOT, NTOUT,NCOUT
  CHARACTER(4) :: OFILE
  INTEGER(4) :: IOUT,JOUT,KOUT,LOUT,k,j,i,n
  REAL :: Q1a,Q1b,Q1c,time,hhh,vvv,vsum,csum
  REAL :: Q2a,Q2b,Q2c,Q3a,Q3b,Q3c
  REAL :: dummy
! LAGRANGIAN MODEL JK
  REAL :: zposreal(ntrac),xposI(ntrac),yposI(ntrac)
   REAL :: depp,DX2,DY2,epsx,epsy,term1,term2,term3,term4
  INTEGER(4) :: NTRA,II,JJ,WEST,SOUTH

IOUT = 0 
UGEO = 0.0

TRIGGER = 1.0

 ! ********** INITILIZE OCEAN MODEL *********

FLOAT_ON = .true.

CALL INITCOH

 ! simulation time

NTTOT = 10*24*3600/delt

 ! output every 3 hours
NTOUT =  3*3600/delt

 ! tracer output every 6 hours

NCOUT = 3*3600/delt

! LAGRANGIAN MODEL JK
 ! float outputs every 1 hour
NTRA = 1*3600/delt

DO n = 1,4

QTNIT1(n) = 550+n
QTNIT2(n) = 650+n
QTNIT3(n) = 750+n

QXNIT1(n) = 580+n
QXNIT2(n) = 680+n

END DO

!========================================================================
! --particle 
OPEN(QTNIT1(1),file = 'dat/Trx.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT2(1),file = 'dat/Try.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT3(1),file = 'dat/Trz.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT1(2),file = 'dat/Trh.dat',form='formatted',status='unknown',recl=1000000)

DO n = 1,ntrac

  II = IT(n)
  JJ = JT(n)
  KK = KT(n)
  xpos(n) = REAL(II-0.5)*gx2(1,1)+XT(n)
  ypos(n) = REAL(JJ-0.5)*gy2(1)+YT(n)
  DX2 = 0.5*GX2(JJ,II)
  DY2 = 0.5*GY2(JJ)
  DZ2 = 0.5*GZSC(KK,JJ,II)

   WEST = 0
  epsx = XT(n)/GX2(JJ,II)
  IF(XT(n)<0.0)THEN
     WEST = -1
     epsx = 1.0+XT(n)/GX2(JJ,II-1)
  END IF

  SOUTH = 0
  epsy = YT(n)/GY2(JJ)
  IF(YT(n)<0.0)THEN
     SOUTH = -1
     epsy = 1.0+YT(n)/GY2(JJ-1)
  END IF

  term1 = h2atc(JJ+SOUTH,II+WEST)*(1.0-epsx)*(1.0-epsy)
  term2 = h2atc(JJ+SOUTH,II+WEST+1)*epsx*(1.0-epsy)
  term3 = h2atc(JJ+SOUTH+1,II+WEST)*(1.0-epsx)*epsy
  term4 = h2atc(JJ+SOUTH+1,II+WEST+1)*epsx*epsy
  depp = term1+term2+term3+term4
  zpos(n) = (1.0-(gz0(KK,JJ,II)-DZ2+ZT(n)) )*depp
  zposreal(n) = zpos(n)
  posH(n) = depp
END DO

WRITE(QTNIT1(1),*)(xpos(n)/1000.0,n=1,ntrac)
WRITE(QTNIT2(1),*)(ypos(n)/1000.0,n=1,ntrac)
WRITE(QTNIT3(1),*)(zposreal(n),n=1,ntrac)
WRITE(QTNIT1(2),*)(posH(n),n=1,ntrac)

fcount = 0
DO  i=1,nc
DO  j=1,nr
  zetafinal(j,i) = 0.0
  DO k=1,nz
    ufinal(k,j,i) = 0.0
    vfinal(k,j,i) = 0.0
    wfinal(k,j,i) = 0.0
    sfinal(k,j,i) = 0.0
    tfinal(k,j,i) = 0.0
    rfinal(k,j,i) = 0.0
    pfinal(k,j,i) = 0.0
    hedfinal(k,j,i) = 0.0
    vedfinal(k,j,i) = 0.0
  END DO
  wfinal(nz+1,j,i) = 0.0
END DO
END DO

 ! ********** TIME ITERATION *********

  DO NTTT = 1, NTTOT

time = time + delt

! call dynamics module

  CALL RUNCOH(NTTT)

IF (MOD(NTTT,NTRA)==0) THEN 


IF(ttime>=2.0)THEN

DO n = 1,ntrac

  II = IT(n)
  JJ = JT(n)
  KK = KT(n)
  xpos(n) = REAL(II-0.5)*gx2(1,1)+XT(n)
  ypos(n) = REAL(JJ-0.5)*gy2(1)+YT(n)
  DX2 = 0.5*GX2(JJ,II)
  DY2 = 0.5*GY2(JJ)
  DZ2 = 0.5*GZSC(KK,JJ,II)

   WEST = 0
  epsx = XT(n)/GX2(JJ,II)
  IF(XT(n)<0.0)THEN
     WEST = -1
     epsx = 1.0+XT(n)/GX2(JJ,II-1)
  END IF

  SOUTH = 0
  epsy = YT(n)/GY2(JJ)
  IF(YT(n)<0.0)THEN
     SOUTH = -1
     epsy = 1.0+YT(n)/GY2(JJ-1)
  END IF

  term1 = h2atc(JJ+SOUTH,II+WEST)*(1.0-epsx)*(1.0-epsy)
  term2 = h2atc(JJ+SOUTH,II+WEST+1)*epsx*(1.0-epsy)
  term3 = h2atc(JJ+SOUTH+1,II+WEST)*(1.0-epsx)*epsy
  term4 = h2atc(JJ+SOUTH+1,II+WEST+1)*epsx*epsy
  depp = term1+term2+term3+term4
  zpos(n) = (1.0- (gz0(KK,JJ,II)-DZ2+ZT(n)) )*depp
  zposreal(n) = zpos(n)
  posH(n) = depp
END DO

WRITE(QTNIT1(1),*)(xpos(n)/1000.0,n=1,ntrac)
WRITE(QTNIT2(1),*)(ypos(n)/1000.0,n=1,ntrac)
WRITE(QTNIT3(1),*)(zposreal(n),n=1,ntrac)
WRITE(QTNIT1(2),*)(posH(n),n=1,ntrac)

END IF


END IF

! output of spatial distribution arrays

IF (MOD(NTTT,NTOUT)==0) THEN 
   IOUT = IOUT + 1
   jout = iout/10
   kout = iout-jout*10
   OFILE = '.000'
   IF (iout<10) THEN
     WRITE(OFILE(4:4),'(I1.1)')iout
   ELSE
   IF (iout<100) THEN
      WRITE(OFILE(3:3),'(I1.1)')jout
      WRITE(OFILE(4:4),'(I1.1)')kout
   ENDIF
   END IF  
   IF (iout>=100) THEN
   lout = iout/100
   jout = (iout-lout*100)/10
   kout = (iout-lout*100-jout*10)
   WRITE(OFILE(2:2),'(I1.1)')lout       
   WRITE(OFILE(3:3),'(I1.1)')jout
   WRITE(OFILE(4:4),'(I1.1)')kout
END IF
  
   !      WRITE(*,*)'OUTPUT TO FILE =',OFILE
 
         open(unit=97,file='dat/us'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,'(250F16.8)')(u2(nz,j,i)+ugeo,i=1,nc)
         END DO
         close(unit=97)

         open(unit=97,file='dat/vs'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,'(250F16.8)')(v2(nz,j,i),i=1,nc)
         END DO
         close(unit=97)

         open(unit=97,file='dat/sb'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,FMT='(250F16.8)')(s(1,j,i),i=1,nc)
         END DO
         close(unit=97)

         open(unit=97,file='dat/tb'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,FMT='(250F16.8)')(t(1,j,i),i=1,nc)
         END DO
         close(unit=97)

         open(unit=97,file='dat/rb'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,FMT='(250F16.8)')(ro(1,j,i)-r0ref,i=1,nc)
         END DO
         close(unit=97)

         open(unit=97,file='dat/pb'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,FMT='(250F16.8)')(rpress(1,j,i),i=1,nc)
         END DO
         close(unit=97)

         open(unit=97,file='dat/ss'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,FMT='(250F16.8)')(s(nz,j,i),i=1,nc)
         END DO
         close(unit=97)

         open(unit=97,file='dat/ts'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,FMT='(250F16.8)')(t(nz,j,i),i=1,nc)
         END DO
         close(unit=97)

         open(unit=97,file='dat/rs'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,FMT='(250F16.8)')(ro(nz,j,i)-r0ref,i=1,nc)
         END DO
         close(unit=97)

         open(unit=97,file='dat/ub'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,'(250F16.8)')(u2(1,j,i)+ugeo,i=1,nc)
         END DO
         close(unit=97)

         open(unit=97,file='dat/vb'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,'(250F16.8)')(v2(1,j,i),i=1,nc)
         END DO
         close(unit=97)

         open(unit=97,file='dat/eta'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,'(250F16.8)')(zeta2(j,i),i=1,nc)
         END DO
         close(unit=97)

! transects

         open(unit=97,file='dat/ttrans1'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(t(k,j,10),j=1,nr)
         END DO
         close(unit=97)

         open(unit=97,file='dat/strans1'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(s(k,j,10),j=1,nr)
         END DO
         close(unit=97)

         open(unit=97,file='dat/utrans1'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(u2(k,j,10)+ugeo,j=1,nr)
         END DO
         close(unit=97)
		 
	 open(unit=97,file='dat/vtrans1'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(v2(k,j,10),j=1,nr)
         END DO
         close(unit=97)

	 open(unit=97,file='dat/rtrans1'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(ro(k,j,10)-r0ref,j=1,nr)
         END DO
         close(unit=97)

         open(unit=97,file='dat/ttrans2'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(t(k,j,20),j=1,nr)
         END DO
         close(unit=97)

         open(unit=97,file='dat/strans2'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(s(k,j,20),j=1,nr)
         END DO
         close(unit=97)

         open(unit=97,file='dat/utrans2'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(u2(k,j,20),j=1,nr)
         END DO
         close(unit=97)
		 
	 open(unit=97,file='dat/vtrans2'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(v2(k,j,20),j=1,nr)
         END DO
         close(unit=97)

	 open(unit=97,file='dat/rtrans2'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(ro(k,j,20)-r0ref,j=1,nr)
         END DO
         close(unit=97)

         open(unit=97,file='dat/ttrans3'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(t(k,j,30),j=1,nr)
         END DO
         close(unit=97)

         open(unit=97,file='dat/strans3'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(s(k,j,30),j=1,nr)
         END DO
         close(unit=97)

         open(unit=97,file='dat/utrans3'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(u2(k,j,30),j=1,nr)
         END DO
         close(unit=97)
		 
	 open(unit=97,file='dat/vtrans3'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(v2(k,j,30),j=1,nr)
         END DO
         close(unit=97)

	 open(unit=97,file='dat/rtrans3'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(ro(k,j,30)-r0ref,j=1,nr)
         END DO
         close(unit=97)


         open(unit=97,file='dat/ttrans4'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(t(k,j,20),j=1,nr)
         END DO
         close(unit=97)

         open(unit=97,file='dat/strans4'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(s(k,j,20),j=1,nr)
         END DO
         close(unit=97)

         open(unit=97,file='dat/utrans4'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(u2(k,j,20),j=1,nr)
         END DO
         close(unit=97)
		 
	 open(unit=97,file='dat/vtrans4'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(v2(k,j,20),j=1,nr)
         END DO
         close(unit=97)

	 open(unit=97,file='dat/rtrans4'//OFILE,form='formatted',status='unknown')
         DO k = 1,nz
         write(97,'(250F16.8)')(ro(k,j,20)-r0ref,j=1,nr)
         END DO
         close(unit=97)

      ENDIF
  
  END DO

close(unit=33)

DO  i=1,nc
DO  j=1,nr
  zetafinal(j,i) =  zetafinal(j,i)/fcount
  DO k=1,nz
    ufinal(k,j,i) = ufinal(k,j,i)/fcount
    vfinal(k,j,i) = vfinal(k,j,i)/fcount
    wfinal(k,j,i) = wfinal(k,j,i)/fcount
    tfinal(k,j,i) = tfinal(k,j,i)/fcount
    sfinal(k,j,i) = sfinal(k,j,i)/fcount
    rfinal(k,j,i) = rfinal(k,j,i)/fcount
    pfinal(k,j,i) = pfinal(k,j,i)/fcount
    hedfinal(k,j,i) = hedfinal(k,j,i)/fcount
    vedfinal(k,j,i) = vedfinal(k,j,i)/fcount
  END DO
  wfinal(nz+1,j,i) = wfinal(nz+1,j,i)/fcount
END DO
END DO

open(unit=97,file='dat/final.dat',form='unformatted',status='unknown')
write(97)ufinal,vfinal,wfinal,hedfinal,vedfinal,tfinal,sfinal,rfinal,pfinal,zetafinal
close(unit=97)


END PROGRAM RUNA
!-------------------------------------------------------------------------------

