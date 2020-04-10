PROGRAM bath

INTEGER, PARAMETER :: nx = 50	
INTEGER, PARAMETER :: ny = 40	
REAL :: PI
REAL :: DX,DY,XX,YY
REAL :: BB, YZERO, XB
REAL :: HC, DH, AA, XZERO
REAL :: FF, FAC, YRE
REAL :: H(ny,nx)
REAL :: dH1(ny,nx),dH2(ny,nx),diff
REAL :: axis,dist,adjust, distx,disty
INTEGER :: J,K

dx = 1000.0
dy = 1000.0

PI = 4.0*ATAN(1.0)
!
! MAIN SLOPE 
!

do i = 1,nx
do j = 1,ny
   H(j,i) = 1000.0-900.0*REAL(j-1)/REAL(ny-1)
end do
end do

!
!
!do j = 35,ny
!  H(j,i) = 200.0
!end do
!end do

do i = 1,nx
 H(ny,i) = 0.0
end do


do j = 1,ny
 a1 = min(max(j-5,0)/8.0,1.0)
 a2 = min(max(34-j,0)/5.0,1.0)
 i = 35
 H(j,i) =  min(H(j,i)+a1*a2*200.0,1000.0)
end do

DO K = 1,NX
DO J = 1,ny
  dh1(j,k) = 0.0
  dh2(j,k) = 0.0
END DO
END DO

! add canyon prototype

!example 1: normal canyon
!
! yo = from 5 to 30
! xo = 17 is constant +-4

!DO K = 1,NX
!DO J = 5,30
!  dist = ABS(REAL(k-25)/2.0) ! distance to the axis
!  if(dist<1.0)THEN
!   dh1(j,k) = 500.0 
!  end if
!  adjust = MIN(REAL(ny-10-j)/5.0,1.0) ! adjustment
!  dh1(j,k) = adjust*dh1(j,k) 
!END DO
!END DO

!DO K = 12,25
!DO J = 1,NY
!  dist = ABS(REAL(j-22)/2.0) ! distance to the axis
!  if(dist<1.0)THEN
!   dh2(j,k) = 500.0 
!  end if
!  adjust = MIN(REAL(k-11)/4.0,1.0) ! adjustment
!  dh2(j,k) = adjust*dh2(j,k) 
!END DO
!END DO

!example 2: inclined canyon via rotation
! CENTRAL point of rotation
!J0 = 30
!K0 = 75
!PI = 4.*atan(1.)
!       write(6,*)'Give rotation angle =>'
!       read(5,*)alpha

!alpha = -30.
!beta = alpha*PI/180.

!dh2(1,1) = dh1(J0,K0)

!DO JJ = 1,NY
!DO KK = 1,NX
!  XX = real(KK-K0)
!  YY = real(JJ-J0)
! new location of grid centre after rotation
!  x = XX*cos(beta)+YY*sin(beta)
!  y = -xx*sin(beta)+YY*cos(beta)
! index for interpolation
!  k = int(x)
!  dx = x-real(k)
!  if(dx.GT.0.)then
!     kr = k+1
!  else
!     kr = k-1
!  end if     
!  j = int(y)
!  dy = y-real(j)
!  if(dy.GT.0.)then
!    jr = j+1
!  else
!    jr = j-1
!  end if     
!! interpolation
!  dx = abs(dx)
!  dxr = 1.-dx
!  dy = abs(dy)
!  dyr = 1.-dy
!
!  j = max(j+j0,1)
!  j = min(j,ny)
!  jr = max(jr+j0,1)
!  jr = min(jr,ny) 
!  k = max(k+k0,1)
!  k = min(k,nx)
!  kr = max(kr+k0,1)
!  kr = min(kr,nx)
!  dh2(JJ,KK) = dx*dy*dh1(j,k)+dx*dyr*dh1(jr,k)+dxr*dy*dh1(j,kr)+dxr*dyr*dh1(jr,kr)   
!end do
!end do

!diffusion

DO j = 1,ny
DO k = 1,nx
 dh1(j,k) = H(j,k)
END DO
END DO


!DIFFUSION

dh = 0.0002

DO n = 1,2000

DO j = 2,ny-1
DO k = 2,nx-1

hc = dh1(j,k)
dh2(j,k) = hc

he = dh1(j,k+1)-hc
hw = dh1(j,k-1)-hc
hn = dh1(j+1,k)-hc
hs = dh1(j-1,k)-hc

!IF(h(j,k+1)<=50.0)he = 0.0
!IF(h(j,k-1)<=50.0)hw = 0.0
IF(h(j+1,k)<=50.0)hn = 0.0
!IF(h(j-1,k)<=50.0)hs = 0.0

!IF(hc>50.0)THEN
 dh2(j,k) = hc+dh*(he+hw+hn+hs)
!END IF

END DO
END DO  

!boundaries
DO k = 1,nx
 dh2(1,k) =  dh2(2,k)
 dh2(ny,k) =  dh2(ny,k)
END DO

DO j = 1,ny
 dh2(j,1) =  dh2(j,2)
 dh2(j,nx) =  dh2(j,nx-1)
END DO

!updating
DO j = 1,ny
DO k = 1,nx
 dh1(j,k) = dh2(j,k)
END DO
END DO  

END DO

DO K = 1,NX
DO J = 1,NY
 H(j,k) =  dh2(j,k)
 H(j,k) =  min(H(j,k),3000.0)
END DO
END DO

 
OPEN(10,file='topo1.dat',form='formatted',recl = 100000)                                                        
DO J = 1,NY
  WRITE(10,'(500F18.8)')(H(J,K),K=1,NX)
END DO                                                                             
CLOSE(10)                                                                          
 
END PROGRAM bath