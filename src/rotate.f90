module rancom
implicit none
save
real(kind=8) :: RANLST(100)
integer(kind=8) :: ISEED3(8),IBFCTR
contains

SUBROUTINE RANDST(ISEED)
real :: diff
integer(kind=8) :: is,iseed
integer :: i
IBFCTR = 256
IS=ISEED
DO I=1,8
   ISEED3(I)=MOD(IS,IBFCTR)
   IS=IS/IBFCTR
enddo
DO I=1,8
   ISEED3(I) = ISEED3(I) + ISEED3(I)
enddo
DO I=1,8
   diff=ISEED3(I)-IBFCTR
   do while(diff>=0) 
      ISEED3(I) = ISEED3(I) - IBFCTR
      ISEED3(I+1) = ISEED3(I+1) + 1
      diff=ISEED3(I)-IBFCTR
   enddo
enddo
ISEED3(1) = ISEED3(1) + 1
DO I=1,8
   diff=ISEED3(I)-IBFCTR
   do while(diff>=0)
      ISEED3(I) = ISEED3(I) - IBFCTR
      ISEED3(I+1) = ISEED3(I+1) + 1
      diff=ISEED3(I)-IBFCTR
   enddo
enddo
DO I=1,100
   RANLST(I)=RAND1(ISEED3)
enddo
end SUBROUTINE RANDST

FUNCTION RAND0(IDUM) 
integer :: idum,j
real(kind=8) :: rand0
J=INT(99E0*RANLST(100))+1
RAND0 = RANLST(100)
RANLST(100)=RANLST(J)
RANLST(J)=RAND1(ISEED3)
END FUNCTION RAND0

FUNCTION RAND1(ISEED)
integer(kind=8) :: ip
integer :: i,j,k
real :: bi,rand1
integer(kind=8) :: iseed(8),IA(8),IC(8),ID(8)
DATA IA/45,127,149,76,45,244,81,88/
DATA BI/3.90625D-3/

ID=0
IC=0
big: DO J=1,8
   middle: DO I=1,9-J
      K=J+I-1
      IP=IA(J)*ISEED(I)
      do while(k<=8)
         IP=IP+ID(K)
         ID(K)=MOD(IP,IBFCTR)
         IP=IP/IBFCTR
         IF (IP==0) cycle middle
         k=k+1
       enddo
   enddo middle
enddo big
iseed=id

RAND1=FLOAT(ISEED(1))
DO I=2,8
   RAND1=FLOAT(ISEED(I))+RAND1*BI
enddo
RAND1=RAND1*BI
END FUNCTION RAND1

end module rancom

program rotate
use constants
use atsymb
use rancom
implicit none
! program to rotate a molecule about its euler angles
real(kind=8) :: rand
integer(kind=8) :: iclock
real, dimension(:), allocatable :: rr
real :: twopi,rnd,phi,csthta,chi,thta,snthta,csphi,snchi,cschi,snphi
real :: rxx, rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz,wt,x,y,z ! xcm1,ycm1,zcm1,xcm2,ycm2,zcm2 removed for generalization
real :: dist,distm,xdum,ydum,zdum,rmin
integer :: i,numberoffragments,n,j,iii,nc,ijk,ii ! ii extra variable needed for generalization
integer,allocatable :: numberofatomsineachfragment(:),com(:) ! n1 , at1, at2 removed for generalization
real,allocatable :: xcm(:),ycm(:),zcm(:) ! Added for generalization
integer :: atomindex
integer, allocatable:: atomindexstart(:),atomindexend(:)
integer :: interactions
real,dimension(:),allocatable :: q,qq,w
real, dimension(3) :: vcm,qcm


character*2,dimension(:),allocatable :: fasymb


! Added by Swaroop for generalization to more than two monomers
read(*,*) numberoffragments,n,dist,distm

allocate(atomindexstart(numberoffragments),atomindexend(numberoffragments))
allocate(numberofatomsineachfragment(numberoffragments),com(numberoffragments))
allocate(xcm(numberoffragments),ycm(numberoffragments),zcm(numberoffragments))

atomindexstart(1)=1
read(*,*) numberofatomsineachfragment(1),com(1)
atomindexend(1)=numberofatomsineachfragment(1)
atomindex=numberofatomsineachfragment(1)


do i=2,numberoffragments
   read(*,*) numberofatomsineachfragment(i), com(i)
   atomindexstart(i)=atomindex+1
   atomindex=atomindex+numberofatomsineachfragment(i)
   atomindexend(i)=atomindex
enddo


allocate(q(3*n),qq(3*n),fasymb(n),w(n))


interactions=0
do i=1,numberoffragments-1
   interactions=interactions+(numberofatomsineachfragment(i)*(n-atomindexend(i)))
enddo


allocate(rr(interactions))

do i=1,n
   read(*,*) fasymb(i),q(3*i-2),q(3*i-1),q(3*i)
   do j = 1 , natom
      if(fasymb(i)==asymb(j)) w(i)=ams(j)
   enddo
enddo


! Generalized calculation of com of all molecules

do ii=1,numberoffragments
   xcm(ii)=0
   ycm(ii)=0
   zcm(ii)=0
   do i=atomindexstart(ii),atomindexend(ii)
      xcm(ii)=xcm(ii)+w(i)*q(3*i-2) 
      ycm(ii)=ycm(ii)+w(i)*q(3*i-1) 
      zcm(ii)=zcm(ii)+w(i)*q(3*i)
      wt=wt+w(i)
   enddo
   xcm(ii)=xcm(ii)/wt
   ycm(ii)=ycm(ii)/wt
   zcm(ii)=zcm(ii)/wt
   do i=atomindexstart(ii),atomindexend(ii)
      if(com(ii)==-1) then
         qq(3*i-2)=q(3*i-2)-xcm(ii)
         qq(3*i-1)=q(3*i-1)-ycm(ii)
         qq(3*i)=q(3*i)-zcm(ii)
      else
         qq(3*i-2)=q(3*i-2)-q(3*com(ii)-2)
         qq(3*i-1)=q(3*i-1)-q(3*com(ii)-1)
         qq(3*i)=q(3*i)-q(3*com(ii))
      endif
   enddo
enddo

twopi=2*pi

! Generalized random rotation of all molecules

CALL SYSTEM_CLOCK(COUNT=iclock)
CALL RANDST(iclock)

do iii=1,100000000
   nc=0
   do ii=1,numberoffragments
      rand=rand0(0)
      PHI=TWOPI*RAND
      rand=rand0(0)
      CSTHTA=2.0D0*RAND-1.0D0
      rand=rand0(0)
      CHI=TWOPI*RAND
      THTA=ACOS(CSTHTA)
      SNTHTA=SIN(THTA)
      SNPHI=SIN(PHI)
      CSPHI=COS(PHI)
      SNCHI=SIN(CHI)
      CSCHI=COS(CHI)
      RXX=CSTHTA*CSPHI*CSCHI-SNPHI*SNCHI
      RXY=-CSTHTA*CSPHI*SNCHI-SNPHI*CSCHI
      RXZ=SNTHTA*CSPHI
      RYX=CSTHTA*SNPHI*CSCHI+CSPHI*SNCHI
      RYY=-CSTHTA*SNPHI*SNCHI+CSPHI*CSCHI
      RYZ=SNTHTA*SNPHI
      RZX=-SNTHTA*CSCHI
      RZY=SNTHTA*SNCHI
      RZZ=CSTHTA
      do i=atomindexstart(ii),atomindexend(ii)
         x=QQ(3*i-2)*RXX+QQ(3*i-1)*RXY+QQ(3*i)*RXZ
         y=QQ(3*i-2)*RYX+QQ(3*i-1)*RYY+QQ(3*i)*RYZ
         z=QQ(3*i-2)*RZX+QQ(3*i-1)*RZY+QQ(3*i)*RZZ
         if(ii==1)then
            q(3*i-2)=x
            q(3*i-1)=y
            q(3*i)=z
         endif
         if(ii>1) then
            q(3*i-2)=x+dist
            q(3*i-1)=y
            q(3*i)=z
            do ijk=atomindexend(ii-1),1,-1
               nc=nc+1
               xdum=q(3*i-2)-q(3*ijk-2)
               ydum=q(3*i-1)-q(3*ijk-1)
               zdum=q(3*i)-q(3*ijk)
               rr(nc)=sqrt(xdum*xdum + ydum*ydum + zdum*zdum)
            enddo
         endif
      enddo
   enddo
   
   rmin=minval(rr)
   if(rmin>distm) exit
enddo


do i=1,n
   print*, fasymb(i),q(3*i-2),q(3*i-1),q(3*i)
enddo

end program rotate 
