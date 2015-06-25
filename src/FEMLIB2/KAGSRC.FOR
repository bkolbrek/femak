	SUBROUTINE KAGSRC(KLOCE,VCORE,VPRNE,VPREE,KNE,VI100Q,
     1		      MODNR,NODNR)
C**********************************************************************C
C     KAGSRC.FOR						       C
C     Tonni F. Johansen 					       C
C     Kompaniet Numerisk Lyd					       C
C     29/11/92							       C
C								       C
C**********************************************************************C
C     To assemble matrix Q as described by Kagawa et.al.	       C
C     J.S.&V. vol.69  1980 p.207-228.				       C
C**********************************************************************C
      EXTERNAL BESSJ0
      DOUBLE PRECISION BESSJ0

      INTEGER MODNR,NODNR
      COMPLEX*16 VI100Q(NODNR,MODNR)
      REAL*8 VCORE(1),VPREE(1)
      INTEGER KLOCE(1),KNE(1)

c      INTEGER NDIM,NNT,NDLN,NDLT,IAX
c      REAL*8 FAC(3)
c      COMMON/COOR/NDIM,NNT,NDLN,NDLT,FAC,IAX


c      INTEGER IFL,IGIMP,IFREQ,IGSRC,IFL1,NI100NOD,NI100MOD,NI100FL
c      REAL*8 I100RAD
c      COMMON/IMPEL/IFL,IGIMP,IFREQ,IGSRC,IFL1,
c     1		   NI100NOD,NI100MOD,NI100FL,I100RAD


c      INTEGER IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
c     &		  IPG,ICODE,IDLE0,INEL0,IPG0
c      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
c     &		  IPG,ICODE,IDLE0,INEL0,IPG0

c      INTEGER M,MR,MP,M1,M2,M3
c      COMMON/ES/M,MR,MP,M1,M2,M3

      include 'femak.cmn'

      INTEGER NODE1,NODE2,NODE3,I,J,IPG1,INI,IG,I1,IK
C---------------------------------------------------
      REAL*8 ZERO,PI,WM(20),DETJ,COEF,Y

      COMPLEX*16 CZERO
C-----------CHARACTERISTIC DIMENSIONS OF THE ELEMENT
C----------(VALID  UP TO 3 DIMENSIONS)
C     DIMENSION VCPG(IPG),VKPG(NDIM*IPG),
      REAL*8 VCPG(3),VKPG(3)
C     DIMENSION VJ(NDIM*NDIM),VJ1(NDIM*NDIM)
      REAL*8 VJ(2),VJ1(2)
C     DIMENSION VNI((1+NDIM)*INEL*IPG)
      REAL*8 VNI(18)
      REAL*8 VKPG1(27),VCPG1(27)

      INTEGER IPGKED
      REAL*8 KSI(3)
      REAL*8 GAMMA(10)
      REAL*8 J0GAMMA(10)
C...........Zero in first order bessel function Gamma
      DATA GAMMA/0.0,3.8317059702,7.0155866698,10.1734681351,
     1	      13.3236919363,16.4706300509,19.6158585105,
     2	      22.7600843806,25.9036720876,29.0468285349/

C............Value of zeroth order bessel function at gamma
      DATA J0GAMMA/1.0,-0.4027593957,0.3001157525,-0.2497048771,
     1	      0.21835994072,-0.1964653715,0.1800633753,-0.1671846005,
     2	      0.1567249863,-0.1480111100/

      DATA KSI/-.774596669241483D0,0.0D0,.774596669241483D0/
      DATA VCPG/0.555555555555556D0,0.888888888888889D0,
     1		0.555555555555556D0/

C------------Number of Gauss points in ksi,eta,dzeta direction
C---------------------------------------------------
      DATA ZERO/0.0D0/,CZERO/(0.0D0,0.0D0)/
      DATA PI/3.141592654/


      IPG=3
C.............Zero the matrix
      DO I=1,NI100NOD
	DO J=1,NI100MOD
	  VI100Q(I,J)=CZERO
	ENDDO
      ENDDO

C
C-----------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C-----------THEIR DERIVATIVES AT G.P.

C      CALL GAUSS(IPGKED,1,VKPG,VCPG,IPG)
C      CALL NI01(VKPG,VNI)
C...........FORMFUNCITONS AND DERIVATIVES
      DO I=1,3
	VNI(6*(I-1)+1)=-0.5*KSI(I)*(1-KSI(I))
	VNI(6*(I-1)+2)=1-KSI(I)*KSI(I)
	VNI(6*(I-1)+3)=0.5*KSI(I)*(1+KSI(I))
	VNI(6*(I-1)+4)=-0.5*(1-2*KSI(I))
	VNI(6*(I-1)+5)=-2*KSI(I)
	VNI(6*(I-1)+6)=0.5*(1+2*KSI(I))
      ENDDO
C...........Loop through all elementsides in hole
      DO I=1,NI100FL
	NODE1=(I-1)*2+1
	NODE2=(I-1)*2+2
	NODE3=(I-1)*2+3

Cnb!! KOPIERT INN FRA ELEM03.FOR STARTER HER

C..............funksjoner etc.
C	CALL NI01(VKPG1,VNI)


	INI=0
	DO IG=1,IPG
C.............. calculate the jacobian matrix VJ , determinant
C.............. will return as zero ,inverse not calculated
C	  I1=INI+INEL+1
C	  CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)

C...........2-DIMENSIONS, J=1 D(X)/D(KSI),J=2 D(Y)/D(KSI)
	  DO J=1,2
	    C=ZERO
	    DO IJ=1,3
	      C=C+VNI((IG-1)*6+3+IJ)*VCORE((I-1)*4+2*(IJ-1)+J)
	    ENDDO
	    VJ(J)=C
	  ENDDO

	  DETJ=SQRT(VJ(1)*VJ(1)+VJ(2)*VJ(2))
C............. compute the weight
	  COEF=VCPG(IG)*DETJ
	  Y=ZERO
	  DO J=1,3
	    Y=Y+VNI((IG-1)*6+J)*VCORE((I-1)*4+2*J)
	  END DO
	  COEF=COEF*2.*PI*Y

C.............Get eigenfunction at Y
	  DO J=1,NI100MOD
	    WM(J)=BESSJ0(GAMMA(J)/I100RAD*Y)/SQRT(PI)/I100RAD/J0GAMMA(J)
	  ENDDO

C............. terms of the KAGAWA Q-MATRIX
	  IK=0
	  DO J=1,NI100MOD
	      VI100Q(NODE1,J)=VI100Q(NODE1,J)+
     1		      VNI((IG-1)*6+1)*WM(J)*COEF
	      VI100Q(NODE2,J)=VI100Q(NODE2,J)+
     1		      VNI((IG-1)*6+2)*WM(J)*COEF
	      VI100Q(NODE3,J)=VI100Q(NODE3,J)+
     1		      VNI((IG-1)*6+3)*WM(J)*COEF
	  ENDDO
C	  INI=INI+(NDIM+1)*INEL
	ENDDO

C..nb!! KOPIERT INN FRA ELEM03.FOR SLUTTER HER



      ENDDO
C...........Finished all element sides in circular hole

      RETURN
      END 
