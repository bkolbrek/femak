      SUBROUTINE EXTFLD100(VI100V,VFREQ,VPREE,VDIMP,KLOCE,
     1     VCORE,VPRNE,KLD,KNE,VCOEXT)
C**********************************************************************C
C     EXTF100.FOR						       C
C     Tonni F. Johansen 					       C
C     Kompaniet Numerisk Lyd					       C
C     27/1/93							       C
C								       C
C**********************************************************************C
C     To calculate external pressure using Rayleighs integral. Element C
C     type 100. As described by Kagawa et.al. J.S.&V. vol.69	       C
C     1980 p.207-228.						       C
C**********************************************************************C

      COMPLEX*16 VI100V(1)
      INTEGER KLD(1),KNE(1),KLOCE(1)
      REAL*8 VFREQ(1),VPREE(1),VCORE(1),VDIMP(1),VPRNE(1),VCOEXT(1)


c      INTEGER NDIM,NNT
c      COMMON/COOR/NDIM,NNT

c      INTEGER NGIMP,NGSRC,NFREQ,NSIMP,NSRC,PRPAR,VELPAR,INTPAR,NPEXT
c      COMMON/AKUS/NGIMP,NGSRC,NFREQ,NSIMP,NSRC,PRPAR,VELPAR,INTPAR,NPEXT

c      INTEGER IFL,IGIMP,IFREQ,IGSRC,IFL1,NI100NOD,NI100MOD,NI100FL
c      REAL*8 I100RAD
c      COMMON/IMPEL/IFL,IGIMP,IFREQ,IGSRC,IFL1,NI100NOD,NI100MOD,
c     1		   NI100FL,I100RAD

c      INTEGER IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
c     &		  IPG,ICODE,IDLE0,INEL0,IPG0
c      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
c     &		  IPG,ICODE,IDLE0,INEL0,IPG0


c      INTEGER M,MR,MP,M1,M2,M3
c      COMMON/ES/M,MR,MP,M1,M2,M3

      include 'femak.cmn'

      COMPLEX*16 CZERO,VGAUSS,GREEN,DELTAP,P
      INTEGER I,J,IJ,IG,NODE1,NODE2,NODE3,NUMTHETA
      REAL*8 VNI(18),VJ(2)
      REAL*8 K,Y,C,KSI(3),VCPG(3),DELTA,THETA
      REAL*8 XM(3),X0(3),NM(3),NP(3),ZERO,PI,DETJ,COEF
      REAL*8 AMPL,DBAMPL,PHASE

      DATA KSI/-.774596669241483D0,0.0D0,.774596669241483D0/
      DATA VCPG/0.555555555555556D0,0.888888888888889D0,
     1		0.555555555555556D0/

C------------Number of Gauss points in ksi,eta,dzeta direction
C---------------------------------------------------
      DATA ZERO/0.0D0/,CZERO/(0.0,0.0)/
      DATA PI/3.141592654/
      DATA NM/0.0,0.0,0.0/,NP/0.0,0.0,0.0/

      W=2*PI*VFREQ(IFREQ)
      K=W/VPREE(2)

C.........Get elementinformation.
      ITPE=0
      REWIND(M6)
      DO WHILE (ITPE.NE.100)
          CALL RDCIMP(M6,KLOCE,VCORE,VPRNE,VPREE,KNE,IGIMP)
      END DO
      IPG=3

      WRITE(MP,1990)
1990  FORMAT(' COMPUTATIONS AT EXTERNAL FIELDPOINTS:',/72('='))

      WRITE(MP,1995)
1995  FORMAT(/1X' POINT NO.',3X,'X',11X,'Y',11X,'Z',18X,'COMPLEX',
     1	     12X,'AMPL',8X,'dB',10X,'PHASE (RAD)'/)
C............Calculate typical size of "element"
      DELTA=I100RAD/(NI100NOD-1)

C...........FORMFUNCITONS AND DERIVATIVES
      DO I=1,3
	VNI(6*(I-1)+1)=-0.5*KSI(I)*(1-KSI(I))
	VNI(6*(I-1)+2)=1-KSI(I)*KSI(I)
	VNI(6*(I-1)+3)=0.5*KSI(I)*(1+KSI(I))
	VNI(6*(I-1)+4)=-0.5*(1-2*KSI(I))
	VNI(6*(I-1)+5)=-2*KSI(I)
	VNI(6*(I-1)+6)=0.5*(1+2*KSI(I))
      ENDDO

C............Loop through all external fieldpoints
      DO IEXT=1,NPEXT
      XM(1)=VCOEXT(IEXT*2-1)
      XM(2)=VCOEXT(IEXT*2)
      XM(3)=0.0
      P=DCMPLX(0.0,0.0)

C...........Loop through all elementsides in hole
      DO I=1,NI100FL
	NODE1=(I-1)*2+1
	NODE2=(I-1)*2+2
	NODE3=(I-1)*2+3

Cnb!! KOPIERT INN FRA ELEM03.FOR STARTER HER

	DO IG=1,IPG

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

C.............Compute velocity in gauss point
	  VGAUSS=VNI((IG-1)*6+1)*VI100V(NODE1)
	  VGAUSS=VGAUSS+VNI((IG-1)*6+2)*VI100V(NODE2)
	  VGAUSS=VGAUSS+VNI((IG-1)*6+3)*VI100V(NODE3)

	  Y=ZERO
	  DO J=1,3
	    Y=Y+VNI((IG-1)*6+J)*VCORE((I-1)*4+2*J)
	  END DO
	  NUMTHETA=MAX0(6,INT(2.*PI*Y/DELTA))

	  COEF=COEF*Y*2*PI/NUMTHETA

C.............Evaluate integral in theta direction
	  DELTAP=CZERO
	  DO IJ=1,NUMTHETA
	    THETA=2*PI*(IJ-1)/NUMTHETA
	    X0(1)=VCORE(1)
	    X0(2)=Y*COS(THETA)
	    X0(3)=Y*SIN(THETA)
	    CALL GREEN3D(XM,X0,NM,NP,K,0,GREEN)
	    DELTAP=DELTAP+2*GREEN
	  END DO

	  P=P+DELTAP*VGAUSS*COEF
	ENDDO

C..nb!! KOPIERT INN FRA ELEM03.FOR SLUTTER HER



      ENDDO
C...........Finished all element sides in circular hole
      P=P*DCMPLX(0.,W*VPREE(1))
      AMPL=SQRT(REAL(P)*REAL(P)+IMAG(P)*IMAG(P))
      DBAMPL=20*LOG10(AMPL/1.0)
      PHASE=ATAN2(IMAG(P),REAL(P))
      WRITE(MP,2000)IEXT,XM(1),XM(2),XM(3),P,AMPL,DBAMPL,PHASE
2000  FORMAT(1X,I5,3E12.5,5X,5E12.5)


      END DO
C............Finished loop through all external fieldpoints


      RETURN
      END 


      SUBROUTINE GREEN3D(XM,XP,NM,NP,XK,GTYPE,GREEN)
C*********************************************************************C
C 8/2-1993 Tonni F. Johansen					       C
C Denne subroutinen beregner "frittfelts" greens funksjon. Det        C
C beregnes normal deriverte eller den dobbelt normalderiverte av      C
C greensfunksjonen.                                                   C
C Med parameteren GTYP bestemmes om det skal finnes  dg(x,x0)/dn0,    C
C dg(x,x0)/dn eller d2g(x,x0)/dndn0                                   C
C INPUT :                                                             C
C -------                                                             C
C     XM    -  Koordinater for "maalepunkt" , REAL XM(3)              C
C     XP    -  Koordinater for "kildepunkt" , REAL XP(3)              C
C     NM    -  Normalen i "maalepunktet",     REAL NM(3)              C
C     NP    -  Normalen i "kildepunktet",     REAL NP(3)              C
C     XK    -  Boelgetallet.                  REAL XM                 C
C     GTYPE -  Bestemmer hvilken type normalederivert :               C
C               GTYPE = 0, beregner g                                 C
C	 >     GTYPE = 1, beregner dg/dn0			      C
C  IKKE	 >     GTYPE = 2, beregner dg/dn			      C
C  MED	 >     GTYPE = 3, beregner d2g/dndn0			      C
C OUTPUT :                                                            C
C --------                                                            C
C     GREEN - Verdien av den greensfunksjonen  ,COMPLEX               C
C                                                                     C
C*********************************************************************C
C NB!!	Denne versjonen bruker j (=-i) som imagin‘r enhet.	      C
C  pr.8/2-1993 er kun greens funksjonen implementert		      C
C*********************************************************************C

      REAL*8 PI
      PARAMETER(PI=3.141592654)

      REAL*8 XM(3),XP(3),NM(3),NP(3),XK
      INTEGER GTYPE
      COMPLEX*16 GREEN
      REAL*8 RV(3),R

C.........Nullstiller Greensfunksjonen
      GREEN=CMPLX(0.,0.)

C..........Finner avstanden mellom maalepunkt og kildepunktene, R(I)
      R=SQRT((XM(1)-XP(1))*(XM(1)-XP(1))+(XM(2)-XP(2))*(XM(2)-XP(2))+
     &       (XM(3)-XP(3))*(XM(3)-XP(3)))


C...............Finner greens funksjon
      GREEN=1/(4*PI*R)*DCMPLX(COS(XK*R),-SIN(XK*R))



100   RETURN
      END
