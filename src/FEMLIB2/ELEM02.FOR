      SUBROUTINE ELEM02(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1      VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C**********************************************************************C
C ELEM02.FOR                                                           C
C  3/1-1995  NTH                                                       C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  8 nodes quadratic element for 2 dimensional elasticity              C
C  Evaluate element information according to ICODE value               C
C     ICODE.EQ.1  return of parameters                                 C
C     ICODE.EQ.2  evaluate interpolation functions and                 C
C                 numerical integration coeffisients                   C
C     ICODE.EQ.3  stiffnes matrix (VKE)                                C
C     ICODE.EQ.4  tangent matrix (VKE)  .....not written.....NB!!      C
C     ICODE.EQ.5  mass matrix (VKE)                                    C
C     ICODE.EQ.6  residuals                                            C
C     ICODE.EQ.7  element load (VFE) (second member)                   C
C     ICODE.EQ.8  evaluate and print stresses                          C
C  Element properties                                                  C
C     VPREE(1)  Young's modulus                                        C
C     VPREE(2)  Poisson's coeffiesient                                 C
C     VPREE(3)  .EQ.0 Plane stress                                     C
C               .EQ.1 Plane strain                                     C
C     VPREE(4)  specific mass                                          C
C**********************************************************************C
C  Programmed after Dhatt & Touzot, p.194-197                          C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
c  $DECLARE

      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VKE1(1),VFE(1)
      COMPLEX*16 VIE(1),VIE1(1),VNDSRC(1),VPRSG(1),VSRE(1),VDLEC(1)
      INTEGER KTYPRSG(1)

C-----------CHARACTERISTIC DIMENSIONS OF THE ELEMENT
C----------(VALID  UP TO 3 DIMENSIONS)
C     DIMENSION VCPG(IPG),VKPG(NDIM*IPG),VDE1(IMATD**2)
      REAL*8 VCPG(9),VKPG(18),VDE1(9)
C     DIMENSION VBE(IMATB*IDLE),VDE(IMATD**2),VJ(NDIM*NDIM),VJ1(NDIM*NDIM)
      REAL*8 VBE(48),VDE(9),VJ(4),VJ1(4)
C     DIMENSION VNIX(INEL*NDIM),VNI((1+NDIM)*INEL*IPG),IPGKED(NDIM)
      REAL*8 VNIX(16),VNI(216)
      INTEGER IPGKED(2)

      SAVE VCPG,VKPG,VNI
      
      REAL*8 DEUX,X05,RADN,EPS,X,Y,DETJ,C,D,EPSX,EPSY,GAMXY,UN,VN
      REAL*8 C1,C2,C3,SIGX,SIGY,TAUXY,FX,FY,DX,DY,XN,YN,TETA,TAUMAX,SIG1
      REAL*8 SIG2
      INTEGER ZERO,I0,IG,I1,I,IDIM1,IDECL,I2,IDL,J,JJ,J0,II,J1,ID,IN,IN1
      INTEGER I3,IMATD
      INTEGER COL,ROW,PNTCO,PNTCO1,PNTRW1
      
C------------Dimension of D-matrix, Number of Gauss points in ksi,eta direct
      DATA IMATD/3/,IPGKED/3,3/
C---------------------------------------------------
      DATA ZERO/0.0D0/,DEUX/2.D0/,X05/0.5D0/,RADN/.572957795130823D2/
      DATA EPS/1.D-6/

      SQRT(X)=DSQRT(X)
      ATAN2(X,Y)=DATAN2(X,Y)

C-------------CHOOSE FUNCTION TO BE EXECUTED
      GO TO (100,200,300,400,500,600,700,800),ICODE

C
C------------RETURN ELEMENT PARAMETERS IN COMMON /RGDT/
100   IDLE0=16
      INEL0=8
      IPG0=9

      CFLAG=0
      RETURN

C
C-----------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C-----------THEIR DERIVATIVES AT G.P.
200   CALL GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
      IF(M.LT.2) GO TO 220
      WRITE(MP,2000)IPG
2000  FORMAT(/I5,'  GAUSS POINTS'/10X,'VCPG',25X,'VKPG')
      I0=1
      DO 210 IG=1,IPG
	I1=I0+NDIM-1
	WRITE(MP,2010)VCPG(IG),(VKPG(I),I=I0,I1)
210	I0=I0+NDIM
2010  FORMAT(1X,F20.15,5X,3F20.15)

220   CALL NI02(VKPG,VNI)
      IF (M.LT.2) RETURN
      I1=3*INEL*IPG
      WRITE(MP,2020) (VNI(I),I=1,I1)
2020  FORMAT(/'FUNCTIONS N AND DERIVATIVES'/ (1X,8E12.5))
      RETURN

C----------COMPUTE ELEMENT STIFFNESS MATRIX
C

C--------INITIALIZE VKE
300   DO 310 I=1,136
        VKE(I)=ZERO
310   CONTINUE

C--------FORM MATRIX D
      CALL D02(VPREE,VDE)
      IF (M.GE.2) WRITE(MP,2030) (VDE(I),I=1,9)
2030  FORMAT(/' MATRIX D'/1X,9E12.5)

C----------LOOP OVER G.P
      I1=1+INEL
      DO 330 IG=1,IPG
C-----------EVALUATE THE JACOBIAN, ITS INVERSE AND ITS DETERMINANT
	CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
	IF (DETJ.LT.EPS) WRITE(MP,2040)IEL,IG,DETJ
2040	FORMAT(' *** ELEM ',I5,' G.P. ',I3,' DET(J)=',E12.5)
	IF (M.GE.2) WRITE(MP,2050)VJ,VJ1,DETJ
2050	FORMAT(/' JACOBIAN=',4E12.5 / ' J INVERS=',
     1	       4E12.5/' DETJ=',E12.5)
C-----------PERFORM D*COEF
	C=VCPG(IG)*DETJ
	DO 320 I=1,9
320	  VDE1(I)=VDE(I)*C
C-----------FORM MATRIX B
	CALL DNIDX(VNI(I1),VJ1,NDIM,INEL,VNIX)
	IF (M.GE.2) WRITE(MP,2060)(VNIX(I),I=1,16)
2060	FORMAT(/' VNIX'/(1X,8E12.5))
	CALL B02(VNIX,INEL,VBE)
	IF (M.GE.2) WRITE(MP,2070)(VBE(I),I=1,48)
2070	FORMAT(/' MATRIX B'/(1X,10E12.5))

	CALL BTDB(VKE,VBE,VDE1,IDLE,IMATD,NSYM)
330	I1=I1+3*INEL

      RETURN

C
C----------EVALUATE ELEMENT TANGENT MATRIX
400   CONTINUE
      RETURN

C
C----------MASS MATRIX

C--------INITIALIZE VKE
500   DO 510 I=1,136
        VKE(I)=ZERO
510   CONTINUE

C------------LOOP OVER THE G.P
      IDIM1=NDIM-1
      IDECL=(NDIM+1)*INEL
      I1=1+INEL
      I2=0
      DO 550 IG=1,IPG
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
	D=VCPG(IG)*DETJ*VPREE(4)

C------------ACCUMULATE MASS TERMS
	IDL=0
	DO 540 J=1,INEL
	  JJ=I2+J
	  J0=1+IDL*(IDL+1)/2
	  DO 530 I=1,J
	    II=I2+I
	    C=VNI(II)*VNI(JJ)*D
	    VKE(J0)=VKE(J0)+C
	    IF (NDIM.EQ.1) GO TO 530
	    J1=J0+IDL+2
	    DO 520 II=1,IDIM1
	      VKE(J1)=VKE(J1)+C
520	      J1=J1+J1+1
530	    J0=J0+NDIM
540	  IDL=IDL+NDIM
	I1=I1+IDECL
550	I2=I2+IDECL



      IF (NSYM.NE.0) THEN
C..........Copy to "help-array"
        DO I=1,IDLE*(IDLE+1)/2
          VKE1(I)=VKE(I)
        ENDDO

C.........Copy upper triangle (in fact rearrange array)
        PNTCO=-IDLE+1
        PNTCO1=1
        DO COL=1,IDLE
          PNTCO1=PNTCO1+COL-1
          PNTCO=PNTCO+IDLE
          PNTRW1=0
          DO ROW=1,IDLE
            IF (COL.GE.ROW) THEN
              VKE(PNTCO+ROW-1)=VKE1(PNTCO1+ROW-1)
              PNTRW1=PNTRW1+1
            ELSE
              PNTRW1=PNTRW1+ROW-1
              VKE(PNTCO+ROW-1)=VKE1(PNTCO1+PNTRW1-1)         
            END IF
          END DO
        END DO
      ENDIF

      RETURN


C    
C-------------EVALUATE THE ELEMENT RESIDUAL
600   CALL D02(VPREE,VDE)

C--------------INITIALIZE THE RESIDUAL VECTOR
      DO 610 ID=1,IDLE
	VFE(ID)=ZERO
610   CONTINUE

C-------------LOOP OVER THE G.P.
      I1=1+INEL
      DO 640 IG=1,IPG
C---------EVALUATE THE JACOBIAN
	CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
C---------EVALUATE FUNCTIONS D(NI)/DX
	CALL DNIDX(VNI(I1),VJ1,NDIM,INEL,VNIX)

C---------EVALUATE STRAINS AND STRESSES
	EPSX=ZERO
	EPSY=ZERO
	GAMXY=ZERO
	ID=1
	DO 620 IN=1,INEL
	  UN=VDLE(ID)
	  VN=VDLE(ID+1)
	  C1=VNIX(IN)
	  IN1=IN+INEL
	  C2=VNIX(IN1)
	  EPSX=EPSX+C1*UN
	  EPSY=EPSY+C2*VN
	  GAMXY=GAMXY+C1*VN+C2*UN
620	  ID=ID+2

	C1=VCPG(IG)*DETJ
	C2=VDE(2)*C1
	C3=VDE(9)*C1
	C1=VDE(1)*C1
	SIGX=C1*EPSX+C2*EPSY
	SIGY=C2*EPSX+C1*EPSY
	TAUXY=C3*GAMXY

C-------------FORM THE RESIDUAL
	ID=1
	DO 630 IN=1,INEL
	  C1=VNIX(IN)
	  IN1=IN+INEL
	  C2=VNIX(IN1)
	  VFE(ID)=VFE(ID)+C1*SIGX+C2*TAUXY
	  VFE(ID+1)=VFE(ID+1)+C2*SIGY+C1*TAUXY
630	  ID=ID+2
640	I1=I1+3*INEL

      RETURN


C
C------------EVALUATE VOLUMIC FORCES, FX, FY PER UNIT VOLUME
C------------(FOR GRAVITY FX=0, FY=-VPREE(4)
C
700   FX=ZERO
      FY=-VPREE(4)
      DO 710 I=1,16
	VFE(I)=ZERO
710   CONTINUE

      I1=1
      IDECL=(NDIM+1)*INEL
      DO 730 IG=1,IPG
	CALL JACOB(VNI(I1+INEL),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
	DX=VCPG(IG)*DETJ
	DY=DX*FY
	DX=DX*FX
	I2=I1
	I3=1
	DO 720 IN=1,INEL
	  VFE(I3)=VFE(I3)+DX*VNI(I2)
	  VFE(I3+1)=VFE(I3+1)+DY*VNI(I2)
	  I2=I2+1
720	  I3=I3+2
730	I1=I1+IDECL

      RETURN

C
C-------------EVALUATE AND PRINT STRESSES AT G.P.
800   WRITE(MP,2080)IEL
2080  FORMAT(//' STRESSES IN ELEMENT:',I5//
     1 '  P.G.',7X,'X',11X,'Y',9X,'EPSX',8X,'EPSY',7X,'GAMXY',8X,'SIGX',
     2 8X,'SIGY',7X,'TAUXY',8X,'TETA'/ 71X ,'SIG1',8X,'SIG2',7X,'TAUMAX'
     3	/)

C----------FORM THE MATRIX D
      CALL D02(VPREE,VDE)
C----------LOOP OVER G.P.
      I1=1+INEL
      I2=0
      DO 820 IG=1,IPG
C---------EVALUATE THE JACOBIAN
	CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
C---------EVALUATE FUNCTIONS D(NI)/DX
	CALL DNIDX(VNI(I1),VJ1,NDIM,INEL,VNIX)

C---------COMPUTE STRAINS AND STRESSES AT G.P.
	EPSX=ZERO
	EPSY=ZERO
	GAMXY=ZERO
	X=ZERO
	Y=ZERO
	ID=1
	DO 810 IN=1,INEL
	  UN=VDLE(ID)
	  VN=VDLE(ID+1)
	  XN=VCORE(ID)
	  YN=VCORE(ID+1)
	  C1=VNIX(IN)
	  IN1=IN+INEL
	  C2=VNIX(IN1)
	  IN1=IN+I2
	  C3=VNI(IN1)
	  EPSX=EPSX+C1*UN
	  EPSY=EPSY+C2*VN
	  GAMXY=GAMXY+C1*VN+C2*UN
	  X=X+C3*XN
	  Y=Y+C3*YN
810	  ID=ID+2

C--------------COMPUTE THE STRESSES
	SIGX=VDE(1)*EPSX+VDE(2)*EPSY
	SIGY=VDE(2)*EPSX+VDE(1)*EPSY
	TAUXY=VDE(9)*GAMXY

C--------------COMPUTE THE PRINCIPAL STRESSES
	TETA=ATAN2(DEUX*TAUXY,SIGX-SIGY)*X05
	TETA=TETA*RADN
	C1=(SIGX+SIGY)*X05
	C2=(SIGX-SIGY)*X05
	TAUMAX=SQRT(C2*C2+TAUXY*TAUXY)
	SIG1=C1+TAUMAX
	SIG2=C1-TAUMAX
	WRITE(MP,2090)IG,X,Y,EPSX,EPSY,GAMXY,SIGX,SIGY,TAUXY,
     1		      TETA,SIG1,SIG2,TAUMAX
2090	FORMAT(1X,I5,8E12.5,5X,F5.1/66X,3E12.5)
	I2=I2+3*INEL
820	I1=I1+3*INEL


      RETURN
      END
