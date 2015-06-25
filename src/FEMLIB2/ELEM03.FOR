      SUBROUTINE ELEM03(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1      VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C********************************************************************C
C ELEM03.FOR                                                         C
C  3/1-1995     NTH                                                  C
C  Ulf R. Kristiansen & Tonni F. Johansen                            C
C********************************************************************C
C  Quadratic element for anisotropic harmonic problems in 1,2 or 3   C
C  dimensions, primarely for acoustic waves in fluid:                C
C     1 dimension: 3 nodes element                                   C
C     2 dimension: 8 nodes isoparametric element                     C
C     3 dimension: 20 nodes isoparametric element                    C
C  Number of integration points: 3 in each direction                 C
C  Number of degrees of freedom per node: 1                          C
C  Element matrix or vector formed by this subbprogram               C
C  according to icode value :                                        C
C     ICODE.EQ.1  return of parameters                               C
C     ICODE.EQ.2  evaluate interpolation functions and               C
C                 numerical integration coeffisients                 C
C     ICODE.EQ.3  element matrix (VKE)                               C
C     ICODE.EQ.4  tangent matrix (VKE)  .....not written.....NB!!    C
C     ICODE.EQ.5  mass matrix (VKE)                                  C
C     ICODE.EQ.6  K.U. product (VFE)                                 C
C     ICODE.EQ.7  element load (VFE)  .....not written.....NB!!      C
C     ICODE.EQ.8  print gradients                                    C
C     ICODE.EQ.9  impedance matrix(VKE)                              C
C     ICODE.EQ.10 source matrix (VSRE), normal velocity              C
C     ICODE.EQ.15 source matrix (VSRE), vector velocity              C
C     ICODE.EQ.20 compute velocities (VIE)                           C
C  Element properties                                                C
C     VPREE(1)  Density of fluid                                     C
C     VPREE(2)  Sound speed in fluid                                 C
C********************************************************************C
C  Based on Dhatt & Touzot, p.194-197                                C
C                                                                    C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
c $declare
      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VKE1(1),VFE(1)
      COMPLEX*16 VIE(1),VIE1(1),VNDSRC(1),VPRSG(1),VSRE(1),VDLEC(1)
      INTEGER KTYPRSG(1)


C-----------CHARACTERISTIC DIMENSIONS OF THE ELEMENT
C----------(VALID  UP TO 3 DIMENSIONS)
C     DIMENSION VCPG(IPG),VKPG(NDIM*IPG),XYZ(NDIM)
      REAL*8 VCPG(27),VKPG(81),XYZ(3)      
C     DIMENSION VJ(NDIM*NDIM),VJ1(NDIM*NDIM)
      REAL*8 VJ(9),VJ1(9)
      COMPLEX*16 VJC(9)
C     DIMENSION VNIX(INEL*NDIM),VNI((1+NDIM)*INEL*IPG),IPGKED(NDIM)
      REAL*8 VNIX(60),VNI(2160)
      INTEGER IPGKED(3)
      REAL*8 VKPG1(81),VCPG1(27)
      INTEGER IVK(9)                              

      SAVE VCPG,VKPG,VNI

      
      REAL*8 ZERO,EPS,PI,DETJ,COEF,X,Y,XXYY,X2,Y2,XS2,XY2,A,B,C,YS2
      INTEGER IKE,I,INI,IG,INI1,J,IK,I1,I2,IJ,IDECL,INI0,IC,I0,IN,IPG1
      INTEGER II,COL,ROW,PNTCO,PNTCO1,PNTRW1

      COMPLEX*16 CC,CZERO,V(8)
      REAL*8 ALMOSTZ,NR(3),KX,KY
C------------Number of Gauss points in ksi,eta,dzeta direction
      DATA IPGKED/3,3,3/
C---------------------------------------------------
      DATA ZERO/0.0D0/,EPS/1.D-6/,ALMOSTZ/1.D-16/,CZERO/(0.0,0.0)/
      DATA PI/3.141592654/
      DATA IVK/2,5,8,29,32,35,56,59,62/
      IKE=IDLE*(IDLE+1)/2
    
C-------------CHOOSE FUNCTION TO BE EXECUTED
      GO TO (100,200,300,400,500,600,700,800,900,1000,9999,9999,9999,
     1 9999,1500,9999,9999,9999,9999,20000),ICODE

C
C------------RETURN ELEMENT PARAMETERS IN COMMON /RGDT/
100   IF (NDIM.EQ.1) THEN
        IDLE0=3
        INEL0=3
        IPG0=3
      ELSEIF (NDIM.EQ.2) THEN
        IDLE0=8
        INEL0=8
        IPG0=9
      ELSEIF (NDIM.EQ.3) THEN
        IDLE0=20
        INEL0=20
        IPG0=27
      ENDIF 
      CFLAG=0
      
      RETURN

C
C-----------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C-----------THEIR DERIVATIVES AT G.P.

200   CALL GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
      CALL NI01(VKPG,VNI)
   
      RETURN

C----------COMPUTE ELEMENT STIFFNESS MATRIX
C

C--------INITIALIZE VKE
300   DO 310 I=1,IKE
        VKE(I)=ZERO
310   CONTINUE

C--------LOOP OVER THE INTEGRATION POINTS
      INI=1+INEL
      DO 330 IG=1,IPG
C----------EVALUATE THE JACOBIAN MATRIX, ITS INVERSE AND ITS DETER,
        CALL JACOB(VNI(INI),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF (DETJ.LT.EPS) WRITE(MP,2000)IEL,IG,DETJ
2000    FORMAT('*** ELEM ',I5,'P.G. ',I3,' DETJ(J)=',E12.5,
     1       ' ,     ABS(DETJ(J)) USED')
        DETJ=ABS(DETJ)

C-----------PERFORM DETJ*WEIGHT
        COEF=VCPG(IG)*DETJ

C-----------IF AXISYMMETRIC PROBLEM
        IF(IAX.EQ.1) THEN
           X=ZERO
           Y=ZERO
c...............initialiser for riktige plasser i vni matrisa
           INI1=(IG-1)*(NDIM+1)*INEL
           DO  J=1,INEL
           X=X+VCORE(2*J-1)*VNI(INI1+J)
           Y=Y+VCORE(2*J)*VNI(INI1+J)
           END DO
           XXYY=X*X+Y*Y
           X2=X/XXYY
           Y2=Y/XXYY
           XS2=X/SQRT(XXYY)
           YS2=Y/SQRT(XXYY)
           XY2=1./(XXYY)
        END IF

C-----------EVALUATE FUNCTIONS D(NI)/D(X)
        CALL DNIDX(VNI(INI),VJ1,NDIM,INEL,VNIX)

C-----------IF AXISYMMETRIC PROBLEM EVAL. D(NI)/D(R)
        IF(IAX.EQ.1) THEN
          DO  I=1,INEL
          A=VNIX(I)
          B=VNIX(INEL+I)
          VNIX(I)=A*XS2+B*YS2
          VNIX(INEL+I)=-A*Y+B*X
          END DO
          GO TO 319
        END IF

C------------ACCUMULATE TERMS OF THE ELEMENT MATRIX
        IK=0
        DO 320 J=1,IDLE
          DO 320 I=1,J
            I1=I
            I2=J
            C=ZERO
            DO 315 IJ=1,NDIM
              C=C+VNIX(I1)*VNIX(I2)
              I1=I1+IDLE
              I2=I2+IDLE
315         CONTINUE
            IK=IK+1
            VKE(IK)=VKE(IK)+C*COEF
320     CONTINUE

C-----------ACCUMULATE TERMS OF AXISYMMETRIC THE ELEMENT MATRIX
319     IF(IAX.EQ.1) THEN
           IK=0
           DO 321 J=1,IDLE
             DO 321 I=1,J
               I1=I
               I2=J
               C=VNIX(I1)*VNIX(I2)
               C=C+XY2*VNIX(I1+IDLE)*VNIX(I2+IDLE)
               IK=IK+1
               VKE(IK)=VKE(IK)+C*COEF*2.*Y*PI
321        CONTINUE
        END IF             

C----------NEXT G.P.
        INI=INI+(NDIM+1)*INEL
330   CONTINUE  

C........If nonsymmetrical problem, copy triangle to lower triangle
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
C----------EVALUATE ELEMENT TANGENT MATRIX
400   CONTINUE


      RETURN

C
C----------MASS MATRIX

C--------INITIALIZE VKE
500   DO 510 I=1,IKE
        VKE(I)=ZERO
510   CONTINUE

C      IF (VPREE(4).EQ.ZERO) RETURN

      INI=0
      DO 530 IG=1,IPG
C----------EVALUATE THE JACOBIAN MATRIX
        I1=INI+INEL+1
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        detj=abs(detj)

C-----------EVALUATE THE COEFFISIENT
        COEF=VCPG(IG)*DETJ
C.................if axi-symmetric problem
        IF(IAX.EQ.1) THEN
         Y=ZERO
         DO  J=1,INEL
          Y=Y+VCORE(2*J)*VNI(INI+J)
         END DO
        COEF=COEF*2.*PI*Y
        END IF

C------------TERMS OF THE MASS MATRIX
        IK=0
        DO 520 J=1,IDLE
          DO 520 I=1,J
            IK=IK+1
            I1=INI+I
            I2=INI+J
            VKE(IK)=VKE(IK)+VNI(I1)*VNI(I2)*COEF/(VPREE(2)*VPREE(2))
520     CONTINUE

C----------NEXT G.P.
        INI=INI+(NDIM+1)*INEL
530   CONTINUE

C........If nonsymmetrical problem, copy triangle to lower triangle
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
600   DO 605 I=1,INEL
        VFE(I)=ZERO
605   CONTINUE

      INI=1+INEL
      DO 640 IG=1,IPG
C---------EVALUATE THE JACOBIAN MATRIX AND THE DERIVATIVES OF N IN X,Y,Z
        CALL JACOB(VNI(INI),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        CALL DNIDX(VNI(INI),VJ1,NDIM,INEL,VNIX)

C---------COMPUTE THE COMMON COEFFISIENT
        COEF=VCPG(IG)*DETJ

C--------VPREE*B*VDLE PRODUCT
        I1=0
        DO 620 I=1,NDIM
          C=ZERO
          DO 610 J=1,INEL
            I1=I1+1
            C=C+VNIX(I1)*VDLE(J)
610       CONTINUE
          VJ(I)=C*COEF*VPREE(I)
620     CONTINUE

C---------(BT)*VJ PRODUCT
        DO 630 I=1,INEL
          I1=I-INEL
          DO 630 J=1,NDIM
            I1=I1+INEL
            VFE(I)=VFE(I)+VNIX(I1)*VJ(J)
630     CONTINUE
        INI=INI+(NDIM+1)*INEL
640   CONTINUE
      RETURN

C
C------------EVALUATE FE
700   CONTINUE
      RETURN

C
C-------------EVALUATE AND PRINT GRADIENTS AT G.P.
800   DO I=1,NDIM
        VIE(I)=0.0
      ENDDO

      WRITE(MP,2010)IEL
2010  FORMAT(//' VELOCITIES IN ELEMENT:',I4//)
      IDECL=(NDIM+1)*INEL
      INI0=1
      INI=1+INEL
      DO 830 IG=1,IPG
        CALL JACOB(VNI(INI),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        CALL DNIDX(VNI(INI),VJ1,NDIM,INEL,VNIX)

C------------EVALUATE THE COORDINATES OF THE G.P.
        DO 803 I=1,NDIM
          XYZ(I)=ZERO
803     CONTINUE

        IC=1
        I0=INI0
        DO 807 IN=1,INEL
          C=VNI(I0)
          DO 805 I=1,NDIM
            XYZ(I)=XYZ(I)+C*VCORE(IC)
            IC=IC+1
805       CONTINUE
          I0=I0+1
807     CONTINUE
        
C-------------EVALUATE THE GRADIENT
        I1=0
        DO 820 I=1,NDIM
          CC=CZERO
          DO 810 J=1,IDLE
            I1=I1+1
            CC=CC+VNIX(I1)*VDLEC(J)
810       CONTINUE
          VJC(I)=CC/CMPLX(0.,-VPREE(1)*W)
820     CONTINUE

C------------PRINT THE GRADIENT
        WRITE(MP,2020)IG,(XYZ(I),I=1,NDIM)
2020    FORMAT(5X,'G.P. :',I3,' COORDINATES :',3E12.5)
        WRITE(MP,2025)(VJC(I),I=1,NDIM)
2025    FORMAT(15X,'VELOCITIES   :',6E12.5)
        
C------------ACCUMULATE VELOCITY
        DO I=1,NDIM
          VIE(I)=VIE(I)+VJC(I)
        ENDDO


        INI0=INI0+IDECL
        INI=INI+IDECL
830   CONTINUE

C------------AVERAGE VELOCITY
        DO I=1,NDIM
          VIE(I)=VIE(I)/IPG
        ENDDO

      WRITE(MP,2030)
2030  FORMAT(//)
      RETURN
C
C ----------EVALUATE THE IMPEDANCE MATRIX
900   CONTINUE


cc190393cc      IF(IFL.NE.IFL1) THEN
C ............gauss punktene
 
       GO TO(911,912,913)NDIM
911     CONTINUE
        WRITE(MP,9110)
9110    FORMAT(5X,'NOT YET ANY ONE DIM ELEMENT')
        RETURN
912     CONTINUE
C.............still IPG temporart til 1dim integrasjon
       IPG1=IPG
       IPG=3
       VCPG1(1)=5./9.
       VCPG1(2)=8./9.
       VCPG1(3)=5./9.
          IF(IFL.GT.3) GO TO 915
          IF(IFL.EQ.2) THEN 
            VKPG1(1)=-1.
          ELSE IF(IFL.EQ.3) THEN
            VKPG1(1)=1. 
          ENDIF    
            VKPG1(3)=VKPG1(1)
            VKPG1(5)=VKPG1(1)
            VKPG1(2)=VKPG(2)
            VKPG1(4)=VKPG(4)
            VKPG1(6)=VKPG(6)    
          GO TO 916
915       IF(IFL.EQ.4) THEN
            VKPG1(2)=-1.
          ELSE IF(IFL.EQ.5) THEN
            VKPG1(2)=1.
          ENDIF
            VKPG1(4)=VKPG1(2)
            VKPG1(6)=VKPG1(2)
            VKPG1(1)=VKPG(1)
            VKPG1(3)=VKPG(7)
            VKPG1(5)=VKPG(13)
916       CONTINUE
        GO TO 917
913     CONTINUE
C...............three dimensional case
C...............define first VCPG and VKPG for on the surface 
C...............integration
      IPG1=IPG
      IPG=9
      VCPG1(1)=(5./9.)*(5./9.)
      VCPG1(2)=(5./9.)*(8./9.)
      VCPG1(3)=VCPG1(1)
      VCPG1(4)=VCPG1(2)
      VCPG1(5)=(8./9.)*(8./9.)
      VCPG1(6)=VCPG1(2)
      VCPG1(7)=VCPG1(1)
      VCPG1(8)=VCPG1(2)
      VCPG1(9)=VCPG1(1)
      IF(IFL.LT.4) THEN
      DO 9131 I=1,27
9131  VKPG1(I)=VKPG(I)
      IF(IFL.EQ.2) II=-1.
      IF(IFL.EQ.3) II=1.
      DO 9132 I=1,25,3
9132  VKPG1(I)=II
      ENDIF
      IF(IFL.GT.3.AND.IFL.LT.6) THEN
      DO 9133 I=1,9
9133  VKPG1(I)=VKPG(I)
      DO 9134 I=1,9
9134  VKPG1(9+I)=VKPG(27+I)
      DO 9135 I=1,9
9135  VKPG1(18+I)=VKPG(54+I)
      IF(IFL.EQ.4) II=-1.
      IF(IFL.EQ.5) II=1.
      DO 9136 I=1,9
C9136  VKPG1(IVK(I))=II
9136  VKPG1((I-1)*3+2)=II
      ENDIF
      IF(IFL.GE.6) THEN
      DO 9137 I=1,9
9137  VKPG1((I-1)*3+1)=VKPG((I-1)*9+1)
      DO 9138 I=1,9
9138  VKPG1((I-1)*3+2)=VKPG((I-1)*9+2)
      DO 9139 I=1,9
9139  VKPG1((I-1)*3+3)=VKPG((I-1)*9+3)
      IF(IFL.EQ.6) II=-1.
      IF(IFL.EQ.7) II=1.
      DO 9140 I=1,9
9140  VKPG1((I-1)*3+3)=II
      ENDIF

cc190393cc       ENDIF
917    CONTINUE
       IFL1=IFL

C..............funksjoner etc. for alle tilfelle av dim
       CALL NI01(VKPG1,VNI)
            DO 909 I=1,INEL*(NDIM+1)*IPG
909    IF(ABS(VNI(I)).LT.(1.D-10)) VNI(I)=ZERO

C.............nullstill matrisa             
       DO 910 I=1,IKE
910    VKE(I)=ZERO
       INI=0
        DO 930 IG=1,IPG
        I1=INI+INEL+1
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF(NDIM.EQ.3) GO TO 940
C...............two dimensional coeff
C.............. calculate the jacobian matrix VJ , determinant
C               will return as zero ,inverse not calculated
        IF(IFL.GT.3) THEN
          DETJ=SQRT(VJ(1)*VJ(1)+VJ(3)*VJ(3))
          ELSE
          DETJ=SQRT(VJ(2)*VJ(2)+VJ(4)*VJ(4))
        ENDIF
C............. compute the weight
        COEF=VCPG1(IG)*DETJ
C.................if axi-symmetric problem
        IF(IAX.EQ.1) THEN
         Y=ZERO
         DO  J=1,INEL
          Y=Y+VCORE(2*J)*VNI(INI+J)
         END DO
         COEF=COEF*2.*PI*Y
        END IF
      GO TO 950
C.............tre dim tilfelle coeff
940   continue
      IF(IFL.LT.4) THEN
      DETJ=(VJ(5)*VJ(9)-VJ(8)*VJ(6))**2
      DETJ=DETJ+(VJ(8)*VJ(3)-VJ(2)*VJ(9))**2
      DETJ=DETJ+(VJ(2)*VJ(6)-VJ(5)*VJ(3))**2
      COEF=SQRT(DETJ)*VCPG1(IG)
      ENDIF
      IF(IFL.GT.3.AND.IFL.LT.6) THEN
      DETJ=(VJ(4)*VJ(9)-VJ(7)*VJ(6))**2
      DETJ=DETJ+(VJ(7)*VJ(3)-VJ(1)*VJ(9))**2
      DETJ=DETJ+(VJ(1)*VJ(6)-VJ(4)*VJ(3))**2
      COEF=SQRT(DETJ)*VCPG1(IG)
      ENDIF
      IF(IFL.GT.5) THEN
      DETJ=(VJ(4)*VJ(8)-VJ(7)*VJ(5))**2
      DETJ=DETJ+(VJ(7)*VJ(2)-VJ(1)*VJ(8))**2
      DETJ=DETJ+(VJ(1)*VJ(5)-VJ(4)*VJ(2))**2
      COEF=SQRT(DETJ)*VCPG1(IG)
      ENDIF
C............. terms of the impedance matrix
950   CONTINUE
      IK=0    
      DO 920 J=1,IDLE
      DO 920 I=1,J
      IK=IK+1
      I1=INI+I
      I2=INI+J
920   VKE(IK)=VKE(IK)+VNI(I1)*VNI(I2)*COEF     
930   INI=INI+(NDIM+1)*INEL

C...............tilbakestill IGP
      IPG=IPG1
      

   
C........If nonsymmetrical problem, copy triangle to lower triangle
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


C-----------EVALUATE THE SOURCE MATRIX -- NORMAL VELOCITY
 
      RETURN

C-----------EVALUATE THE SOURCE VECTOR
1000  IF(IFL.NE.IFL1) THEN
C ............gauss points
 
        GO TO(1011,1012,1013)NDIM
1011    CONTINUE
1012    CONTINUE
C.............put IPG temporarely to 1-dim integration
        IPG1=IPG
        IPG=3
        VCPG1(1)=5./9.
        VCPG1(2)=8./9.
        VCPG1(3)=5./9.
        IF(IFL.LE.3) THEN
          IF(IFL.EQ.2) THEN 
            VKPG1(1)=-1.
          ELSE IF(IFL.EQ.3) THEN
            VKPG1(1)=1. 
          ENDIF    
          
          VKPG1(3)=VKPG1(1)
          VKPG1(5)=VKPG1(1)
          VKPG1(2)=VKPG(2)
          VKPG1(4)=VKPG(4)
          VKPG1(6)=VKPG(6)    
        ELSE
          IF(IFL.EQ.4) THEN
            VKPG1(2)=-1.
          ELSE IF(IFL.EQ.5) THEN
            VKPG1(2)=1.
          ENDIF
          
          VKPG1(4)=VKPG1(2)
          VKPG1(6)=VKPG1(2)
          VKPG1(1)=VKPG(1)
          VKPG1(3)=VKPG(7)
          VKPG1(5)=VKPG(13)
        ENDIF
1013    CONTINUE
      ENDIF
      IFL1=IFL

C..............functions etc.
      CALL NI01(VKPG1,VNI)
      DO 1009 I=1,INEL*(NDIM+1)*IPG
1009    IF(ABS(VNI(I)).LT.(1.D-10)) VNI(I)=ZERO

C.............initiate matrix
      DO 1010 I=1,IDLE
1010    VSRE(I)=CZERO

C-------------Find the velocities.      
      IF (NDIM.EQ.2) THEN
        V(1)=VNDSRC((IFREQ-1)*4+1)
        V(2)=VNDSRC((IFREQ-1)*4+2)
        V(3)=VNDSRC((IFREQ-1)*4+3)
      ENDIF
      
      INI=0
      DO 1030 IG=1,IPG
C.............. calculate the jacobian matrix VJ , determinant
C               will return as zero ,inverse not calculated
        I1=INI+INEL+1
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF(IFL.GT.3) THEN
          DETJ=SQRT(VJ(1)*VJ(1)+VJ(3)*VJ(3))
        ELSE
          DETJ=SQRT(VJ(2)*VJ(2)+VJ(4)*VJ(4))
        ENDIF
        
C---------------COMPUTE THE SOURCE TERM IN THE GAUSS POINT
        IF (NDIM.EQ.2) THEN
          IF (IFL.EQ.2) THEN
            CC=VNI(INI+7)*V(1)+VNI(INI+8)*V(2)+VNI(INI+1)*V(3)
          ELSEIF (IFL.EQ.3) THEN
            CC=VNI(INI+3)*V(1)+VNI(INI+4)*V(2)+VNI(INI+5)*V(3)
          ELSEIF (IFL.EQ.4) THEN
            CC=VNI(INI+1)*V(1)+VNI(INI+2)*V(2)+VNI(INI+3)*V(3)
          ELSEIF (IFL.EQ.5) THEN
            CC=VNI(INI+5)*V(1)+VNI(INI+6)*V(2)+VNI(INI+7)*V(3)
          ENDIF
        ENDIF
        
C............. compute the weight
        CC=VCPG1(IG)*DETJ*CC*VPREE(1)*CMPLX(0.,-W)
        
                
C.................if axi-symmetric problem
        IF(IAX.EQ.1) THEN
          Y=ZERO
          DO  J=1,INEL
            Y=Y+VCORE(2*J)*VNI(INI+J)
          END DO
          CC=CC*2.*PI*Y
        END IF
        
        
C............. find the terms of the source matrix
        DO J=1,inel
          I2=INI+J
          VSRE(J)=VSRE(J)+VNI(I2)*CC
        ENDDO
1030    INI=INI+(NDIM+1)*INEL
      
C................... reset IGP
      IPG=IPG1

      RETURN


C-----------EVALUATE THE SOURCE VECTOR - any velocity
1500  IF(IFL.NE.IFL1) THEN
C ............gauss points
 
        GO TO(1511,1512,1513)NDIM
1511    CONTINUE
1512    CONTINUE
C.............put IPG temporarely to 1-dim integration
        IPG1=IPG
        IPG=3
        VCPG1(1)=5./9.
        VCPG1(2)=8./9.
        VCPG1(3)=5./9.
        IF(IFL.LE.3) THEN
          IF(IFL.EQ.2) THEN 
            VKPG1(1)=-1.
          ELSE IF(IFL.EQ.3) THEN
            VKPG1(1)=1. 
          ENDIF    
          
          VKPG1(3)=VKPG1(1)
          VKPG1(5)=VKPG1(1)
          VKPG1(2)=VKPG(2)
          VKPG1(4)=VKPG(4)
          VKPG1(6)=VKPG(6)    
        ELSE
          IF(IFL.EQ.4) THEN
            VKPG1(2)=-1.
          ELSE IF(IFL.EQ.5) THEN
            VKPG1(2)=1.
          ENDIF
          
          VKPG1(4)=VKPG1(2)
          VKPG1(6)=VKPG1(2)
          VKPG1(1)=VKPG(1)
          VKPG1(3)=VKPG(7)
          VKPG1(5)=VKPG(13)
        ENDIF
1513    CONTINUE
      ENDIF
      IFL1=IFL

C..............functions etc.
      CALL NI01(VKPG1,VNI)
      DO I=1,INEL*(NDIM+1)*IPG
        IF(ABS(VNI(I)).LT.(1.D-10)) VNI(I)=ZERO
      ENDDO

C.............initiate matrix
      DO I=1,IDLE
        VSRE(I)=CZERO
      ENDDO

C-------------Find the velocities.      
      IF (NDIM.EQ.2) THEN
        V(1)=VNDSRC((IFREQ-1)*4+1)
        V(2)=VNDSRC((IFREQ-1)*4+2)
        V(3)=VNDSRC((IFREQ-1)*4+3)
        KX=REAL(VNDSRC((IFREQ-1)*4+4))
        KY=AIMAG(VNDSRC((IFREQ-1)*4+4))
      ENDIF
      
      INI=0
      DO IG=1,IPG
C.............. calculate the jacobian matrix VJ , determinant
C               will return as zero ,inverse not calculated
        I1=INI+INEL+1
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF(IFL.GT.3) THEN
          DETJ=SQRT(VJ(1)*VJ(1)+VJ(3)*VJ(3))
          IF(IFL.EQ.4) THEN
            NR(1)=VJ(3)/DETJ
            NR(2)=-VJ(1)/DETJ
          ELSE
            NR(1)=-VJ(3)/DETJ
            NR(2)=VJ(1)/DETJ
          ENDIF  
        ELSE
          DETJ=SQRT(VJ(2)*VJ(2)+VJ(4)*VJ(4))
          IF(IFL.EQ.2) THEN
            NR(1)=-VJ(4)/DETJ
            NR(2)=VJ(2)/DETJ
          ELSE
            NR(1)=VJ(4)/DETJ
            NR(2)=-VJ(2)/DETJ
          ENDIF  
        ENDIF
        
C---------------COMPUTE THE SOURCE TERM IN THE GAUSS POINT
        IF (NDIM.EQ.2) THEN
          IF (IFL.EQ.2) THEN
            CC=VNI(INI+7)*V(1)+VNI(INI+8)*V(2)+VNI(INI+1)*V(3)
          ELSEIF (IFL.EQ.3) THEN
            CC=VNI(INI+3)*V(1)+VNI(INI+4)*V(2)+VNI(INI+5)*V(3)
          ELSEIF (IFL.EQ.4) THEN
            CC=VNI(INI+1)*V(1)+VNI(INI+2)*V(2)+VNI(INI+3)*V(3)
          ELSEIF (IFL.EQ.5) THEN
            CC=VNI(INI+5)*V(1)+VNI(INI+6)*V(2)+VNI(INI+7)*V(3)
          ENDIF
        ENDIF
        
C............. compute the weight
        CC=VCPG1(IG)*DETJ*CC*VPREE(1)*CMPLX(0.,-W)
        
                
C.................if axi-symmetric problem
        IF(IAX.EQ.1) THEN
          Y=ZERO
          DO  J=1,INEL
            Y=Y+VCORE(2*J)*VNI(INI+J)
          END DO
          CC=CC*2.*PI*Y
        END IF
        
        
C............. find the terms of the source matrix
        DO J=1,inel
          I2=INI+J
          VSRE(J)=VSRE(J)+VNI(I2)*CC*(KX*NR(1)+KY*NR(2))
        ENDDO

        INI=INI+(NDIM+1)*INEL
      ENDDO
      
C................... reset IGP
      IPG=IPG1

      RETURN



C...........Calculate velocity in nodes
20000  CONTINUE

C--------INITIALIZE VKE
      DO I=1,IKE
        VIE(I)=ZERO
      ENDDO

C-----------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C-----------THEIR DERIVATIVES AT NODES

      VKPG1(1)=-1
      VKPG1(2)=-1
      VKPG1(3)=0
      VKPG1(4)=-1
      VKPG1(5)=1
      VKPG1(6)=-1
      VKPG1(7)=1
      VKPG1(8)=0
      VKPG1(9)=1
      VKPG1(10)=1
      VKPG1(11)=0
      VKPG1(12)=1
      VKPG1(13)=-1
      VKPG1(14)=1
      VKPG1(15)=-1
      VKPG1(16)=0
      
      IPG0=IPG
      IPG=8
      CALL NI01(VKPG1,VNI)
      IPG=IPG0
      
      INI=1+INEL
      DO IG=1,INEL
C----------EVALUATE THE JACOBIAN MATRIX, ITS INVERSE AND ITS DETER,
        CALL JACOB(VNI(INI),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF (DETJ.LT.EPS) WRITE(MP,20010)IEL,IG,DETJ
20010    FORMAT('*** ELEM ',I5,'P.G. ',I3,' DETJ(J)=',E12.5,
     1       ' ,     ABS(DETJ(J)) USED')
C        DETJ=ABS(DETJ)


        IF (DETJ.GT.ALMOSTZ) THEN
C...........Must check whether DETJ is zero or not (collapsing nodes)
C-----------EVALUATE FUNCTIONS D(NI)/D(X)
          CALL DNIDX(VNI(INI),VJ1,NDIM,INEL,VNIX)

C-----------Compute velocity in nodes, put results in VIE
C-----------1:reset VIE
          DO I=1,NDIM
            VIE(NDIM*(IG-1)+I)=0.0
          ENDDO        

C-----------2:calculate gradient in node IG
          DO I=1,INEL
            DO J=1,NDIM
              VIE(NDIM*(IG-1)+J)=VIE(NDIM*(IG-1)+J)+
     1               VNIX(I+8*(J-1))*VDLEC(I)
            ENDDO
          ENDDO

C-----------3:From the gradient, find the velocity on node IG
            DO J=1,NDIM
              VIE(NDIM*(IG-1)+J)=VIE(NDIM*(IG-1)+J)/CMPLX(0.,VPREE(1)*W)
            ENDDO
        
        ENDIF
C........Finished check on DETJ.GT.0.0
        
        INI=INI+(NDIM+1)*INEL

      ENDDO      

      RETURN

9999  RETURN
      END

