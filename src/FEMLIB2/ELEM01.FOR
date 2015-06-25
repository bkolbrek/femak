      SUBROUTINE ELEM01(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1      VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C********************************************************************C
C ELEM01.FOR                                                         C
C  19/1-1995                                                         C
C  Tonni F. Johansen, NTH                                            C
C********************************************************************C
C  Quadratic element for anisotropic harmonic problems in 1,2 or 3   C
C  dimensions:                                                       C
C     1 dimension: 3 nodes element                                   C
C     2 dimension: 8 nodes isoparametric element                     C
C     3 dimension: 20 nodes isoparametric element                    C
C  Number of integration points: 2 in each direction                 C
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
C     ICODE.EQ.30  volume of the element                             C
C     ICODE.EQ.31 average over the element                           C
C  Element properties                                                C
C     VPREE(1)  coeffiesients DX                                     C
C     VPREE(2)  coeffiesients DY                                     C
C     VPREE(3)  coeffiesients DZ                                     C
C     VPREE(4)  specific heat capacity                               C
C********************************************************************C
C  Programmed after Dhatt & Touzot, p.194-197                        C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
C  $DECLARE

      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VKE1(1),VFE(1)
      COMPLEX*16 VIE(1),VIE1(1),VNDSRC(1),VPRSG(1),VSRE(1),VDLEC(1)
      INTEGER KTYPRSG(1)

C-----------CHARACTERISTIC DIMENSIONS OF THE ELEMENT
C----------(VALID  UP TO 3 DIMENSIONS)
C     DIMENSION VCPG(IPG),VKPG(NDIM*IPG),XYZ(NDIM)
      REAL*8 VCPG(9),VKPG(27),XYZ(3)      
C     DIMENSION VJ(NDIM*NDIM),VJ1(NDIM*NDIM)
      REAL*8 VJ(9),VJ1(9)
C     DIMENSION VNIX(INEL*NDIM),VNI((1+NDIM)*INEL*IPG),IPGKED(NDIM)
      REAL*8 VNIX(60),VNI(2160)
      INTEGER IPGKED(3)
      
      SAVE VCPG,VKPG,VNI

      REAL*8 ZERO,EPS,DETJ,COEF,C
      INTEGER IKE,I,INI,IG,IK,J,I1,I2,IJ,IDECL,INI0,IC,I0,IN
C------------Number of Gauss points in ksi,eta,dzeta direction
      DATA IPGKED/3,3,3/
C---------------------------------------------------
      DATA ZERO/0.0D0/,EPS/1.D-6/

      IKE=IDLE*(IDLE+1)/2
      
C-------------CHOOSE FUNCTION TO BE EXECUTED
      GO TO (100,200,300,400,500,600,700,800,9999,9999,9999,9999,9999,
     1       9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,
     2       9999,9999,9999,9999,9999,900,1000),ICODE

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
2000    FORMAT('*** ELEM ',I5,'P.G. ',I3,' DETJ(J)=',E12.5)

C-----------PERFORM DETJ*WEIGHT
        COEF=VCPG(IG)*DETJ
C-----------EVALUATE FUNCTIONS D(NI)/D(X)
        CALL DNIDX(VNI(INI),VJ1,NDIM,INEL,VNIX)

C------------ACCUMULATE TERMS OF THE ELEMENT MATRIX
        IK=0
        DO 320 J=1,IDLE
          DO 320 I=1,J
            I1=I
            I2=J
            C=ZERO
            DO 315 IJ=1,NDIM
              C=C+VNIX(I1)*VNIX(I2)*VPREE(IJ)
              I1=I1+IDLE
              I2=I2+IDLE
315         CONTINUE
            IK=IK+1
            VKE(IK)=VKE(IK)+C*COEF
320     CONTINUE

C----------NEXT G.P.
        INI=INI+(NDIM+1)*INEL
330   CONTINUE  
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

      IF (VPREE(4).EQ.ZERO) RETURN

      INI=0
      DO 530 IG=1,IPG
C----------EVALUATE THE JACOBIAN MATRIX
        I1=INI+INEL+1
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)

C-----------EVALUATE THE COEFFISIENT
        COEF=VCPG(IG)*DETJ*VPREE(4)

C------------TERMS OF THE MASS MATRIX
        IK=0
        DO 520 J=1,IDLE
          DO 520 I=1,J
            IK=IK+1
            I1=INI+I
            I2=INI+J
            VKE(IK)=VKE(IK)+VNI(I1)*VNI(I2)*COEF
520     CONTINUE

C----------NEXT G.P.
        INI=INI+(NDIM+1)*INEL
530   CONTINUE  
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
800   WRITE(MP,2010)IEL
2010  FORMAT(//' GRADIENTS IN ELEMENT:',I4//)
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
          C=ZERO
          DO 810 J=1,IDLE
            I1=I1+1
            C=C+VNIX(I1)*VDLE(J)
810       CONTINUE
          VJ(I)=C*VPREE(I)
820     CONTINUE

C------------PRINT THE GRADIENT
        WRITE(MP,2020)IG,(XYZ(I),I=1,NDIM)
2020    FORMAT(5X,'G.P. :',I3,' COORDINATES :',3E12.5)
        WRITE(MP,2025)(VJ(I),I=1,NDIM)
2025    FORMAT(15X,'GRADIENTS   :',3E12.5)
        
        INI0=INI0+IDECL
        INI=INI+IDECL
830   CONTINUE

      WRITE(MP,2030)
2030  FORMAT(//)

      RETURN

C----------COMPUTE THE AREA OF THE ELEMENT
C

C--------INITIALIZE VKE
900   VKE(1)=ZERO

C--------LOOP OVER THE INTEGRATION POINTS
      INI=1+INEL
      DO 930 IG=1,IPG
C----------EVALUATE THE JACOBIAN MATRIX, ITS INVERSE AND ITS DETER,
        CALL JACOB(VNI(INI),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF (DETJ.LT.EPS) WRITE(MP,2000)IEL,IG,DETJ

C-----------PERFORM DETJ*WEIGHT
        COEF=VCPG(IG)*DETJ
        VKE(1)=VKE(1)+COEF

C----------NEXT G.P.
        INI=INI+(NDIM+1)*INEL
930   CONTINUE  
      RETURN

C----------COMPUTE THE AVERAGE OVER THE ELEMENT
C

C--------INITIALIZE VKE
1000   VKE(1)=ZERO

C--------LOOP OVER THE INTEGRATION POINTS
      INI=0
      DO 1030 IG=1,IPG
C----------EVALUATE THE JACOBIAN MATRIX, ITS INVERSE AND ITS DETER,
        I1=INI+INEL+1
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF (DETJ.LT.EPS) WRITE(MP,2000)IEL,IG,DETJ

C-----------PERFORM DETJ*WEIGHT
        COEF=VCPG(IG)*DETJ

C------------ACCUMULATE TERMS
        DO 1020 J=1,IDLE
            VKE(1)=VKE(1)+VNI(INI+J)*VDLE(J)*COEF
1020     CONTINUE

C----------NEXT G.P.
        INI=INI+(NDIM+1)*INEL
1030   CONTINUE  

      RETURN

9999  CONTINUE
      RETURN

      END
