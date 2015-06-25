      SUBROUTINE ELEM11(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1      VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C********************************************************************C
C ELEM11.FOR                                                         C
C  3/1-1995     NTH                                                  C
C   Ulf R. Kristiansen & Tonni F. Johansen                           C
C********************************************************************C
C  Quadratic wave envelope element for  1,2 or 3                     C
C  dimensions:                                                       C
C     1 dimension: 3 nodes element                                   C
C     2 dimension: 8 nodes element  only this is implemented         C
C     3 dimension: 20 nodes element                                  C
C  The element is described by Astley in International Journal for   C
C  numerical methods in fluid dynamics vol3 1983                     C
C  The element is copied from elem04, icode=3 and 5 only valid       C
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
C     ICODE.EQ.10 source matrix (VKE)                                C
C  Element properties                                                C
C     VPREE(1)  Density of fluid                                     C
C     VPREE(2)  Sound speed in fluid                                 C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'

C $DECLARE

      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VKE1(1),VFE(1)
      COMPLEX*16 VIE(1),VIE1(1),VNDSRC(1),VPRSG(1),VSRE(1),VDLEC(1)
      INTEGER KTYPRSG(1)


      COMPLEX*16 H0,H0CC,H0I,H0ICC,H0J,H0JCC,H0D,H0DCC,BXI,BXJ,BYI,BYJ

C-----------CHARACTERISTIC DIMENSIONS OF THE ELEMENT
C----------(VALID  UP TO 3 DIMENSIONS)
C     DIMENSION VCPG(IPG),VKPG(NDIM*IPG),XYZ(NDIM)
      REAL*8 VCPG(27),VKPG(81),XYZ(3)
C     DIMENSION VJ(NDIM*NDIM),VJ1(NDIM*NDIM)
      REAL*8 VJ(9),VJ1(9)
C     DIMENSION VNIX(INEL*NDIM),VNI((1+NDIM)*INEL*IPG),IPGKED(NDIM)
      REAL*8 VNIX(60),VNI(2160)
      INTEGER IPGKED(3)
      REAL*8 VKPG1(81),VCPG1(27)
      INTEGER IVK(9)

      SAVE VCPG,VKPG,VNI


      INTEGER I,IKE,INI,INII,IG,INI1,J,IK,I1,J1,I2,INI0,IC,I0,IN
      INTEGER IPG1,II,IDECL
      REAL*8 A,RI,RJ,R,ZERO,EPS,PI,DETJ,COEF,X,Y,DRDX,DRDY,XXYY,XS2,YS2
      REAL*8 B,ARG,FTGN,C



C------------Number of Gauss points in ksi,eta,dzeta direction
      DATA IPGKED/3,3,3/
C---------------------------------------------------
      DATA ZERO/0.0D0/,EPS/1.D-6/
      DATA PI/3.141592654/
      DATA IVK/2,5,8,29,32,35,56,59,62/
      IKE=IDLE*IDLE

C-------------CHOOSE FUNCTION TO BE EXECUTED
      GO TO (100,200,300,400,500,600,700,800,900,1000),ICODE

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

      CFLAG=1
      RETURN

C
C-----------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C-----------THEIR DERIVATIVES AT G.P.

200        CALL GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
           CALL NI01(VKPG,VNI)

      RETURN

C----------COMPUTE ELEMENT STIFFNESS MATRIX
C

C--------INITIALIZE VIE
300   DO 310 I=1,IKE
        VIE(I)=DCMPLX(ZERO,ZERO)
310   CONTINUE


c           startparam for (det,dnidx) og formfunk
      INI=1+INEL
      INII=0

C--------LOOP OVER THE INTEGRATION POINTS

      DO 330 IG=1,IPG

C----------EVALUATE THE JACOBIAN MATRIX,USING TRAD SHAPE FUNCTIONS

        CALL JACOB(VNI(INI),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF (DETJ.LT.EPS) WRITE(MP,2000)IEL,IG,DETJ
2000    FORMAT('*** ELEM ',I5,'P.G. ',I3,' DETJ(J)=',E12.5,
     1       ' ,     ABS(DETJ(J)) USED')
        DETJ=ABS(DETJ)

C-----------PERFORM DETJ*WEIGHT
        COEF=VCPG(IG)*DETJ

        CALL DNIDX(VNI(INI),VJ1,NDIM,INEL,VNIX)

        IF(IAX.EQ.1) GO TO 319


C............evaluate x,y,r at gausspoint
         X=ZERO
         Y=ZERO
         INI1=(IG-1)*(NDIM+1)*INEL
         DO J=1,INEL
         X=X+VCORE(2*J-1)*VNI(INI1+J)
         Y=Y+VCORE(2*J)*VNI(INI1+J)
         ENDDO
         R=DSQRT(X*X+Y*Y)
         DRDX=X/DSQRT(X*X+Y*Y)
         DRDY=Y/DSQRT(X*X+Y*Y)
C
C------------ACCUMULATE TERMS OF THE ELEMENT MATRIX
        IK=0
        A=W/VPREE(2)
        DO 320 J=1,IDLE
          DO 320 I=1,IDLE
          IK=IK+1
          I1=INII+I
          J1=INII+J
          RI=DSQRT(VCORE(2*I-1)**2+VCORE(2*I)**2)
          RJ=DSQRT(VCORE(2*J-1)**2+VCORE(2*J)**2)
          H0=DCMPLX(BESJ0(R*A),BESY0(R*A))
          H0D=-A*DCMPLX(BESJ1(A*R),BESY1(A*R))
          H0I=DCMPLX(BESJ0(RI*A),BESY0(RI*A))
          H0J=DCMPLX(BESJ0(RJ*A),BESY0(RJ*A))
          H0CC=DCONJG(H0)
          H0DCC=DCONJG(H0D)
          H0ICC=DCONJG(H0I)
          H0JCC=DCONJG(H0J)

C.............X part
          BXI=(1./(H0I))*(VNIX(I)*H0+VNI(I1)*H0D*DRDX)
          BXJ=(1./(H0JCC))*(VNIX(J)*H0CC+VNI(J1)*H0DCC*DRDX)

C.............Y part
          BYI=(1./(H0I))*(VNIX(I+INEL)*H0+VNI(I1)*H0D*DRDY)
          BYJ=(1./(H0JCC))*(VNIX(J+INEL)*H0CC+VNI(J1)*H0DCC*DRDY)



          VIE(IK)=VIE(IK)+COEF*(BXI*BXJ+BYI*BYJ)

320     CONTINUE
        GO TO 322


C................... axisymmetric system
319     CONTINUE
           X=ZERO
           Y=ZERO
c...............initialiser for riktige plasser i vni matrisa
           INI1=(IG-1)*(NDIM+1)*INEL
           DO  J=1,INEL
           X=X+VCORE(2*J-1)*VNI(INI1+J)
           Y=Y+VCORE(2*J)*VNI(INI1+J)
           END DO
           XXYY=X*X+Y*Y
           R=DSQRT(XXYY)
           XS2=X/SQRT(XXYY)
           YS2=Y/SQRT(XXYY)



C-----------EVAL. D(NI)/D(R) og D(NI)/D(THETA)
          DO  I=1,INEL
          A=VNIX(I)
          B=VNIX(INEL+I)
          VNIX(I)=A*XS2+B*YS2
          VNIX(INEL+I)=-A*Y+B*X
          END DO
          IK=0
           DO 321 J=1,IDLE
             DO 321 I=1,IDLE
               IK=IK+1
               I1=INII+I
               J1=INII+J
               RI=DSQRT(VCORE(2*I-1)**2+VCORE(2*I)**2)
               RJ=DSQRT(VCORE(2*J-1)**2+VCORE(2*J)**2)
               ARG=W/VPREE(2)
               FTGN=-1.
c... test
               ftgn=1.

               H0=DCMPLX(DCOS(ARG*R)/R,FTGN*DSIN(ARG*R)/R)
               H0I=DCMPLX(DCOS(ARG*RI)/RI,FTGN*DSIN(ARG*RI)/RI)
               H0J=DCMPLX(DCOS(ARG*RJ)/RJ,FTGN*DSIN(ARG*RJ)/RJ)
               H0CC=DCONJG(H0)
               H0ICC=DCONJG(H0I)
               H0JCC=DCONJG(H0J)

C...............r part
               BXI=DCMPLX(1./R,-FTGN*ARG)
               BXI=(H0/H0I)*(VNIX(I)-VNI(I1)*BXI)
               BXJ=DCMPLX(1./R,FTGN*ARG)
               BXJ=(H0CC/H0JCC)*(VNIX(J)-VNI(J1)*BXJ)
C...............theta part
c... tatt bort 19/9 som minitest
c               BYI=(H0/H0I)*(VNIX(I1+IDLE))
c               BYJ=(H0CC/H0JCC)*(VNIX(J1+IDLE))
c................enda et forsoek20/9
               BYI=(H0/H0I)*(VNIX(I+IDLE))
               BYJ=(H0CC/H0JCC)*(VNIX(J+IDLE))


c


               VIE(IK)=VIE(IK)+COEF*2.*Y*PI*(BXI*BXJ+(1./(R*R))*BYI*BYJ)
321        CONTINUE




322     CONTINUE
C----------NEXT G.P.
        INI=INI+(NDIM+1)*INEL
        INII=INII+(NDIM+1)*INEL
330   CONTINUE

c        do ik=1,64
c         print*,' stivhhetsmatrise',ik,vie(ik)
c        enddo

      RETURN

C
C----------EVALUATE ELEMENT TANGENT MATRIX
400   CONTINUE


      RETURN

C
C----------MASS MATRIX

C--------INITIALIZE VKE
500   DO 510 I=1,IKE
        VIE(I)=CMPLX(ZERO,ZERO)
510   CONTINUE

      A=W/VPREE(2)


      INI=0
      DO 530 IG=1,IPG
C----------EVALUATE THE JACOBIAN MATRIX USING TRAD SHAPE FUNCTIONS
        I1=INI+INEL+1
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        detj=abs(detj)

C-----------EVALUATE THE COEFFISIENT
        COEF=VCPG(IG)*DETJ

c.............if axisymmetric problem
         IF(IAX.EQ.1) GO TO 519

C............evaluate x,y,r at gausspoint
         X=ZERO
         Y=ZERO
         INI1=(IG-1)*(NDIM+1)*INEL
         DO J=1,INEL
         X=X+VCORE(2*J-1)*VNI(INI1+J)
         Y=Y+VCORE(2*J)*VNI(INI1+J)
         ENDDO
         R=DSQRT(X*X+Y*Y)


C------------TERMS OF THE MASS MATRIX
        IK=0
        DO 520 J=1,IDLE
          DO 520 I=1,IDLE
            IK=IK+1
            I1=INI+I
            I2=INI+J
            RI=DSQRT(VCORE(2*I-1)**2+VCORE(2*I)**2)
            RJ=DSQRT(VCORE(2*J-1)**2+VCORE(2*J)**2)

            H0=DCMPLX(BESJ0(R*A),BESY0(R*A))
            H0I=DCMPLX(BESJ0(RI*A),BESY0(RI*A))
            H0J=DCMPLX(BESJ0(RJ*A),BESY0(RJ*A))

            H0CC=DCONJG(H0)
            H0ICC=DCONJG(H0I)
            H0JCC=DCONJG(H0J)
            BXI=VNI(I1)*H0/H0I
            BXJ=VNI(I2)*H0CC/H0JCC
c  ...den neste linja var fallt ut paa et eller annet tidspunkt
            VIE(IK)=VIE(IK)+COEF*BXI*BXJ*(-W*W)/(VPREE(2)*VPREE(2))

C


520     CONTINUE
        GO TO 522
C.............axisymmetric case
519     CONTINUE
C............evaluate x,y,r at gausspoint
         X=ZERO
         Y=ZERO
         INI1=(IG-1)*(NDIM+1)*INEL
         DO J=1,INEL
         X=X+VCORE(2*J-1)*VNI(INI1+J)
         Y=Y+VCORE(2*J)*VNI(INI1+J)
         ENDDO
         R=DSQRT(X*X+Y*Y)


C------------TERMS OF THE MASS MATRIX
        IK=0
        DO 521 J=1,IDLE
          DO 521 I=1,IDLE
            IK=IK+1
            I1=INI+I
            I2=INI+J
            RI=DSQRT(VCORE(2*I-1)**2+VCORE(2*I)**2)
            RJ=DSQRT(VCORE(2*J-1)**2+VCORE(2*J)**2)

            FTGN=-1.
c  ... test16/9
            ftgn=1.
            H0=DCMPLX(DCOS(A*R)/R,FTGN*DSIN(A*R)/R)
            H0I=DCMPLX(DCOS(A*RI)/RI,FTGN*DSIN(A*RI)/RI)
            H0J=DCMPLX(DCOS(A*RJ)/RJ,FTGN*DSIN(A*RJ)/RJ)
            H0CC=DCONJG(H0)
            H0ICC=DCONJG(H0I)
            H0JCC=DCONJG(H0J)
            BXI=VNI(I1)*H0/H0I
            BXJ=VNI(I2)*H0CC/H0JCC
C...... en test16/9
c            BXJ=VNI(I1)*H0/H0I
c            BXI=VNI(I2)*H0CC/H0JCC

            VIE(IK)=VIE(IK)+2.*PI*Y*COEF*BXI*BXJ
     +       *(-W*W)/(VPREE(2)*VPREE(2))



521     CONTINUE

522   CONTINUE
C----------NEXT G.P.
        INI=INI+(NDIM+1)*INEL
530   CONTINUE
c        do ik=1,64
c         print*,' massematrise',ik,vie(ik)
c        enddo

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
c    for debugger og spesialversjon
c       IFL=3
c       iflx=ifl
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
9136  VKPG1(IVK(I))=II
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
910    VIE(I)=DCMPLX(0.0,0.0)
       INI=0
       DO 930 IG=1,IPG
      I1=INI+INEL+1
      CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
C............evaluate x,y,r at gausspoint
c ...........lagt inn nov17 93
         X=ZERO
         Y=ZERO
         INI1=(IG-1)*(NDIM+1)*INEL
         DO J=1,INEL
         X=X+VCORE(2*J-1)*VNI(INI1+J)
         Y=Y+VCORE(2*J)*VNI(INI1+J)
         ENDDO
         R=DSQRT(X*X+Y*Y)

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
      DO 920 I=1,IDLE
      IK=IK+1
      I1=INI+I
      I2=INI+J
c ..............lagt inn nov17 93
            if(ifl.eq.4) then
            RI=DSQRT(VCORE(2*I-1)**2+VCORE(2*I)**2)
            RJ=DSQRT(VCORE(2*J-1)**2+VCORE(2*J)**2)

            H0=DCMPLX(BESJ0(R*A),BESY0(R*A))
            H0I=DCMPLX(BESJ0(RI*A),BESY0(RI*A))
            H0J=DCMPLX(BESJ0(RJ*A),BESY0(RJ*A))
c    miniscule verdier for test
c            h0=dcmplx(1.,0.)
c            h0i=dcmplx(1.,0.)
c            h0j=dcmplx(1.,0.)

            H0CC=DCONJG(H0)
            H0ICC=DCONJG(H0I)
            H0JCC=DCONJG(H0J)

            BXI=VNI(I1)*H0CC/H0ICC
            BXJ=VNI(I2)*H0/H0J
C......EN TEST
            BXI=VNI(I1)*H0/H0I
            BXJ=VNI(I2)*H0CC/H0JCC
            VIE(IK)=VIE(IK)+COEF*BXI*BXJ
            else

      VIE(IK)=VIE(IK)+VNI(I1)*VNI(I2)*COEF
            endif
920   continue

930   INI=INI+(NDIM+1)*INEL
c   midlertidig for utskrift
c      do  i=1,ike
c      print*,'ik,vIe(ik)',i,vie(i)
c      end do

C...............tilbakestill IGP
      IPG=IPG1

      RETURN


C-----------EVALUATE THE SOURCE MATRIX
1000  CONTINUE

      IF(IFL.NE.IFL1) THEN
C ............gauss punktene

       GO TO(1011,1012,1013)NDIM
1011     CONTINUE
1012     CONTINUE
C.............still IPG temporart til 1dim integrasjon
       IPG1=IPG
       IPG=3
       VCPG1(1)=8./9.
       VCPG1(2)=5./9.
       VCPG1(3)=8./9.
          IF(IFL.GT.3) GO TO 1015
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
          GO TO 1016
1015       IF(IFL.EQ.4) THEN
            VKPG1(2)=-1.
          ELSE IF(IFL.EQ.5) THEN
            VKPG1(2)=1.
          ENDIF
            VKPG1(4)=VKPG1(2)
            VKPG1(6)=VKPG1(2)
            VKPG1(1)=VKPG(1)
            VKPG1(3)=VKPG(7)
            VKPG1(5)=VKPG(13)
1016       CONTINUE
1013     CONTINUE
       ENDIF
       IFL1=IFL

C..............funksjoner etc.
       CALL NI01(VKPG1,VNI)
       DO 1009 I=1,INEL*(NDIM+1)*IPG
1009    IF(ABS(VNI(I)).LT.(1.D-10)) VNI(I)=ZERO

C.............nullstill matrisa
       DO 1010 I=1,IKE
1010    VFE(I)=ZERO

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

C............. terms of the source matrix
      DO 1020 J=1,IDLE
      I2=INI+J
1020  VFE(J)=VFE(J)+VNI(I2)*COEF
1030  INI=INI+(NDIM+1)*INEL

C................... tilbakestill IGP
      IPG=IPG1

      RETURN
      END
c************************************************
      function  besj0(x)
      implicit real*8(a-h,o-z)
      real*8 y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6
     +     ,s1,s2,s3,s4,s5,s6,AX,Z,XX
      data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     +  -.2073370639d-5,.2093887211d-6/,q1,q2,q3,q4,q5/-.1562499995d-1,
     +  .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     + 651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,
     +  s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,
     +  9494680.718d0,59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.) then
      y=x**2
      besj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     +     /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
      ax=abs(x)
      z=8./ax
      y=z**2
      xx=ax-.785398164
      besj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     +    *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      end
c***********************************************************************
      function besy0(x)
      implicit real*8(a-h,o-z)

      real*8 y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6
     +     ,s1,s2,s3,s4,s5,s6,Z,XX
      data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     +  -.2073370639d-5,.2093887211d-6/,q1,q2,q3,q4,q5/-.1562499995d-1,
     +  .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      data r1,r2,r3,r4,r5,r6/-2957821389.d0,7062834065.d0,-512359803.6d0
     +     ,10879881.29d0,-86327.92757d0,228.4622733d0/ ,
     +     s1,s2,s3,s4,s5,s6/40076544269.d0,745249964.8d0,
     +     7189466.438d0,47447.26470d0,226.1030244d0,1.d0/
      if(x.lt.8.) then
      y=x**2
      besy0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     +/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))+.636619772*besj0(x)*log(x)
      else
      z=8./x
      y=z**2
      xx=x-.785398164
      besy0=sqrt(.636619772/x)*(sin(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     +    *p5))))+z*cos(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      end

c***********************************************************************
      function besj1(x)
      implicit real*8(a-h,o-z)
      REAL*8 X
      REAL*8 AX,Z,XX
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6
      REAL*8 S1,S2,S3,S4,S5,S6

      DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,242396853.1D0
     1	  ,-2972611.439D0,15704.48260D0,-30.16036606D0/
      DATA S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0,
     1	   18583304.74D0,99447.43394D0,376.9991397D0,1.D0/
      DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,.2457520174D-5
     1	  ,-.240337019D-6/,Q1,Q2,Q3,Q4,Q5/.04687499995D0,
     2	   -.2002690873D-3,.8449199096D-5,-.88228987D-6,.105787412D-6/



      IF (ABS(X).LT.8) THEN
	Y=X*X
	BESJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))/
     1	       (S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
	AX=ABS(X)
	Z=8./AX
	Y=Z*Z
	XX=AX-2.356194491
	BESJ1=SQRT(.636619772/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*
     1		    P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
     2		    *SIGN(1.0D0,X)
      ENDIF

      RETURN
      END
c************************************************************************
      function besy1(x)
      implicit real*8(a-h,o-z)
      REAL*8 X
      REAL*8 Z,XX
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6
      REAL*8 S1,S2,S3,S4,S5,S6,S7
      DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,.2457520174D-5
     1	  ,-.240337019D-6/,Q1,Q2,Q3,Q4,Q5/.04687499995D0,
     2	   -.2002690873D-3,.8449199096D-5,-.88228987D-6,.105787412D-6/
      DATA R1,R2,R3,R4,R5,R6/-.4900604943D13,.1275274390D13,
     1    -.5153438139D11,.7349264551D9,-.4237922726D7,.8511937935D4/,
     2     S1,S2,S3,S4,S5,S6,S7/.2499580570D14,.4244419664D12,
     3  .3733650367D10,.2245904002D8,.1020426050D6,.3549632885D3,1.D0/



      IF (ABS(X).LT.8) THEN
	Y=X*X
	BESY1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))/
     1	     (S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*(S6+Y*S7))))))+.636619772
     2     *(BESJ1(X)*LOG(X)-1./X)
      ELSE

	Z=8./X
	Y=Z*Z
	XX=X-2.356194491
	BESY1=SQRT(.636619772/X)*(SIN(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*
     1		    P5))))+Z*COS(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))

      ENDIF

      RETURN
      END


