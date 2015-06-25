      SUBROUTINE ELEM20(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1      VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C********************************************************************C
C ELEM20.FOR                                                         C
C  30/5-1996                  N.T.H                                  C
C  Tonni F. Johansen                                                 C
C********************************************************************C
C  8 nodes quadratic element for 2 dimensional and axisymmetrical    C
C  elasticity, using "Rayleigh-loss"                                 C
C********************************************************************C
C  20 nodes quadratic element for 3 dimensional                      C
C  elasticity, using "Rayleigh-loss"                                 C
C********************************************************************C
C  Evaluate element information according to ICODE value             C
C     ICODE.EQ.1  return of parameters                               C
C     ICODE.EQ.2  evaluate interpolation functions and               C
C                 numerical integration coeffisients                 C
C     ICODE.EQ.3  stiffnes matrix (VKE)                              C
C     ICODE.EQ.4  tangent matrix (VKE)  .....not written.....NB!!    C
C     ICODE.EQ.5  mass matrix (VKE)                                  C
C     ICODE.EQ.6  residuals                                          C
C     ICODE.EQ.7  element load (VFE) (second member)                 C
C     ICODE.EQ.8  evaluate and print stresses                        C
C  Element properties                                                C
C     VPREE(1) Young's modulus                                       C
C     VPREE(2) Poisson's coeffiesient                                C
C     VPREE(3) Specific mass                                        C
C     VPREE(4) Rayleigh coeffisients alpha                           C
C     VPREE(5) Rayleigh coeffisients beta                            C
C********************************************************************C
C  Similar to ELEM02, programmert after Dhatt & Touzot, p.194-197    C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
c  $DECLARE

      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VKE1(1),VFE(1)
      COMPLEX*16 VIE(1),VIE1(1),VNDSRC(1),VPRSG(1),VSRE(1),VDLEC(1)
      INTEGER KTYPRSG(1)


C-----------CHARACTERISTIC DIMENSIONS OF THE ELEMENT
C----------(VALID  UP TO 3 DIMENSIONS)
C     DIMENSION VCPG(IPG),VKPG(NDIM*IPG)
      REAL*8 VCPG(27),VKPG(81)
C     DIMENSION VBE(IMATD*IDLE),VDE(IMATD**2),VJ(NDIM*NDIM),VJ1(NDIM*NDIM)
      REAL*8 VBE(360),VDE(36),VJ(9),VJ1(9),VDE1(36)
C     DIMENSION VNIX(INEL*NDIM),VNI((1+NDIM)*INEL*IPG),IPGKED(NDIM)
      REAL*8 VNIX(60),VNI(2160)
      REAL*8 VKPG1(81),VCPG1(27)
      INTEGER IPGKED(3)


      COMPLEX*16 V(8),CZERO,D
      REAL*8 PI,R,ZERO,DEUX,X05,RADN,EPS,X,Y,DETJ,C,EPSX,EPSY,GAMXY
      REAL*8 UN,VN,C1,C2,C3,SIGX,SIGY,TAUXY,FX,FY,DX,DY,XN,YN,TETA
      REAL*8 TAUMAX,SIG1,SIG2
      INTEGER IMATD,IKE,I0,IG,I1,I,IDECL,IDIM1,I2,IDL,J,JJ,J0,II,J1,ID
      INTEGER IN,IN1,I3
      INTEGER COL,ROW,PNTCO,PNTCO1,PNTRW1

C------------Dimension of D-matrix, Number of Gauss points in ksi,eta direct
C      DATA IMATD/3/,IPGKED/3,3/
      DATA IPGKED/3,3,3/
C---------------------------------------------------
      DATA ZERO/0.0D0/,DEUX/2.D0/,X05/0.5D0/,RADN/.572957795130823D2/
      DATA EPS/1.D-6/,CZERO/(0.0,0.0)/
      DATA PI/3.141592654/

c      SQRT(X)=DSQRT(X)
c      ATAN2(X,Y)=DATAN2(X,Y)

      IF (NSYM.EQ.0) THEN
        IKE=IDLE*(IDLE+1)/2
      ELSE
        IKE=IDLE*IDLE
      ENDIF

C--------------CHOOSE FUNCTION TO BE EXECUTED
      GO TO (100,200,300,400,500,600,700,800,900,1000),ICODE

C
C------------RETURN ELEMENT PARAMETERS IN COMMON /RGDT/
100   IF (NDIM.EQ.2) THEN
        IDLE0=16
        INEL0=8
        IPG0=9
      ELSEIF (NDIM.EQ.3) THEN
        IDLE0=60
        INEL0=20
        IPG0=27
      ENDIF 


      CFLAG=1
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
210     I0=I0+NDIM
2010  FORMAT(1X,F20.15,5X,3F20.15)

220   CALL NI01(VKPG,VNI)
      IF (M.LT.2) RETURN
      I1=(NDIM+1)*INEL*IPG
      WRITE(MP,2020) (VNI(I),I=1,I1)
2020  FORMAT(/'FUNCTIONS N AND DERIVATIVES'/ (1X,8E12.5))
      RETURN


C----------COMPUTE ELEMENT STIFFNESS MATRIX
C
C--------INITIALIZE VKE
300   DO I=1,IKE
        VKE(I)=ZERO
      ENDDO

      IF (NDIM.EQ.2) THEN
        IF (IAX.EQ.0) THEN
          IMATD=3
        ELSE
          IMATD=4
        ENDIF
      ELSEIF (NDIM.EQ.3) THEN
        IMATD = 6
      ENDIF

C--------FORM MATRIX D
      CALL D20(VPREE,VDE)
      IF (M.GE.2) WRITE(MP,2030) (VDE(I),I=1,IMATD*IMATD)
2030  FORMAT(/' MATRIX D'/1X,9E12.5)

C----------LOOP OVER G.P
      IDECL=(NDIM+1)*INEL
      I1=1+INEL
      DO IG=1,IPG
C-----------EVALUATE THE JACOBIAN, ITS INVERSE AND ITS DETERMINANT
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF (DETJ.LT.EPS) WRITE(MP,2040)IEL,IG,DETJ
2040  FORMAT(' *** ELEM ',I5,' G.P. ',I3,' DET(J)=',E12.5)
        IF (M.GE.2) WRITE(MP,2050)VJ,VJ1,DETJ
2050  FORMAT(/' JACOBIAN=',4E12.5 / ' J INVERS=',
     1       4E12.5/' DETJ=',E12.5)
C-----------PERFORM D*COEF
        C=VCPG(IG)*DETJ
        IF (IAX.EQ.0) THEN
          DO I=1,IMATD*IMATD
             VDE1(I)=VDE(I)*C
          ENDDO
        ELSE
          R=0.0
          DO I=1,INEL
            R=R+VNI(IDECL*(IG-1)+I)*VCORE((I-1)*NDIM+2)
          ENDDO
          C=C*2*PI*R
          DO I=1,IMATD*IMATD
            VDE1(I)=VDE(I)*C
          ENDDO
        ENDIF  
C-----------FORM MATRIX B
        CALL DNIDX(VNI(I1),VJ1,NDIM,INEL,VNIX)
        IF (M.GE.2) WRITE(MP,2060)(VNIX(I),I=1,INEL*NDIM)
2060  FORMAT(/' VNIX'/(1X,10E12.5))
      IF (NDIM.EQ.2) THEN
        IF (IAX.EQ.0) THEN
          CALL B02(VNIX,INEL,VBE)
          IF (M.GE.2) WRITE(MP,2070)(VBE(I),I=1,3*NDIM*INEL)
        ELSE
          CALL BAXI(VNIX,VNI((IG-1)*IDECL+1),INEL,R,VBE)
          IF (M.GE.2) WRITE(MP,2070)(VBE(I),I=1,4*NDIM*INEL)
        ENDIF
      ELSEIF (NDIM.EQ.3) THEN
        CALL B3D(VNIX,INEL,VBE)
        IF (M.GE.2) WRITE(MP,2070)(VBE(I),I=1,6*NDIM*INEL)
      ENDIF 
2070  FORMAT(/' MATRIX B'/(1X,10E12.5))

        CALL BTDB(VKE,VBE,VDE1,IDLE,IMATD,NSYM)
        I1=I1+IDECL

C.........End loop over g.p.
      ENDDO

C...........Perform multiplication with Rayleigh parameter.
      DO I=1,IKE
        VIE(I)=VKE(I)*DCMPLX(1.0D0,VPREE(4)*W)
      ENDDO

C........If nonsymmetrical problem, copy triangle to lower triangle
      IF (NSYM.NE.0) THEN
C..........Copy to "help-array"
        DO I=1,IDLE*(IDLE+1)/2
          VIE1(I)=VIE(I)
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
              VIE(PNTCO+ROW-1)=VIE1(PNTCO1+ROW-1)
              PNTRW1=PNTRW1+1
            ELSE
              PNTRW1=PNTRW1+ROW-1
              VIE(PNTCO+ROW-1)=VIE1(PNTCO1+PNTRW1-1)         
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
500   DO I=1,IKE
        VKE(I)=ZERO
      ENDDO

C------------LOOP OVER THE G.P
      IDIM1=NDIM-1
      IDECL=(NDIM+1)*INEL
      I1=1+INEL
      I2=0
      DO IG=1,IPG
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        C=VCPG(IG)*DETJ*VPREE(3)

        IF (IAX.NE.0) THEN
          R=0.0
          DO I=1,INEL
            R=R+VNI(IDECL*(IG-1)+I)*VCORE((I-1)*NDIM+2)
          ENDDO
          C=C*2*PI*R
        ENDIF  

C------------ACCUMULATE MASS TERMS
        IDL=0
        DO J=1,INEL
          JJ=I2+J
          J0=1+IDL*(IDL+1)/2
          DO 530 I=1,J
            II=I2+I
            C=VNI(II)*VNI(JJ)*C
            VKE(J0)=VKE(J0)+C
            IF (NDIM.EQ.1) GO TO 530
            J1=J0+IDL+2
            DO II=1,IDIM1
              VKE(J1)=VKE(J1)+C
c              J1=J1+J1+1
              J1=J1+IDL+3
            ENDDO
530         J0=J0+NDIM
          IDL=IDL+NDIM
        ENDDO
        I1=I1+IDECL
        I2=I2+IDECL
      ENDDO

C...........Perform multiplication with Rayleigh parameter, and -W2
      DO I=1,IKE
        VIE(I)=-W2*VKE(I)*DCMPLX(1.0D0,VPREE(5)/W)
      ENDDO

C........If nonsymmetrical problem, copy triangle to lower triangle
      IF (NSYM.NE.0) THEN
C..........Copy to "help-array"
        DO I=1,IDLE*(IDLE+1)/2
          VIE1(I)=VIE(I)
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
              VIE(PNTCO+ROW-1)=VIE1(PNTCO1+ROW-1)
              PNTRW1=PNTRW1+1
            ELSE
              PNTRW1=PNTRW1+ROW-1
              VIE(PNTCO+ROW-1)=VIE1(PNTCO1+PNTRW1-1)         
            END IF
          END DO
        END DO
      ENDIF

      RETURN

C    
C-------------EVALUATE THE ELEMENT RESIDUAL
600   CALL D02(VPREE,VDE)

C--------------INITIALIZE THE RESIDUAL VECTOR
      DO ID=1,IDLE
        VFE(ID)=ZERO
      ENDDO                               

C-------------LOOP OVER THE G.P.
      I1=1+INEL
      DO IG=1,IPG
C---------EVALUATE THE JACOBIAN
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
C---------EVALUATE FUNCTIONS D(NI)/DX
        CALL DNIDX(VNI(I1),VJ1,NDIM,INEL,VNIX)

C---------EVALUATE STRAINS AND STRESSES
        EPSX=ZERO
        EPSY=ZERO
        GAMXY=ZERO
        ID=1
        DO IN=1,INEL
          UN=VDLE(ID)
          VN=VDLE(ID+1)
          C1=VNIX(IN)
          IN1=IN+INEL
          C2=VNIX(IN1)
          EPSX=EPSX+C1*UN
          EPSY=EPSY+C2*VN
          GAMXY=GAMXY+C1*VN+C2*UN
          ID=ID+2
        ENDDO

        C1=VCPG(IG)*DETJ
        C2=VDE(2)*C1
        C3=VDE(9)*C1
        C1=VDE(1)*C1
        SIGX=C1*EPSX+C2*EPSY
        SIGY=C2*EPSX+C1*EPSY
        TAUXY=C3*GAMXY

C-------------FORM THE RESIDUAL
        ID=1
        DO IN=1,INEL
          C1=VNIX(IN)
          IN1=IN+INEL
          C2=VNIX(IN1)
          VFE(ID)=VFE(ID)+C1*SIGX+C2*TAUXY
          VFE(ID+1)=VFE(ID+1)+C2*SIGY+C1*TAUXY
          ID=ID+2
        ENDDO
        I1=I1+(NDIM+1)*INEL
      ENDDO

      RETURN


C
C------------EVALUATE VOLUMIC FORCES, FX, FY PER UNIT VOLUME
C------------(FOR GRAVITY FX=0, FY=-VPREE(4)
C
700   FX=ZERO
      FY=-VPREE(4)
      DO I=1,16
        VFE(I)=ZERO
      ENDDO

      I1=1
      IDECL=(NDIM+1)*INEL
      DO IG=1,IPG
        CALL JACOB(VNI(I1+INEL),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        DX=VCPG(IG)*DETJ
        DY=DX*FY
        DX=DX*FX
        I2=I1
        I3=1
        DO IN=1,INEL
          VFE(I3)=VFE(I3)+DX*VNI(I2)
          VFE(I3+1)=VFE(I3+1)+DY*VNI(I2)
          I2=I2+1
          I3=I3+2
        ENDDO
        I1=I1+IDECL
      ENDDO

      RETURN

C
C-------------EVALUATE AND PRINT STRESSES AT G.P.
800   WRITE(MP,2080)IEL
2080  FORMAT(//' STRESSES IN ELEMENT:',I5//
     1 '  P.G.',7X,'X',11X,'Y',9X,'EPSX',8X,'EPSY',7X,'GAMXY',8X,'SIGX',
     2 8X,'SIGY',7X,'TAUXY',8X,'TETA'/ 71X ,'SIG1',8X,'SIG2',7X,'TAUMAX'
     3 /)

C----------FORM THE MATRIX D
      CALL D02(VPREE,VDE)
C----------LOOP OVER G.P.
      I1=1+INEL
      I2=0
      DO IG=1,IPG
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
        DO IN=1,INEL
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
          ID=ID+2
        ENDDO

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
     1      TETA,SIG1,SIG2,TAUMAX
2090  FORMAT(1X,I5,8E12.5,5X,F5.1/66X,3E12.5)
        I2=I2+(NDIM+1)*INEL
        I1=I1+(NDIM+1)*INEL
      ENDDO

      RETURN

900   CONTINUE
      RETURN


1000  IF(IFL.NE.IFL1) THEN
C ............gauss points

        IF (NDIM.EQ.2) THEN
C.............still IPG temporart til 1dim integrasjon
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

        ELSEIF (NDIM.EQ.3) THEN
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
            DO I=1,27
              VKPG1(I)=VKPG(I)
            ENDDO
            IF(IFL.EQ.2) II=-1.
            IF(IFL.EQ.3) II=1.
            DO I=1,25,3
              VKPG1(I)=II
            ENDDO
          ELSEIF(IFL.GT.3.AND.IFL.LT.6) THEN
            DO I=1,9
              VKPG1(I)=VKPG(I)
            ENDDO
            DO I=1,9
              VKPG1(9+I)=VKPG(27+I)
            ENDDO
            DO I=1,9
              VKPG1(18+I)=VKPG(54+I)
            ENDDO
            IF(IFL.EQ.4) II=-1.
            IF(IFL.EQ.5) II=1.
            DO I=1,9
C              VKPG1(IVK(I))=II
              VKPG1((I-1)*3+2)=II
            ENDDO
          ELSEIF(IFL.GE.6) THEN
            DO I=1,9
              VKPG1((I-1)*3+1)=VKPG((I-1)*9+1)
            ENDDO
            DO I=1,9
              VKPG1((I-1)*3+2)=VKPG((I-1)*9+2)
            ENDDO
            DO I=1,9
              VKPG1((I-1)*3+3)=VKPG((I-1)*9+3)
            ENDDO
            IF(IFL.EQ.6) II=-1.
            IF(IFL.EQ.7) II=1.
            DO I=1,9
              VKPG1((I-1)*3+3)=II
            ENDDO
          ENDIF

        ENDIF
      ENDIF

      IFL1=IFL

C..............funksjoner etc.
      CALL NI01(VKPG1,VNI)

      DO I=1,INEL*(NDIM+1)*IPG
        IF(ABS(VNI(I)).LT.(1.D-10)) VNI(I)=ZERO
      ENDDO

C.............source vector to zero
      DO I=1,IDLE
        VSRE(I)=CZERO
      ENDDO

C-------------Find the velocities.      
      IF (NDIM.EQ.2) THEN
        V(1)=VNDSRC((IFREQ-1)*4+1)
        V(2)=VNDSRC((IFREQ-1)*4+2)
        V(3)=VNDSRC((IFREQ-1)*4+3)
      ELSEIF (NDIM.EQ.3) THEN
        V(1)=VNDSRC((IFREQ-1)*10+1)
        V(2)=VNDSRC((IFREQ-1)*10+2)
        V(3)=VNDSRC((IFREQ-1)*10+3)
        V(4)=VNDSRC((IFREQ-1)*10+4)
        V(5)=VNDSRC((IFREQ-1)*10+5)
        V(6)=VNDSRC((IFREQ-1)*10+6)
        V(7)=VNDSRC((IFREQ-1)*10+7)
        V(8)=VNDSRC((IFREQ-1)*10+8)
      ENDIF

      IDECL=(NDIM+1)*INEL
      INI=0
      DO 1030 IG=1,IPG
C.............. calculate the jacobian matrix VJ , determinant
C               will return as zero ,inverse not calculated
        I1=INI+INEL+1
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF (NDIM.EQ.2) THEN
          IF(IFL.GT.3) THEN
            DETJ=VJ(1)*VJ(1)+VJ(3)*VJ(3)
          ELSE
            DETJ=VJ(2)*VJ(2)+VJ(4)*VJ(4)
          ENDIF
        ELSEIF (NDIM.EQ.3) THEN
          IF(IFL.LT.4) THEN
            DETJ=(VJ(5)*VJ(9)-VJ(8)*VJ(6))**2
            DETJ=DETJ+(VJ(8)*VJ(3)-VJ(2)*VJ(9))**2
            DETJ=DETJ+(VJ(2)*VJ(6)-VJ(5)*VJ(3))**2
          ELSEIF(IFL.GT.3.AND.IFL.LT.6) THEN
            DETJ=(VJ(4)*VJ(9)-VJ(7)*VJ(6))**2
            DETJ=DETJ+(VJ(7)*VJ(3)-VJ(1)*VJ(9))**2
            DETJ=DETJ+(VJ(1)*VJ(6)-VJ(4)*VJ(3))**2
          ELSEIF(IFL.GT.5) THEN
            DETJ=(VJ(4)*VJ(8)-VJ(7)*VJ(5))**2
            DETJ=DETJ+(VJ(7)*VJ(2)-VJ(1)*VJ(8))**2
            DETJ=DETJ+(VJ(1)*VJ(5)-VJ(4)*VJ(2))**2
          ENDIF
        ENDIF
        DETJ=SQRT(DETJ)

C---------------COMPUTE THE SOURCE TERM IN THE GAUSS POINT
        IF (NDIM.EQ.2) THEN
          IF (IFL.EQ.2) THEN
            D=VNI(INI+7)*V(1)+VNI(INI+8)*V(2)+VNI(INI+1)*V(3)
          ELSEIF (IFL.EQ.3) THEN
            D=VNI(INI+3)*V(1)+VNI(INI+4)*V(2)+VNI(INI+5)*V(3)
          ELSEIF (IFL.EQ.4) THEN
            D=VNI(INI+1)*V(1)+VNI(INI+2)*V(2)+VNI(INI+3)*V(3)
          ELSEIF (IFL.EQ.5) THEN
            D=VNI(INI+5)*V(1)+VNI(INI+6)*V(2)+VNI(INI+7)*V(3)
          ENDIF
        ELSEIF (NDIM.EQ.3) THEN
          IF (IFL.EQ.2) THEN
            D=VNI(INI+1)*V(1)+VNI(INI+9)*V(2)+VNI(INI+13)*V(3)+
     1       VNI(INI+20)*V(4)+VNI(INI+19)*V(5)+VNI(INI+12)*V(6)+
     2       VNI(INI+7)*V(7)+VNI(INI+8)*V(8)
          ELSEIF (IFL.EQ.3) THEN
            D=VNI(INI+3)*V(1)+VNI(INI+4)*V(2)+VNI(INI+5)*V(3)+
     1       VNI(INI+11)*V(4)+VNI(INI+17)*V(5)+VNI(INI+16)*V(6)+
     2       VNI(INI+15)*V(7)+VNI(INI+10)*V(8)
          ELSEIF (IFL.EQ.4) THEN
            D=VNI(INI+1)*V(1)+VNI(INI+2)*V(2)+VNI(INI+3)*V(3)+
     1       VNI(INI+10)*V(4)+VNI(INI+15)*V(5)+VNI(INI+14)*V(6)+
     2       VNI(INI+13)*V(7)+VNI(INI+9)*V(8)
          ELSEIF (IFL.EQ.5) THEN
            D=VNI(INI+5)*V(1)+VNI(INI+6)*V(2)+VNI(INI+7)*V(3)+
     1       VNI(INI+12)*V(4)+VNI(INI+19)*V(5)+VNI(INI+18)*V(6)+
     2       VNI(INI+17)*V(7)+VNI(INI+11)*V(8)
          ELSEIF (IFL.EQ.6) THEN
            D=VNI(INI+1)*V(1)+VNI(INI+8)*V(2)+VNI(INI+7)*V(3)+
     1       VNI(INI+6)*V(4)+VNI(INI+5)*V(5)+VNI(INI+4)*V(6)+
     2       VNI(INI+3)*V(7)+VNI(INI+2)*V(8)
          ELSEIF (IFL.EQ.7) THEN
            D=VNI(INI+13)*V(1)+VNI(INI+14)*V(2)+VNI(INI+15)*V(3)+
     1       VNI(INI+16)*V(4)+VNI(INI+17)*V(5)+VNI(INI+18)*V(6)+
     2       VNI(INI+19)*V(7)+VNI(INI+20)*V(8)
          ENDIF
        ENDIF

C............. compute the weight
        D=VCPG1(IG)*DETJ*D

C.................if axi-symmetric problem
        IF(IAX.EQ.1) THEN
          Y=ZERO
          DO  J=1,INEL
            Y=Y+VCORE(2*J)*VNI(INI+J)
          END DO
          D=D*2.*PI*Y
        END IF


C............. terms of the source matrix, frame only
        DO J=1,INEL
          I2=INI+J
          IF (NDIM.EQ.2) THEN
            VSRE(2*(J-1)+1)=VSRE(2*(J-1)+1)+
     1           VNI(I2)*D*REAL(VNDSRC((IFREQ-1)*4+4))
            VSRE(2*(J-1)+2)=VSRE(2*(J-1)+2)+
     1           VNI(I2)*D*IMAG(VNDSRC((IFREQ-1)*4+4))
          ELSEIF (NDIM.EQ.3) THEN
            VSRE(3*(J-1)+1)=VSRE(3*(J-1)+1)+
     1           VNI(I2)*D*REAL(VNDSRC((IFREQ-1)*10+9))
            VSRE(3*(J-1)+2)=VSRE(3*(J-1)+2)+
     1           VNI(I2)*D*IMAG(VNDSRC((IFREQ-1)*10+9))
            VSRE(3*(J-1)+3)=VSRE(3*(J-1)+3)+
     1           VNI(I2)*D*REAL(VNDSRC((IFREQ-1)*10+10))
          ENDIF
        ENDDO

1030    INI=INI+IDECL

C................... tilbakestill IGP
      IPG=IPG1

      RETURN
      END

