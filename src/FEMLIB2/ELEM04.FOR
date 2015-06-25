      SUBROUTINE ELEM04(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1      VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C********************************************************************C
C ELEM04.FOR                                                         C
C  3/1-1995     NTH                                                  C
C  Tonni F. Johansen  & Ulf R. Kristiansen                           C
C********************************************************************C
C  Quadratic element for porous materials using equivalent fluid     C
C  model 1,2 or 3 dimensions:                                        C
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
C     ICODE.EQ.8  print gradients     .....not written.....NB!!      C
C     ICODE.EQ.9  impedance matrix(VKE)  .....not written.....NB!!   C
C     ICODE.EQ.10 source matrix (VKE)    .....not written.....NB!!   C
C  Element properties                                                C
C     VPREE(1)  Model used to describe equivalent fluid              C
C       VPREE(1).EQ.1  Craggs' model                                 C
C       VPREE(1).EQ.2  Johnson & Allard's model                      C
C       VPREE(1).EQ.10 Craggs' model, directly after paper           C
C                                                                    C
C     If VPREE(1).EQ.1 (Craggs' model)                               C
C       VPREE(2)  Pore fluid density                                 C
C       VPREE(3)  Pore fluid wave velocity                           C
C       VPREE(4)  Flow resistivity (Rayls/m)                         C
C       VPREE(5)  Porosity                                           C
C       VPREE(6)  Structural factor, Ks                              C
C                                                                    C
C     If VPREE(1).EQ.2 (Johnson-Allard's model)                      C
C       VPREE(2)  Density of fluid                                   C
C       VPREE(3)  Input: Wave velocity of fluid,                     C
C                 Internal: Bulk Modulus.                            C
C       VPREE(4)  Prandtl's number                                   C
C       VPREE(5)  Viscosity of fluid                                 C
C       VPREE(6)  Gamma for the fluid, type of gas (for liquid ??)   C
C       VPREE(7) Porosity                                            C
C       VPREE(8) Tortuosity,Structure factor                         C
C       VPREE(9) Flow resistivity                                    C
C       VPREE(10) Lambda, parameter for viscous losses               C
C       VPREE(11) Lambda-dot, parameter for thermal losses.          C
C                                                                    C
C     If VPREE(10).EQ.1 (Craggs' model), directly after paper        C
C       VPREE(2)  Pore fluid density                                 C
C       VPREE(3)  Pore fluid wave velocity                           C
C       VPREE(4)  Flow resistivity (Rayls/m)                         C
C       VPREE(5)  Porosity                                           C
C       VPREE(6)  Structural factor, Ks                              C
C                                                                    C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
c$DECLARE      
      
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

      SAVE VCPG,VKPG,VNI

      INTEGER I,IKE,INI,IG,INI1,J,IK,I1,I2,IJ
      REAL*8 RH2,WR2,T,K,K2
      COMPLEX*16 FAK1,FAK2,RHOD,KDOT,GJ
      REAL*8 ZERO,EPS,PI,DETJ,COEF,X,Y,XXYY,Y2,X2,XS2,YS2,XY2,A,B,C,c0

      INTEGER IDECL,INI0,IC,I0,IN
      INTEGER COL,ROW,PNTCO,PNTCO1,PNTRW1
      COMPLEX*16 CC,CZERO
      
C------------Number of Gauss points in ksi,eta,dzeta direction
      DATA IPGKED/3,3,3/
C---------------------------------------------------
      DATA ZERO/0.0D0/,EPS/1.D-6/,CZERO/(0.0,0.0)/
      DATA PI/3.141592654/

      IKE=IDLE*(IDLE+1)/2
    
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
        

C........Choose model of the equivalnt fluid
      IF (NINT(VPREE(1)).EQ.1) THEN
C........Craggs' model

c        RH2=VPREE(2)*VPREE(2)
c        WR2=W2*RH2
c        T=WR2+VPREE(4)*VPREE(4)
c        FAK1=CMPLX(WR2/T,W*VPREE(2)*VPREE(4)/T)
        FAK1=CMPLX(0.,VPREE(2)*W)/
     1                   CMPLX(VPREE(4),VPREE(2)*VPREE(6)/VPREE(5)*W)
      ELSEIF (NINT(VPREE(1)).EQ.2) THEN
C........Johnson-Allard's model

C---------REPLACE SOUND VELOCITY BY BULKMODULUS FOR FLUID.
C        VPREE(2)=VPREE(2)*VPREE(2)*VPREE(1)

C------------EVALUATE THE FUNCTION Gj FOR VISCOUS LOSSES
        CALL GJ04(VPREE,1,GJ)
        RHOD=VPREE(8)*VPREE(2)*(1/VPREE(7)+VPREE(9)/
     1         CMPLX(0.,VPREE(2)*W*VPREE(8))*GJ)

        FAK1=VPREE(2)/RHOD
      ELSEIF (NINT(VPREE(1)).EQ.10) THEN
C........Craggs' model, old implementation

        RH2=VPREE(2)*VPREE(2)
        WR2=W2*RH2
        T=WR2+VPREE(4)*VPREE(4)
        FAK1=CMPLX(WR2/T,W*VPREE(2)*VPREE(4)/T)
      ENDIF

      DO I=1,IKE
        VIE(I)=CMPLX(VKE(I),0.)*FAK1
      END DO


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
            VKE(IK)=VKE(IK)+VNI(I1)*VNI(I2)*COEF
520     CONTINUE

C----------NEXT G.P.
        INI=INI+(NDIM+1)*INEL
530   CONTINUE

C........Choose model of the equivalnt fluid
      IF (NINT(VPREE(1)).EQ.1) THEN
C........Craggs' model
C        RH2=VPREE(2)*VPREE(2)
C        WR2=W2*RH2
C        T=WR2+VPREE(4)*VPREE(4)
C        FAK1=CMPLX(WR2/T,W*VPREE(2)*VPREE(4)/T)
        FAK1=CMPLX(0.,VPREE(2)*W)/
     1                   CMPLX(VPREE(4),VPREE(2)*VPREE(6)/VPREE(5)*W)

C nb LOOKMHERE TO CHEXK
        K=W/VPREE(3)
        K2=K*K
C        FAK2=FAK1*CMPLX(-K2*VPREE(6)*VPREE(5),
C     1             K*VPREE(4)*VPREE(5)/(VPREE(2)*VPREE(3)))
        FAK2=FAK1*CMPLX(-K2*VPREE(6),
     1             K*VPREE(4)*VPREE(5)/(VPREE(2)*VPREE(3)))

      ELSEIF (NINT(VPREE(1)).EQ.2) THEN
C........Johnson-Allard's model

C---------REPLACE SOUND VELOCITY BY BULKMODULUS FOR FLUID.
        c0=vpree(3)
        VPREE(3)=VPREE(3)*VPREE(3)*VPREE(2)

C------------EVALUATE THE FUNCTION Gj FOR VISCOUS LOSSES
        CALL GJ04(VPREE,1,GJ)
        RHOD=VPREE(8)*VPREE(2)*(1/VPREE(7)+VPREE(9)/
     1         CMPLX(0.,VPREE(2)*W*VPREE(8))*GJ)

C------------EVALUATE THE FUNCTION Gj FOR VISCOUS LOSSES
        CALL GJ04(VPREE,2,GJ)
        
c        KDOT=VPREE(3)/VPREE(7)/(VPREE(6)-(VPREE(6)-1)/
c     1     (1-CMPLX(0.,1.)*8*VPREE(5)/
c     2     (VPREE(11)*VPREE(11)*VPREE(4)*W*VPREE(2))*GJ))
        KDOT=VPREE(3)/VPREE(7)/(VPREE(6)-(VPREE(6)-1)/
     1     (1-CMPLX(0.,1.)*8*VPREE(5)/
     2     (VPREE(11)*VPREE(11)*VPREE(4)*W*VPREE(2))*GJ))


C        VELDOT2=(KDOT/RHOD)
        FAK1=VPREE(2)/RHOD
        FAK2=-FAK1*W2*RHOD/KDOT
        vpree(3)=c0

      ELSEIF (NINT(VPREE(1)).EQ.10) THEN
C........Craggs' model, old implementation
        RH2=VPREE(2)*VPREE(2)
        WR2=W2*RH2
        T=WR2+VPREE(4)*VPREE(4)
        FAK1=CMPLX(WR2/T,W*VPREE(2)*VPREE(4)/T)

        K=W/VPREE(3)
        K2=K*K
        FAK2=FAK1*CMPLX(-K2*VPREE(6)*VPREE(5),
     1             K*VPREE(4)*VPREE(5)/(VPREE(2)*VPREE(3)))

      ENDIF

      DO I=1,IKE
        VIE(I)=CMPLX(VKE(I),0.)*FAK2
      END DO
 
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
600   CONTINUE
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

C------------Find equivalent density of material

      IF (NINT(VPREE(1)).EQ.1) THEN
C........Craggs' model
        RHOD=CMPLX(VPREE(2)*VPREE(6)/VPREE(5),-VPREE(4)/W)


      ELSEIF (NINT(VPREE(1)).EQ.2) THEN
C........Johnson-Allard's model

C---------REPLACE SOUND VELOCITY BY BULKMODULUS FOR FLUID.
        c0=vpree(3)
        VPREE(3)=VPREE(3)*VPREE(3)*VPREE(2)

C------------EVALUATE THE FUNCTION Gj FOR VISCOUS LOSSES
        CALL GJ04(VPREE,1,GJ)
        RHOD=VPREE(8)*VPREE(2)*(1/VPREE(7)+VPREE(9)/
     1         CMPLX(0.,VPREE(2)*W*VPREE(8))*GJ)

        vpree(3)=c0

      ENDIF
C--------------Finished evaluate equivalent density of material




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
          VJC(I)=CC/CMPLX(0.,-W)/RHOD
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
C
C ----------EVALUATE THE IMPEDANCE MATRIX
900   CONTINUE
      RETURN


C-----------EVALUATE THE SOURCE MATRIX
1000  CONTINUE
      RETURN


      END

