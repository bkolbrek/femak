      SUBROUTINE ELEM05(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1      VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C********************************************************************C
C ELEM05.FOR                                                         C
C  23/9-1994     NTH                                                 C
C  Tonni F. Johansen  & Ulf R. Kristiansen                           C
C********************************************************************C
C  Quadratic element for porous materials using Biot-theory          C
C  2, 3 dimensions or axisymmetry                                    C
C     2 dimension: 8 nodes isoparametric element                     C
C     axisymmetry: 8 nodes isoparametric element                     C
C     3 dimension: 20 nodes isoparametric element                    C
C  Number of integration points: 3 in each direction                 C
C  Number of degrees of freedom per node: 4 or 6                     C
C  Element matrix or vector formed by this subbprogram               C
C  according to icode value :                                        C
C     ICODE.EQ.1  return of parameters                               C
C     ICODE.EQ.2  evaluate interpolation functions and               C
C                 numerical integration coeffisients                 C
C     ICODE.EQ.3  stiffness matrix (VKE)                             C
C     ICODE.EQ.4  tangent matrix (VKE)  .....not written.....NB!!    C
C     ICODE.EQ.5  mass matrix (VIE)                                  C
C     ICODE.EQ.6  K.U. product (VFE)                                 C
C     ICODE.EQ.7  element load (VFE)  .....not written.....NB!!      C
C     ICODE.EQ.8  print gradients                                    C
C     ICODE.EQ.9  impedance matrix(VIE)                              C
C     ICODE.EQ.10 source matrix (VSRE)                               C
C     ICODE.EQ.11 loss matrix (VIE)                                  C
C  Element properties                                                C
C     VPREE(1)  Density of fluid                                     C
C     VPREE(2)  Input: Wave velocity of fluid,                       C
C               Internal: Bulk Modulus.                              C
C     VPREE(3)  Prandtl's number                                     C
C     VPREE(4)  Viscosity of fluid                                   C
C     VPREE(5)  Gamma for the fluid, type of gas (for liquid ??)     C
C     VPREE(6)  Density of frame.                                    C
C     VPREE(7)  Shear modulus of frame, real part.                   C
C     VPREE(8)  Shear modulus of frame, imaginary part.              C
C     VPREE(9)  Poisson ratio of frame, real part.                   C
C     VPREE(10) Poisson ratio of frame, imaginary part.              C
C     VPREE(11) Porosity                                             C
C     VPREE(12) Tortuosity                                           C
C     VPREE(13) Flow resistivity                                     C
C     VPREE(14) Lambda, parameter for viscous losses                 C
C     VPREE(15) Lambda-dot, parameter for thermal losses.            C
C********************************************************************C
C Theory and implementation of this element is described in FEMAK    C
C programmer's manual and user's manual                              C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)


      include 'femak.cmn'
c   $DECLARE      
      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VKE1(1),VFE(1)
      COMPLEX*16 VIE(1),VIE1(1),VSRE(1),VNDSRC(1),VPRSG(1),VDLEC(1)
      INTEGER KTYPRSG(1)
      

C-----------CHARACTERISTIC DIMENSIONS OF THE ELEMENT
C----------(VALID  UP TO 3 DIMENSIONS)
C     DIMENSION VCPG(IPG),VKPG(NDIM*IPG),XYZ(NDIM),VDEC1(IMATD**2)
      REAL*8 VCPG(27),VKPG(81)
      COMPLEX*16 VDEC1(36)
C     DIMENSION VBE(IMATB*IDLE),VDEC(IMATD**2),VJ(NDIM*NDIM),VJ1(NDIM*NDIM)
      REAL*8 VBE(360),VJ(9),VJ1(9)
C NB!! change VBE(48) to VBE(64) check for troubles.!!!

      COMPLEX*16 VDEC(36)
C     DIMENSION VNIX(INEL*NDIM),VNI((1+NDIM)*INEL*IPG),IPGKED(NDIM)
      REAL*8 VNIX(60),VNI(2160)
      REAL*8 VKPG1(81),VCPG1(27)
      INTEGER IPGKED(3)
      INTEGER IVK(9)                              

      SAVE VCPG,VKPG,VNI


      INTEGER I,J,II,JJ,JJJ,I0,IG,I1,IDIM1,IDECL,I2,J0,IDL,J1,IPG1,INI
      INTEGER IKE,IIUP,ID,IDN,IN,IN1
      REAL*8 PI,DETJ,EPS,X,Y,Z,ZERO,R
      COMPLEX*16 BIOTM,CZERO,GJ,TB(9),BTDBC(9),D,C
      COMPLEX*16 V(8)
      COMPLEX*16 UX,UY,UZ,WX,WY,WZ,XN,YN,ZN,EPSX,EPSY,EPSZ
      COMPLEX*16 GAMXY,GAMYZ,GAMZX,DIVW,C1,C2,C3,C4
      COMPLEX*16 SIGX,SIGY,SIGZ,TAUXY,TAUYZ,TAUZX,PORP
      INTEGER COL,ROW,PNTCO,PNTCO1,PNTRW1


C------------Number of Gauss points in ksi,eta,dzeta direction
      DATA IPGKED/3,3,3/
c      DATA IPGKED/6,6,6/
C---------------------------------------------------
      DATA EPS/1.D-6/,CZERO/(0.0,0.0)/,ZERO/0/
      DATA IVK/2,5,8,29,32,35,56,59,62/
      DATA PI/3.141592654/
      
      
      IKE=IDLE*(IDLE+1)/2

C-------------CHOOSE FUNCTION TO BE EXECUTED
      GO TO (100,200,300,400,500,600,700,800,900,1000,1100),ICODE

C
C------------RETURN ELEMENT PARAMETERS IN COMMON /RGDT/
100   IF (NDIM.EQ.2) THEN
        IDLE0=32
        INEL0=8
        IPG0=9
      ELSEIF (NDIM.EQ.3) THEN
        IDLE0=120
        INEL0=20
        IPG0=27
      ENDIF 

C---------REPLACE SOUND VELOCITY BY BULKMODULUS FOR FLUID.
      VPREE(2)=VPREE(2)*VPREE(2)*VPREE(1)
      CFLAG=1

      RETURN

C
C-----------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C-----------THEIR DERIVATIVES AT G.P.

200   CALL GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
      IF(M.GE.2) THEN
        WRITE(MP,2000)IPG
2000  FORMAT(/I5,'  GAUSS POINTS'/10X,'VCPG',25X,'VKPG')
        I0=1
        DO IG=1,IPG
          I1=I0+NDIM-1
          WRITE(MP,2010)VCPG(IG),(VKPG(I),I=I0,I1)
          I0=I0+NDIM
        ENDDO
      ENDIF
2010  FORMAT(1X,F20.15,5X,3F20.15)

c      CALL NI02(VKPG,VNI)
      CALL NI01(VKPG,VNI)
      IF (M.GE.2) THEN
        I1=3*INEL*IPG
        WRITE(MP,2020) (VNI(I),I=1,I1)
      ENDIF
2020  FORMAT(/'FUNCTIONS N AND DERIVATIVES'/ (1X,8E12.5))
   
      RETURN

C----------COMPUTE ELEMENT STIFFNESS MATRIX
C
C--------INITIALIZE VKE
300   DO I=1,IKE
        VIE(I)=CZERO
      ENDDO

C--------calculate coefficients. 
C------------EVALUATE THE FUNCTION Gj FOR THERMAL LOSSES
      CALL GJ05(VPREE,2,GJ)
      
      BIOTM=VPREE(2)/VPREE(11)/(VPREE(5)-(VPREE(5)-1)/
     1     (1-CMPLX(0.,1.)*8*VPREE(4)/
     2     (VPREE(15)*VPREE(15)*VPREE(3)*W*VPREE(1))*GJ))

C--------FORM MATRIX D
      CALL D05(VPREE,BIOTM,VDEC)
      IF (M.GE.2) THEN

        IF ((NDIM.EQ.2).AND.(IAX.EQ.0)) THEN
          WRITE(MP,2030) (VDEC(I),I=1,9)
        ELSEIF ((NDIM.EQ.2).AND.(IAX.NE.0)) THEN
          WRITE(MP,2031) (VDEC(I),I=1,16)
        ELSEIF (NDIM.EQ.3) THEN
          WRITE(MP,2031) (VDEC(I),I=1,36)
        ENDIF
      ENDIF
2030  FORMAT(/' MATRIX D'/1X,10E12.5)
2031  FORMAT(/' MATRIX D'/1X,10E12.5)

C----------LOOP OVER G.P
      IDIM1=NDIM-1
      IDECL=(NDIM+1)*INEL
      I1=1+INEL
      I2=0

      DO IG=1,IPG
C-----------EVALUATE THE JACOBIAN, ITS INVERSE AND ITS DETERMINANT
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF (DETJ.LT.EPS) WRITE(MP,2040)IEL,IG,DETJ
2040  FORMAT(' *** ELEM ',I5,' G.P. ',I3,' DET(J)=',E12.5)

        IF (M.GE.2) THEN                                 
          IF (NDIM.EQ.2) THEN
            WRITE(MP,2050)(VJ(I),I=1,4),(VJ1(I),I=1,4),DETJ
          ELSEIF (NDIM.EQ.3) THEN
            WRITE(MP,2052)(VJ(I),I=1,9),(VJ1(I),I=1,9),DETJ
          ENDIF
        ENDIF  
2050  FORMAT(/' JACOBIAN=',4E12.5 / ' J INVERS=',
     1       4E12.5/' DETJ=',E12.5)
2052  FORMAT(/' JACOBIAN=',9E12.5 / ' J INVERS=',
     1       9E12.5/' DETJ=',E12.5)

C-----------PERFORM D*COEF
C-----------If axisymmetry; find radius
        IF (IAX.NE.0) THEN
          R=0.0
          DO I=1,INEL
            R=R+VNI(IDECL*(IG-1)+I)*VCORE((I-1)*NDIM+2)
          ENDDO
        ENDIF

        D=VCPG(IG)*DETJ
        IF (IAX.NE.0) THEN
          D=D*R*2*PI
        ENDIF  
        
        IF ((NDIM.EQ.2).AND.(IAX.NE.0)) THEN
          DO I=1,16
            VDEC1(I)=VDEC(I)*D
          ENDDO
        ELSEIF ((NDIM.EQ.2).AND.(IAX.EQ.0)) THEN
          DO I=1,9
            VDEC1(I)=VDEC(I)*D
          ENDDO
        ELSEIF (NDIM.EQ.3) THEN
          DO I=1,36
            VDEC1(I)=VDEC(I)*D
          ENDDO
        ENDIF
        
        
C-----------FORM MATRIX B
        CALL DNIDX(VNI(I1),VJ1,NDIM,INEL,VNIX)
        IF (M.GE.2) WRITE(MP,2060)(VNIX(I),I=1,16)
2060  FORMAT(/' VNIX'/(1X,8E12.5))
        IF ((NDIM.EQ.2).AND.(IAX.EQ.0)) THEN
          CALL B02(VNIX,INEL,VBE)
          IF (M.GE.2) WRITE(MP,2070)(VBE(I),I=1,3*(NDIM-1)*NDIM*INEL)
        ELSEIF ((NDIM.EQ.2).AND.(IAX.NE.0)) THEN
          CALL BAXI(VNIX,VNI((IG-1)*IDECL+1),INEL,R,VBE)
          IF (M.GE.2) WRITE(MP,2070)(VBE(I),I=1,4*NDIM*INEL)
        ELSEIF (NDIM.EQ.3) THEN
          CALL B3D(VNIX,INEL,VBE)
          IF (M.GE.2) WRITE(MP,2070)(VBE(I),I=1,6*NDIM*INEL)
        ENDIF
2070  FORMAT(/' MATRIX B'/(1X,10E12.5))


        IDL=0
C----------LOOP THROUGH ALL "SMALL MATRICES", COLOUMNS
        DO J=1,INEL
          JJ=I2+J
          J0=1+IDL*(IDL+1)/2
C------------LOOP THROUGH ALL "SMALL MATRICES", LINES ON AND ABOVE THE DIAGONAL
          DO I=1,J

C--------------Assemble Kuu first.

C--------------Find BtDB for one "small matrix"
            DO II=1,NDIM
C---------------First Bt*D
              IF ((IAX.NE.0).AND.(NDIM.EQ.2)) THEN
C..................Axisymmetric version                  
                TB(1)=VBE(8*(I-1)+4*(II-1)+1)*VDEC1(1)+
     1                VBE(8*(I-1)+4*(II-1)+2)*VDEC1(2)+
     2                VBE(8*(I-1)+4*(II-1)+3)*VDEC1(3)
                TB(2)=VBE(8*(I-1)+4*(II-1)+1)*VDEC1(5)+
     1                VBE(8*(I-1)+4*(II-1)+2)*VDEC1(6)+
     2                VBE(8*(I-1)+4*(II-1)+3)*VDEC1(7)
                TB(3)=VBE(8*(I-1)+4*(II-1)+1)*VDEC1(9)+
     1                VBE(8*(I-1)+4*(II-1)+2)*VDEC1(10)+
     2                VBE(8*(I-1)+4*(II-1)+3)*VDEC1(11)
                TB(4)=VBE(8*(I-1)+4*(II-1)+4)*VDEC1(16)
              ELSEIF ((IAX.EQ.0).AND.(NDIM.EQ.2)) THEN
C..................2D version                  
                TB(1)=VBE(6*(I-1)+3*(II-1)+1)*VDEC1(1)+
     1                VBE(6*(I-1)+3*(II-1)+2)*VDEC1(2)
                TB(2)=VBE(6*(I-1)+3*(II-1)+1)*VDEC1(4)+
     1                VBE(6*(I-1)+3*(II-1)+2)*VDEC1(5)
                TB(3)=VBE(6*(I-1)+3*(II-1)+3)*VDEC1(9)
              ELSEIF (NDIM.EQ.3) THEN
C..................3D version                  
                TB(1)=VBE(18*(I-1)+6*(II-1)+1)*VDEC1(1)+
     1                VBE(18*(I-1)+6*(II-1)+2)*VDEC1(2)+
     2                VBE(18*(I-1)+6*(II-1)+3)*VDEC1(3)
                TB(2)=VBE(18*(I-1)+6*(II-1)+1)*VDEC1(7)+
     1                VBE(18*(I-1)+6*(II-1)+2)*VDEC1(8)+
     2                VBE(18*(I-1)+6*(II-1)+3)*VDEC1(9)
                TB(3)=VBE(18*(I-1)+6*(II-1)+1)*VDEC1(13)+
     1                VBE(18*(I-1)+6*(II-1)+2)*VDEC1(14)+
     2                VBE(18*(I-1)+6*(II-1)+3)*VDEC1(15)
                TB(4)=VBE(18*(I-1)+6*(II-1)+4)*VDEC1(22)
                TB(5)=VBE(18*(I-1)+6*(II-1)+5)*VDEC1(29)
                TB(6)=VBE(18*(I-1)+6*(II-1)+6)*VDEC1(36)
              ENDIF
C-----------------Then (Bt*D)*B
              IF ((NDIM.EQ.2).AND.(IAX.NE.0)) THEN
C..................Axisymmetric version                  
                DO JJ=1,2
                  BTDBC((JJ-1)*2+II)=CZERO
                  DO JJJ=1,4
                    BTDBC((JJ-1)*2+II)=BTDBC((JJ-1)*2+II)+
     1                               TB(JJJ)*VBE(8*(J-1)+4*(JJ-1)+JJJ)
                  ENDDO
                ENDDO
              ELSEIF ((IAX.EQ.0).AND.(NDIM.EQ.2)) THEN
C..................2D version                  
                DO JJ=1,2
                  BTDBC((JJ-1)*2+II)=CZERO
                  DO JJJ=1,3
                    BTDBC((JJ-1)*2+II)=BTDBC((JJ-1)*2+II)+
     1                               TB(JJJ)*VBE(6*(J-1)+3*(JJ-1)+JJJ)
                  ENDDO
                ENDDO
              ELSEIF (NDIM.EQ.3) THEN
C..................3D version                  
                DO JJ=1,NDIM
                  BTDBC((JJ-1)*3+II)=CZERO
                  DO JJJ=1,6
                    BTDBC((JJ-1)*3+II)=BTDBC((JJ-1)*3+II)+
     1                               TB(JJJ)*VBE(18*(J-1)+6*(JJ-1)+JJJ)
                  ENDDO
                ENDDO
               ENDIF
            ENDDO


            J1=J0
C...............Put BTDBC into VIE            
            DO JJ=1,NDIM
C..............Loop through colomns(JJ)
              IF (I.EQ.J) THEN
                IIUP=JJ
              ELSE
                IIUP=NDIM
              ENDIF
              DO II=1,IIUP
C..............Loop through rows(II)
                VIE(J1+II-1)=VIE(J1+II-1)+BTDBC((JJ-1)*NDIM+II)
              ENDDO
              J1=J1+IDL+JJ
            ENDDO



C--------------Point at upper left corner of Kwu and assemble Kwu if necessary.
            IF (I.NE.J) THEN
              J0=J0+NDIM
C              VIE(J0)=VIE(J0)+VNIX(J)*VNIX(I)*D*BIOTM
              IF (NDIM.EQ.2) THEN
                VIE(J0)=VIE(J0)+VNIX(J)*VNIX(I)*D*BIOTM
c                VIE(J0+1)=VIE(J0+1)+VNIX(J)*VNIX(I+INEL)*D*BIOTM
c               VIE(J0+IDL+1)=VIE(J0+IDL+1)+VNIX(J+INEL)*VNIX(I)*D*BIOTM
                IF (IAX.EQ.0) THEN  
                  VIE(J0+1)=VIE(J0+1)+VNIX(J)*VNIX(I+INEL)*D*BIOTM
                  VIE(J0+IDL+1)=VIE(J0+IDL+1)+
     1                  VNIX(J+INEL)*VNIX(I)*D*BIOTM
                  VIE(J0+IDL+2)=VIE(J0+IDL+2)+
     1                 VNIX(J+INEL)*VNIX(I+INEL)*D*BIOTM
                ELSE
                  VIE(J0+1)=VIE(J0+1)+
     1                 VNIX(J)*(VBE((I-1)*8+7)+VNIX(I+INEL))*D*BIOTM
                  VIE(J0+IDL+1)=VIE(J0+IDL+1)+
     1                  (VBE((J-1)*8+7)+VNIX(J+INEL))*VNIX(I)*D*BIOTM
                  VIE(J0+IDL+2)=VIE(J0+IDL+2)+
     1                 (VBE((J-1)*8+7)+VNIX(J+INEL))*
     2                 (VBE((I-1)*8+7)+VNIX(I+INEL))*D*BIOTM
                ENDIF
              ELSEIF (NDIM.EQ.3) THEN
                J1=J0
                DO JJ=1,NDIM
C..............Loop through colomns(JJ)
                  DO II=1,NDIM
C..............Loop through rows(II)
                    VIE(J1+II-1)=VIE(J1+II-1)+
     1                   VNIX(J+(JJ-1)*INEL)*VNIX(I+(II-1)*INEL)*D*BIOTM
                  ENDDO
                  J1=J1+IDL+JJ
                ENDDO
              ENDIF
            ELSE
C---------------Take into acount that Kuu(2,1) and Kwu is not in upper triangle
              J0=J0+NDIM
            ENDIF
          
C--------------Point at upper left corner of Kuw and assemble Kuw.
C            J0=J0+(NDIM-1)*IDL+NDIM*(NDIM+1)/2+NDIM
            IF(NDIM.EQ.2) THEN
              J0=J0+NDIM*IDL+NDIM-1
            ELSEIF(NDIM.EQ.3) THEN
              J0=J0+NDIM*IDL+NDIM
            ENDIF
            IF (NDIM.EQ.2) THEN
              VIE(J0)=VIE(J0)+VNIX(J)*VNIX(I)*D*BIOTM
c              VIE(J0+1)=VIE(J0+1)+VNIX(J)*VNIX(I+INEL)*D*BIOTM
c              VIE(J0+IDL+3)=VIE(J0+IDL+3)+VNIX(J+INEL)*VNIX(I)*D*BIOTM
c              VIE(J0+1)=VIE(J0+1)+VNIX(J)*VNIX(I+INEL)*D*BIOTM
c              VIE(J0+IDL+3)=VIE(J0+IDL+3)+VNIX(J+INEL)*VNIX(I)*D*BIOTM
c              VIE(J0+IDL+4)=VIE(J0+IDL+4)+
c     1              VNIX(J+INEL)*VNIX(I+INEL)*D*BIOTM
              IF (IAX.EQ.0) THEN  
                VIE(J0+1)=VIE(J0+1)+VNIX(J)*VNIX(I+INEL)*D*BIOTM
                VIE(J0+IDL+3)=VIE(J0+IDL+3)+
     1                  VNIX(J+INEL)*VNIX(I)*D*BIOTM
                VIE(J0+IDL+4)=VIE(J0+IDL+4)+
     1                 VNIX(J+INEL)*VNIX(I+INEL)*D*BIOTM
              ELSE
c  .....Note difference in use of dN/dy and dN/dr is 2D correct?????
                VIE(J0+1)=VIE(J0+1)+
     1                VNIX(J)*(VBE((I-1)*8+7)+VNIX(I+INEL))*D*BIOTM
                VIE(J0+IDL+3)=VIE(J0+IDL+3)+
     1                  (VBE((J-1)*8+7)+VNIX(J+INEL))*VNIX(I)*D*BIOTM
                VIE(J0+IDL+4)=VIE(J0+IDL+4)+
     1                 (VBE((J-1)*8+7)+VNIX(J+INEL))*
     2                 (VBE((I-1)*8+7)+VNIX(I+INEL))*D*BIOTM
              ENDIF
              ELSEIF (NDIM.EQ.3) THEN
                J1=J0
                DO JJ=1,NDIM
C..............Loop through colomns(JJ)
                  DO II=1,NDIM
C..............Loop through rows(II)
C..............usikker på beregning av adresse her
                    VIE(J1+II-1)=VIE(J1+II-1)+
     1                   VNIX(J+(JJ-1)*INEL)*VNIX(I+(II-1)*INEL)*D*BIOTM
                  ENDDO
                  J1=J1+IDL+NDIM+JJ
                ENDDO
            ENDIF
          
C--------------Point at upper left corner of Kww and assemble Kww.
            J0=J0+NDIM

            IF (NDIM.EQ.2) THEN
              VIE(J0)=VIE(J0)+VNIX(J)*VNIX(I)*D*BIOTM
              IF (I.NE.J) THEN
                IF (IAX.EQ.0) THEN
                  VIE(J0+1)=VIE(J0+1)+VNIX(J)*VNIX(I+INEL)*D*BIOTM
                ELSE  
                  VIE(J0+1)=VIE(J0+1)+
     1                  VNIX(J)*(VNIX(I+INEL)+VBE((I-1)*8+7))*D*BIOTM
                ENDIF  
              ENDIF

              IF (IAX.EQ.0) THEN
                VIE(J0+IDL+3)=VIE(J0+IDL+3)+VNIX(J+INEL)*VNIX(I)*D*BIOTM
                VIE(J0+IDL+4)=VIE(J0+IDL+4)+
     1              VNIX(J+INEL)*VNIX(I+INEL)*D*BIOTM
              ELSE
                VIE(J0+IDL+3)=VIE(J0+IDL+3)+
     1                (VNIX(J+INEL)+VBE((J-1)*8+7))*VNIX(I)*D*BIOTM
                VIE(J0+IDL+4)=VIE(J0+IDL+4)+
     1              (VNIX(I+INEL)+VBE((I-1)*8+7))*
     2              (VNIX(J+INEL)+VBE((J-1)*8+7))*D*BIOTM
              
              ENDIF

            ELSEIF (NDIM.EQ.3) THEN  
              J1=J0
              DO JJ=1,NDIM
C..............Loop through colomns(JJ)
                IF (I.EQ.J) THEN
                  IIUP=JJ
                ELSE
                  IIUP=NDIM
                ENDIF
                DO II=1,IIUP
C..............Loop through rows(II)
                  VIE(J1+II-1)=VIE(J1+II-1)+
     1                VNIX(J+(JJ-1)*INEL)*VNIX(I+(II-1)*INEL)*D*BIOTM
                ENDDO
                J1=J1+IDL+JJ+NDIM
              ENDDO


     
            ENDIF

            IF (NDIM.EQ.2) THEN
              J0=J0-(4*J-4.5+NDIM*0.5)*NDIM
            ELSEIF (NDIM.EQ.3) THEN
              J0=J0-(IDL*(NDIM-1)+NDIM*(NDIM+1)/2+NDIM+(J-2)*6)
            ENDIF
          
          ENDDO
          IDL=IDL+NDIM*2
          
        ENDDO
C        CALL BTDBC(VKE,VBE,VDE1,IDLE,IMATD,NSYM)
        I1=I1+IDECL
        I2=I2+IDECL
      ENDDO
C----------END OF LOOP OF GAUSS.POINTS.                  


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


500   DO I=1,IKE
        VIE(I)=CZERO
      ENDDO

C------------LOOP OVER THE G.P
      IDIM1=NDIM-1
      IDECL=(NDIM+1)*INEL
      I1=1+INEL
      I2=0
      DO IG=1,IPG
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)

C--------------COMMON FACTOR -w^2 * gauss veight * jacobi detr.
C-----------If axisymmetry; find radius
        IF (IAX.NE.0) THEN
          R=0.0
          DO I=1,INEL
            R=R+VNI(IDECL*(IG-1)+I)*VCORE((I-1)*NDIM+2)
          ENDDO
        ENDIF

        D=-W2*VCPG(IG)*DETJ
        IF (IAX.NE.0) THEN
          D=D*R*2*PI
        ENDIF  

C------------ACCUMULATE MASS TERMS
        IDL=0

C------------LOOP THROUGH ALL "SMALL MATRICES", COLOUMNS
        DO J=1,INEL
          JJ=I2+J
          J0=1+IDL*(IDL+1)/2


C------------LOOP THROUGH ALL "SMALL MATRICES", LINES ON AND ABOVE THE DIAGONAL
          DO I=1,J

C--------------Assemble Muu first.
            II=I2+I
            C=VNI(II)*VNI(JJ)*D
            VIE(J0)=VIE(J0)+C*(VPREE(6)+VPREE(11)*VPREE(1))
            IF (NDIM.NE.1) THEN
              J1=J0+IDL+2
              DO II=1,IDIM1
                VIE(J1)=VIE(J1)+C*(VPREE(6)+VPREE(11)*VPREE(1))
                J1=J1+IDL+II+2
              ENDDO
            ENDIF

C--------------Point at upper left corner of Mwu and assemble Mwu if necessary.
            IF (I.NE.J) THEN
              J0=J0+NDIM
              VIE(J0)=VIE(J0)+C*VPREE(1)
              IF (NDIM.NE.1) THEN
                J1=J0+IDL+2
                DO II=1,IDIM1
                  VIE(J1)=VIE(J1)+C*VPREE(1)
                  J1=J1+IDL+II+2
                ENDDO
              ENDIF
            ELSE
C---------------Take into acount that muu(2,1) and mwu is not in upper triangle
              J0=J0+NDIM
            ENDIF
          
C--------------Point at upper left corner of Muw and assemble Muw.
            IF(NDIM.EQ.2) THEN
              J0=J0+NDIM*IDL+NDIM-1
            ELSEIF(NDIM.EQ.3) THEN
              J0=J0+NDIM*IDL+NDIM
            ENDIF
            VIE(J0)=VIE(J0)+C*VPREE(1)
            IF (NDIM.NE.1) THEN
              J1=J0+IDL+NDIM+2
              DO II=1,IDIM1
                VIE(J1)=VIE(J1)+C*VPREE(1)
                J1=J1+IDL+II+NDIM+2
              ENDDO
            ENDIF
          
C--------------Point at upper left corner of Mww and assemble Mww.
            J0=J0+NDIM
            VIE(J0)=VIE(J0)+C*VPREE(1)*VPREE(12)/VPREE(11)
            IF (NDIM.NE.1) THEN
              J1=J0+IDL+NDIM+2
              DO II=1,IDIM1
                VIE(J1)=VIE(J1)+C*VPREE(1)*VPREE(12)/VPREE(11)
                J1=J1+IDL+II+NDIM+2
              ENDDO
            ENDIF

            IF (NDIM.EQ.2) THEN
              J0=J0-(4*J-4.5+NDIM*0.5)*NDIM
            ELSEIF (NDIM.EQ.3) THEN
              J0=J0-(IDL*(NDIM-1)+NDIM*(NDIM+1)/2+NDIM+(J-2)*6)
c              J0=J0-(6*J-1+(NDIM+1)/2.0)*NDIM
            ENDIF
          
          ENDDO
          IDL=IDL+NDIM*2
          
        ENDDO
        I1=I1+IDECL
        I2=I2+IDECL
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
600   RETURN

C
C------------EVALUATE FE
700   CONTINUE
      RETURN

C
C-------------EVALUATE AND PRINT GRADIENTS AT G.P.
800   IF (NDIM.EQ.2) THEN
        WRITE(MP,2080)IEL
      ELSE
        WRITE(MP,2081)IEL
      ENDIF

2080  FORMAT(//' STRESSES IN ELEMENT:',I5//
     1 '  P.G.',7X,'X',11X,'Y',15X,'EPSX',20X,'EPSY',19X,'GAMXY',20X,
     2 'DIVW',20X,'SIGX',20X,'SIGY',19X,'TAUXY',20X,'POR.PRES.'/)

2081  FORMAT(//' STRESSES IN ELEMENT:',I5//
     1 '  P.G.',7X,'X',11X,'Y',11X,'Z',15X,'EPSX',20X,'EPSY',20X,
     2 'EPSZ',19X,'GAMXY',20X,'GAMYZ',20X,'GAMZX',20X,'DIVW',20X/
     3 11X,'SIGX',20X,'SIGY',20X,'SIGZ',19X,'TAUXY',19X,'TAUYZ',19X,
     4 'TAUZX',19X,'POR.PRES.'/)

C--------calculate coefficients. 
C------------EVALUATE THE FUNCTION Gj FOR THERMAL LOSSES
      CALL GJ05(VPREE,2,GJ)
      
      BIOTM=VPREE(2)/VPREE(11)/(VPREE(5)-(VPREE(5)-1)/
     1     (1-CMPLX(0.,1.)*8*VPREE(4)/
     2     (VPREE(15)*VPREE(15)*VPREE(3)*W*VPREE(1))*GJ))

C--------FORM MATRIX D
      CALL D05(VPREE,BIOTM,VDEC)

C----------LOOP OVER G.P.
      I1=1+INEL
      I2=0
      DO IG=1,IPG
C---------EVALUATE THE JACOBIAN
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
C---------EVALUATE FUNCTIONS D(NI)/DX
        CALL DNIDX(VNI(I1),VJ1,NDIM,INEL,VNIX)

C---------COMPUTE STRAINS AND STRESSES AT G.P.
        IF (NDIM.EQ.2) THEN
          EPSX=CZERO
          EPSY=CZERO
          GAMXY=CZERO
          DIVW=CZERO
          X=ZERO
          Y=ZERO
        ELSE
          EPSX=CZERO
          EPSY=CZERO
          EPSZ=CZERO
          GAMXY=CZERO
          GAMYZ=CZERO
          GAMZX=CZERO
          DIVW=CZERO
          X=ZERO
          Y=ZERO
          Z=ZERO
        ENDIF
        ID=1
        IDN=1
        DO IN=1,INEL
c..........Find nodal dof's and coordiantes
          IF (NDIM.EQ.2) THEN
            UX=VDLEC(ID)
            UY=VDLEC(ID+1)
            WX=VDLEC(ID+2)
            WY=VDLEC(ID+3)
            XN=VCORE(IDN)
            YN=VCORE(IDN+1)
          ELSE
            UX=VDLEC(ID)
            UY=VDLEC(ID+1)
            UZ=VDLEC(ID+2)
            WX=VDLEC(ID+3)
            WY=VDLEC(ID+4)
            WZ=VDLEC(ID+5)
            XN=VCORE(IDN)
            YN=VCORE(IDN+1)
            ZN=VCORE(IDN+2)
          ENDIF

c............Find partial derivatives
          IF (NDIM.EQ.2) THEN
            C1=VNIX(IN)
            IN1=IN+INEL
            C2=VNIX(IN1)
            IN1=IN+I2
            C3=VNI(IN1)
          ELSE
            C1=VNIX(IN)
            IN1=IN+INEL
            C2=VNIX(IN1)
            IN1=IN1+INEL
            C3=VNIX(IN1)
            IN1=IN+I2
            C4=VNI(IN1)
          ENDIF

C...........Find strains 
          IF (NDIM.EQ.2) THEN
            EPSX=EPSX+C1*UX
            EPSY=EPSY+C2*UY
            GAMXY=GAMXY+C1*UY+C2*UX
            DIVW=DIVW+C1*WX+C2*WY
          ELSE
            EPSX=EPSX+C1*UX
            EPSY=EPSY+C2*UY
            EPSZ=EPSZ+C3*UZ
            GAMXY=GAMXY+C1*UY+C2*UX
            GAMYZ=GAMYZ+C3*UY+C2*UZ
            GAMZX=GAMZX+C1*UZ+C3*UX
            DIVW=DIVW+C1*WX+C2*WY+C3*WZ
          ENDIF

C...........Find global coordiantes of g.p.
          IF (NDIM.EQ.2) THEN
            X=X+C3*XN
            Y=Y+C3*YN
          ELSE
            X=X+C4*XN
            Y=Y+C4*YN
            Z=Z+C4*ZN
          ENDIF

          ID=ID+2*NDIM
          IDN=IDN+NDIM
        ENDDO

C--------------COMPUTE THE STRESSES
        IF (NDIM.EQ.2) THEN
          SIGX=VDEC(1)*EPSX+VDEC(2)*EPSY+BIOTM*DIVW
          SIGY=VDEC(2)*EPSX+VDEC(1)*EPSY+BIOTM*DIVW
          TAUXY=VDEC(9)*GAMXY
          PORP=-BIOTM*(EPSX+EPSY+DIVW)
        ELSE
          SIGX=VDEC(1)*EPSX+VDEC(2)*EPSY+VDEC(3)*EPSZ+BIOTM*DIVW
          SIGY=VDEC(2)*EPSX+VDEC(1)*EPSY+VDEC(3)*EPSZ+BIOTM*DIVW
          SIGZ=VDEC(2)*EPSX+VDEC(2)*EPSY+VDEC(1)*EPSZ+BIOTM*DIVW
          TAUXY=VDEC(36)*GAMXY
          TAUYZ=VDEC(36)*GAMYZ
          TAUZX=VDEC(36)*GAMZX
          PORP=-BIOTM*(EPSX+EPSY+EPSZ+DIVW)
        ENDIF


        IF (NDIM.EQ.2) THEN
          WRITE(MP,2090)IG,X,Y,EPSX,EPSY,GAMXY,DIVW,SIGX,SIGY,TAUXY,PORP
        ELSE
          WRITE(MP,2091)IG,X,Y,Z,EPSX,EPSY,EPSZ,GAMXY,GAMYZ,GAMZX,DIVW,
     Q     SIGX,SIGY,SIGY,TAUXY,TAUYZ,TAUZX,PORP
        ENDIF
2090    FORMAT(1X,I5,2E12.5,16E12.5)
2091    FORMAT(1X,I5,3E12.5,14E12.5/,14E12.5)
        I2=I2+(NDIM+1)*INEL
        I1=I1+(NDIM+1)*INEL
      ENDDO


      RETURN
C
C ----------EVALUATE THE IMPEDANCE MATRIX
900   RETURN


C-----------EVALUATE THE SOURCE MATRIX
1000  IF(IFL.NE.IFL1) THEN
C ............gauss punktene
 
C        GO TO(1011,1012,1013)NDIM
        
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

C..............functions etc.
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
      
      INI=0
      DO 1030 IG=1,IPG
C.............. calculate the jacobian matrix VJ , determinant
C               will return as zero ,inverse not calculated
        I1=INI+INEL+1
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF (NDIM.EQ.2) THEN
          IF(IFL.GT.3) THEN
            DETJ=SQRT(VJ(1)*VJ(1)+VJ(3)*VJ(3))
          ELSE
            DETJ=SQRT(VJ(2)*VJ(2)+VJ(4)*VJ(4))
          ENDIF
        ELSEIF (NDIM.EQ.3) THEN
          IF(IFL.LT.4) THEN
            DETJ=(VJ(5)*VJ(9)-VJ(8)*VJ(6))**2
            DETJ=DETJ+(VJ(8)*VJ(3)-VJ(2)*VJ(9))**2
            DETJ=DETJ+(VJ(2)*VJ(6)-VJ(5)*VJ(3))**2
            DETJ=SQRT(DETJ)
C            COEF=SQRT(DETJ)*VCPG1(IG)
          ELSEIF(IFL.GT.3.AND.IFL.LT.6) THEN
            DETJ=(VJ(4)*VJ(9)-VJ(7)*VJ(6))**2
            DETJ=DETJ+(VJ(7)*VJ(3)-VJ(1)*VJ(9))**2
            DETJ=DETJ+(VJ(1)*VJ(6)-VJ(4)*VJ(3))**2
            DETJ=SQRT(DETJ)
C            COEF=SQRT(DETJ)*VCPG1(IG)
          ELSEIF(IFL.GT.5) THEN
            DETJ=(VJ(4)*VJ(8)-VJ(7)*VJ(5))**2
            DETJ=DETJ+(VJ(7)*VJ(2)-VJ(1)*VJ(8))**2
            DETJ=DETJ+(VJ(1)*VJ(5)-VJ(4)*VJ(2))**2
            DETJ=SQRT(DETJ)
          ENDIF
        ENDIF
        
C---------------COMPUTE THE SOURCE TERM IN THE GAUSS POINT
        IF (NDIM.EQ.2) THEN
          IF (IFL.EQ.2) THEN
            C=VNI(INI+7)*V(1)+VNI(INI+8)*V(2)+VNI(INI+1)*V(3)
          ELSEIF (IFL.EQ.3) THEN
            C=VNI(INI+3)*V(1)+VNI(INI+4)*V(2)+VNI(INI+5)*V(3)
          ELSEIF (IFL.EQ.4) THEN
            C=VNI(INI+1)*V(1)+VNI(INI+2)*V(2)+VNI(INI+3)*V(3)
          ELSEIF (IFL.EQ.5) THEN
            C=VNI(INI+5)*V(1)+VNI(INI+6)*V(2)+VNI(INI+7)*V(3)
          ENDIF
        ELSEIF (NDIM.EQ.3) THEN
          IF (IFL.EQ.2) THEN
            C=VNI(INI+1)*V(1)+VNI(INI+9)*V(2)+VNI(INI+13)*V(3)+
     1       VNI(INI+20)*V(4)+VNI(INI+19)*V(5)+VNI(INI+12)*V(6)+
     2       VNI(INI+7)*V(7)+VNI(INI+8)*V(8)
          ELSEIF (IFL.EQ.3) THEN
            C=VNI(INI+3)*V(1)+VNI(INI+4)*V(2)+VNI(INI+5)*V(3)+
     1       VNI(INI+11)*V(4)+VNI(INI+17)*V(5)+VNI(INI+16)*V(6)+
     2       VNI(INI+15)*V(7)+VNI(INI+10)*V(8)
          ELSEIF (IFL.EQ.4) THEN
            C=VNI(INI+1)*V(1)+VNI(INI+2)*V(2)+VNI(INI+3)*V(3)+
     1       VNI(INI+10)*V(4)+VNI(INI+15)*V(5)+VNI(INI+14)*V(6)+
     2       VNI(INI+13)*V(7)+VNI(INI+9)*V(8)
          ELSEIF (IFL.EQ.5) THEN
            C=VNI(INI+5)*V(1)+VNI(INI+6)*V(2)+VNI(INI+7)*V(3)+
     1       VNI(INI+12)*V(4)+VNI(INI+19)*V(5)+VNI(INI+18)*V(6)+
     2       VNI(INI+17)*V(7)+VNI(INI+11)*V(8)
          ELSEIF (IFL.EQ.6) THEN
            C=VNI(INI+1)*V(1)+VNI(INI+8)*V(2)+VNI(INI+7)*V(3)+
     1       VNI(INI+6)*V(4)+VNI(INI+5)*V(5)+VNI(INI+4)*V(6)+
     2       VNI(INI+3)*V(7)+VNI(INI+2)*V(8)
          ELSEIF (IFL.EQ.7) THEN
            C=VNI(INI+13)*V(1)+VNI(INI+14)*V(2)+VNI(INI+15)*V(3)+
     1       VNI(INI+16)*V(4)+VNI(INI+17)*V(5)+VNI(INI+18)*V(6)+
     2       VNI(INI+19)*V(7)+VNI(INI+20)*V(8)
          ENDIF
        ENDIF
        
C............. compute the weight
        C=VCPG1(IG)*DETJ*C
        
                
C.................if axi-symmetric problem
        IF(IAX.EQ.1) THEN
          Y=ZERO
          DO  J=1,INEL
            Y=Y+VCORE(2*J)*VNI(INI+J)
          END DO
          C=C*2.*PI*Y
        END IF
        
        

C............. terms of the source matrix, frame only
        DO J=1,inel
          I2=INI+J
          IF (NDIM.EQ.2) THEN
            VSRE(4*(J-1)+1)=VSRE(4*(J-1)+1)+
     1           VNI(I2)*C*REAL(VNDSRC((IFREQ-1)*4+4))
            VSRE(4*(J-1)+2)=VSRE(4*(J-1)+2)+
     1           VNI(I2)*C*DIMAG(VNDSRC((IFREQ-1)*4+4))
          ELSEIF (NDIM.EQ.3) THEN
            VSRE(6*(J-1)+1)=VSRE(6*(J-1)+1)+
     1           VNI(I2)*C*REAL(VNDSRC((IFREQ-1)*10+9))
            VSRE(6*(J-1)+2)=VSRE(6*(J-1)+2)+
     1           VNI(I2)*C*DIMAG(VNDSRC((IFREQ-1)*10+9))
            VSRE(6*(J-1)+3)=VSRE(6*(J-1)+3)+
     1           VNI(I2)*C*REAL(VNDSRC((IFREQ-1)*10+10))
          ENDIF
        ENDDO
1030    INI=INI+(NDIM+1)*INEL
      
C........COPY SOURCE TERMS FOR frame displacement too relative displacement

      DO J=1,inel
        IF (NDIM.EQ.2) THEN
          VSRE(4*(J-1)+3)=VSRE(4*(J-1)+1)
          VSRE(4*(J-1)+4)=VSRE(4*(J-1)+2)
        ELSEIF (NDIM.EQ.3) THEN
          VSRE(6*(J-1)+4)=VSRE(6*(J-1)+1)
          VSRE(6*(J-1)+5)=VSRE(6*(J-1)+2)
          VSRE(6*(J-1)+6)=VSRE(6*(J-1)+3)
        ENDIF
      ENDDO      

C................... tilbakestill IGP
      IPG=IPG1

      RETURN


C----------EVALUTE THE LOSS MATRIX
1100  DO I=1,IKE
        VIE(I)=CZERO
      ENDDO

C------------EVALUATE THE FUNCTION Gj FOR VISCOUS LOSSES
      CALL GJ05(VPREE,1,GJ)

C------------LOOP OVER THE G.P
      IDIM1=NDIM-1
      IDECL=(NDIM+1)*INEL
      I1=1+INEL
      I2=0
      DO IG=1,IPG
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)

C--------------COMMON FACTOR: jw * gauss veight * jacobi detr.
C-----------If axisymmetry; find radius
        IF (IAX.NE.0) THEN
          R=0.0
          DO I=1,INEL
            R=R+VNI(IDECL*(IG-1)+I)*VCORE((I-1)*NDIM+2)
          ENDDO
        ENDIF

        D=DCMPLX(0.0,W)*VCPG(IG)*DETJ

        IF (IAX.NE.0) THEN
          D=D*R*2*PI
        ENDIF  

C------------ACCUMULATE LOSS TERMS
        IDL=0

C------------LOOP THROUGH ALL "SMALL MATRICES", COLOUMNS
        DO J=1,INEL
          JJ=I2+J
          IF (NDIM.EQ.2) THEN
            J0=6+0.5*IDL*(IDL+5)
          ELSEIF (NDIM.EQ.3) THEN
            J0=4+0.5*(IDL+3)*(IDL+4)
          ENDIF

C------------LOOP THROUGH ALL "SMALL MATRICES", LINES ON AND ABOVE THE DIAGONAL
          DO I=1,J

            II=I2+I
            C=VNI(II)*VNI(JJ)*D
C--------------Point at upper left corner of Cww and assemble Cww.
            VIE(J0)=VIE(J0)+C*VPREE(13)*Gj
            IF (NDIM.NE.1) THEN
              J1=J0+IDL+NDIM+2
              DO II=1,IDIM1
                VIE(J1)=VIE(J1)+C*VPREE(13)*Gj
                J1=J1+IDL+II+NDIM+2
C                J1=J1+J1+1
              ENDDO
            ENDIF
            J0=J0+2*NDIM
          
          ENDDO
          IDL=IDL+NDIM*2
          
        ENDDO
        I1=I1+IDECL
        I2=I2+IDECL
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

      END

