      SUBROUTINE ELEM09(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1      VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C**********************************************************************C
C ELEM09.FOR                                                           C
C  23/5-1995     Kompaniet Numerisk lyd                                C
C  Bruno Brouard                                                       C
C**********************************************************************C
C  Quadratic element for perforated plate 2 or 3                       C
C  dimensions:                                                         C
C     2 dimension: 6 nodes isoparametric element                       C
C     3 dimension: 16 nodes isoparametric element                      C
C  Number of integration points: 3 **(ndim -1)                         C
C  Number of degrees of freedom per node: 1                            C
C  Element matrix or vector formed by this subbprogram                 C
C  according to icode value :                                          C
C     ICODE.EQ.1  return of parameters                                 C
C     ICODE.EQ.2  evaluate interpolation functions and                 C
C                 numerical integration coeffisients                   C
C     ICODE.EQ.3  element matrix (VKE)  .....not written.....NB!!      C
C     ICODE.EQ.4  tangent matrix (VKE)  .....not written.....NB!!      C
C     ICODE.EQ.5  mass matrix (VKE)     .....not written.....NB!!      C
C     ICODE.EQ.6  K.U. product (VFE)    .....not written.....NB!!      C
C     ICODE.EQ.7  element load (VFE)  .....not written.....NB!!        C
C     ICODE.EQ.8  print gradients      .....not written.....NB!!       C
C     ICODE.EQ.9  impedance matrix(VKE)  .....not written.....NB!!     C
C     ICODE.EQ.10 source matrix (VKE)     .....not written.....NB!!    C
C  Element properties                                                  C
C     VPREE(1) density of the fluid                                    C
C     VPREE(2) viscosity of the fluid                                  C
C     VPREE(3) open area ratio of the plate                            C
C     VPREE(4) radius of the holes of the plates                       C
C     VPREE(5) atmospheric pressure                                    C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

c      COMMON/COOR/NDIM,NNT,NDLN,NDLT,FAC(3),IAX
c      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
c     &            IPG,ICODE,IDLE0,INEL0,IPG0
c
c      COMMON/ES/M,MR,MP
c
c      INTEGER IFL,IGIMP,IFREQ,IGSRC,IFL1,NI100NOD,NI100MOD,NI100FL
c      REAL*8 I100RAD,W,W2
c      COMMON/IMPEL/IFL,IGIMP,IFREQ,IGSRC,IFL1,NI100NOD,NI100MOD,
c     1		   NI100FL,I100RAD,W,W2
c
      include 'femak.cmn'
c      
      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VKE1(1),VFE(1)
      COMPLEX*16 VIE(1),VIE1(1),VSRE(1),VNDSRC(1),VPRSG(1),VDLEC(1)
      INTEGER KTYPRSG(1)
c      
C-----------CHARACTERISTIC DIMENSIONS OF THE ELEMENT
C----------(VALID  UP TO 3 DIMENSIONS)
C     DIMENSION VCPG(IPG),VKPG(NDIM*IPG)
      DIMENSION VCPG(9),VKPG(27)
C     DIMENSION VJ(NDIM*NDIM),VJ1(NDIM*NDIM)
      DIMENSION VJ(9),VJ1(9)
C     DIMENSION VNIX(INEL*NDIM),VNI((1+NDIM)*INEL*IPG),IPGKED(NDIM)
      DIMENSION VNIX(48),VNI(576),IPGKED(3)
      DIMENSION VCPG1(9),VKPG1(27)
c      
      DIMENSION FOK(8)
      INTEGER I
      REAL*8 KSI,PSI
      COMPLEX*16 FAK,BESS0,BESS1,CARRE,X,CARREP
      REAL*8  NR(3),GKSI(3),GETA(3)
c      
C------------Number of Gauss points in ksi,eta,dzeta direction
      DATA IPGKED/1,3,3/
C---------------------------------------------------
      DATA ZERO/0.0D0/,EPS/1.D-6/,CZERO/(0.0,0.0)/
      DATA PI/3.141592654/
      DATA   FOK /-1.4092, 0.0, 0.33818, 0.0, 0.06793, -0.02287,
     1             0.03015, -0.01641/
c
      IKE=IDLE*(IDLE+1)/2
c    
C-------------CHOOSE FUNCTION TO BE EXECUTED
      GO TO (100,200,300,400,500,600,700,800,900,1000),ICODE
c
C
C------------RETURN ELEMENT PARAMETERS IN COMMON /RGDT/
100   IF (NDIM.EQ.2) THEN
        IDLE0=6
        INEL0=6
        IPG0=3
      ELSEIF (NDIM.EQ.3) THEN
        IDLE0=16
        INEL0=16
        IPG0=9
      ENDIF 
      CFLAG=1
c
      RETURN
c
C
C-----------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C-----------THEIR DERIVATIVES AT G.P.
c
200   CONTINUE
      CALL GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
c
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
c
      CALL NI03(VKPG,VNI)
c
      IF (M.GE.2) THEN
        I1=(NDIM+1)*INEL*IPG
        WRITE(MP,2020) (VNI(I),I=1,I1)
      ENDIF
2020  FORMAT(/'FUNCTIONS N AND DERIVATIVES'/ (1X,8E12.5))
c   
      RETURN
c
c
C----------COMPUTE ELEMENT STIFFNESS MATRIX
c
c
C--------INITIALIZE VIE
300   DO 310 I=1,IKE
        VIE(I)=CZERO
310   CONTINUE
c
C--------LOOP OVER THE INTEGRATION POINTS
      INI=1+INEL
      DO 330 IG=1,IPG
C----------EVALUATE THE JACOBIAN MATRIX, ITS INVERSE AND ITS DETER,
        CALL JACOB(VNI(INI),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF (DETJ.LT.EPS) WRITE(MP,3030)IEL,IG,DETJ
3030    FORMAT('*** ELEM ',I5,'P.G. ',I3,' DETJ(J)=',E12.5,
     1       ' ,     ABS(DETJ(J)) USED')
        DETJ=ABS(DETJ)
c
C-----------If axisymmetry; find radius
        IF (IAX.NE.0) THEN
          R=0.0
          DO I=1,INEL
            R=R+VNI(3*INEL*(IG-1)+I)*VCORE(2*I)
          ENDDO
        ENDIF

C----------------CALCULATION OF THE MEAN THICNESS OF THE PLATE
C----------------FROM THE MIDDLE FACE NODES
        IF (NDIM.EQ.2) THEN
          THICK = SQRT( (VCORE(3)- VCORE(9))**2 +
     1                (VCORE(4)- VCORE(10))**2 )
        ELSEIF (NDIM.EQ.3) THEN
          THICK = SQRT( (VCORE(4)- VCORE(28))**2 +
     1                (VCORE(5)- VCORE(29))**2 +
     2                (VCORE(6)- VCORE(30))**2 )
        ENDIF 

C----------------PERFORM COEF
        COEF=VCPG(IG)*DETJ* 4/THICK/THICK
        IF (IAX.NE.0) THEN
          COEF = COEF*R*2*PI
        ENDIF
c
C------------ACCUMULATE TERMS OF THE ELEMENT MATRIX
        IK=0
        DO 320 J=1,IDLE
          DO 320 I=1,J
           IF (NDIM .EQ. 2) THEN
             C=VNI(INI+I-1)*VNI(INI+J-1)
           ELSE
             C=VNI(INI+2*INEL+I-1)*VNI(INI+2*INEL+J-1)
           ENDIF
           IK=IK+1
           VIE(IK)=VIE(IK)+CMPLX(C*COEF,0.D0)
320     CONTINUE
c
c
C----------NEXT G.P.
        INI=INI+(NDIM+1)*INEL
330   CONTINUE  
c
C CALCULATION OF THE MELLING TRANSFERT IMPEDANCE
C CALCUL DE LA FONCTION DE FOK
c      KSI = SQRT(VPREE(3) / PI)     essai de correction
      KSI = SQRT(VPREE(3))
c      
      PSI = 1.0
      DO I=1,8 
        PSI = PSI + FOK(I) * KSI**FLOAT(I)
      ENDDO
c
CC        X = SQRT(CMPLX(0.0,-1.0)) * RH * SQRT(W* RHO /MU)
      X = SQRT(CMPLX(.0,-1.))*VPREE(4)*SQRT(W*VPREE(1)/VPREE(2))
      CARREP = 1./ (1 - 2 * BESS1(X)/ BESS0(X) / X  )
      X = X / SQRT(2.)
      CARRE = 1./ (1 - 2 * BESS1(X)/ BESS0(X) / X  )
CC        Z = JJ*W*RHO*CARRE/PHI *(LP + 16*RH/3./PI*PSI)
      FAK = VPREE(1)/VPREE(3)*
     1     (THICK*CARRE+16*VPREE(4)*CARREP/3./PI*PSI)
c
      FAK = VPREE(1) * THICK / FAK
c                                   
      DO I=1,IKE
        VIE(I)=VIE(I)*FAK
      END DO
c        
      RETURN
c
C
C----------EVALUATE ELEMENT TANGENT MATRIX
400   CONTINUE
c
c
      RETURN
c
C
C----------MASS MATRIX
c
C--------INITIALIZE VKE
500   DO 510 I=1,IKE
        VKE(I)= ZERO
510   CONTINUE
c
      INI=0
      DO 530 IG=1,IPG
C----------EVALUATE THE JACOBIAN MATRIX
        I1=INI+INEL+1
        CALL JACOB(VNI(I1),VCORE,NDIM,INEL,VJ,VJ1,DETJ)
        IF (DETJ.LT.EPS) WRITE(MP,3030)IEL,IG,DETJ
        DETJ=ABS(DETJ)
c
C-----------If axisymmetry; find radius
        IF (IAX.NE.0) THEN
          R=0.0
          DO I=1,INEL
            R=R+VNI(INEL*(NDIM+1)*(IG-1)+I)*VCORE((I-1)*NDIM+2)
          ENDDO
        ENDIF
c
C----------------PERFORM COEF
        COEF=VCPG(IG)*DETJ
        IF (IAX.NE.0) THEN
          COEF = COEF*R*2*PI
        ENDIF
c
c
C------------TERMS OF THE MASS MATRIX
        IK=0
        DO 520 J=1,IDLE
          DO 520 I=1,J
            IK=IK+1
            I1=INI+I
            I2=INI+J
            VKE(IK)=VKE(IK)+VNI(I1)*VNI(I2)*COEF
520     CONTINUE
c
C----------NEXT G.P.
        INI=INI+(NDIM+1)*INEL
530   CONTINUE
c
      FAK = CMPLX(VPREE(1) * W * W / VPREE(5),0.D0)
      DO I=1,IKE
        VIE(I)=CMPLX(VKE(I),0.D0)*FAK
      END DO
c 
      RETURN
c
c
C    
C-------------EVALUATE THE ELEMENT RESIDUAL
600   CONTINUE
c
      RETURN
c
C
C------------EVALUATE FE
700   CONTINUE
c
      RETURN
c
C
C-------------EVALUATE AND PRINT GRADIENTS AT G.P.
800   CONTINUE
c
      RETURN
C
C ----------EVALUATE THE IMPEDANCE MATRIX
900   CONTINUE
c
      RETURN
c
C-----------EVALUATE THE SOURCE MATRIX
1000  CONTINUE
c     
      RETURN
      END
c
