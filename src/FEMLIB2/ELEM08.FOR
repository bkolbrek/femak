      SUBROUTINE ELEM08(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1      VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
     
C      SUBROUTINE ELEM08(VCORE,VPRNE,VPREE,VDLE,VKE,VFE)
C**********************************************************************C
C ELEM08.FOR                                                           C
C  23/1-1995  Kompaniet Numerisk lyd                                   C
C  Marc Dhainaut                                                       C
C**********************************************************************C
C  NDIM = 2 : BEAM ELEMENT   3 nodes                                   C
C  NDIM = 3 : PLATE ELEMENT  9 nodes                                   C
C  9 nodes quadratic element for 2 dimensional elasticity              C
C  8 nodes is possible : - change all the CALL NI08 by CALL NI02       C
C                        - change IDLE0=27 by 24 and INEL0=9 by 8      C
C----------------------------------------------------------------------C
C  Evaluate element information according to ICODE value               C
C     ICODE.EQ.1  return of parameters                                 C
C     ICODE.EQ.2  evaluate interpolation functions and                 C
C                 numerical integration coefficients                   C
C     ICODE.EQ.3  stiffnes matrix (VKE)                                C
C     ICODE.EQ.4  tangent matrix (VKE)  .....not written.....NB!!      C
C     ICODE.EQ.5  mass matrix (VKE)                                    C
C     ICODE.EQ.6  residuals             .....not written.....NB!!      C
C     ICODE.EQ.7  element load (VFE)    .....not written.....NB!!      C
C     ICODE.EQ.8  evaluate and print    .....not written.....NB!!      C
C----------------------------------------------------------------------C
C  Element properties                                                  C
C     VPREE(1) : Young's modulus                                       C
C     VPREE(2) : Poisson's coeffiecient                                C
C     VPREE(3).EQ.0 : Plane stress                                     C
C             .EQ.1 : Plane strain                                     C                     
C     VPREE(4) : specific mass                                         C
C     VPREE(5) : thickness                                             C
C----------------------------------------------------------------------C
C  Symmetry of the element                                             C
c     The matrices (Ke, Me) of the element are always symmetric        C
C**********************************************************************C
C  Programmert etter Dhatt & Touzot, s.194-197                         C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
c      COMMON/COOR/NDIM
c      COMMON/ASSE/NSYM
c      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
c     &            IPG,ICODE,IDLE0,INEL0,IPG0
c      COMMON/ES/M,MR,MP

      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VKE1(1),VFE(1)
      COMPLEX*16 VIE(1),VIE1(1),VNDSRC(1),VPRSG(1),VSRE(1),VDLEC(1)
      INTEGER KTYPRSG(1)
c      DIMENSION VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VFE(1)

C-----------CHARACTERISTIC DIMENSIONS OF THE ELEMENT
C----------(VALID  UP TO 3 DIMENSIONS)
C     DIMENSION VCPG(IPG),VKPG(NDIM*IPG),VDE1(IMATD**2)
      DIMENSION VCPG(9),VKPG(18),VDE1(9)
C     DIMENSION VBBE(IMATB*IDLE),VBSE(NDIM*IDLE),VDE(IMATD**2),VJ(NDIM*NDIM),VJ1(NDIM*NDIM)
      DIMENSION VBBE(81),VBSE(54),VDE(9),VJ(4),VJ1(4)
C     DIMENSION VNIX(INEL*NDIM),VNI((1+NDIM)*INEL*IPG),IPGKED(NDIM)
      DIMENSION VNIX(18),VNI(243),IPGKED(2)
C     DIMENSION VCORE1(3*INEL)
      DIMENSION VCORE1(27)
C
      SAVE VCPG,VKPG,VNI


      REAL*8 E,XNU,RHO,H,G,DB2,DB3,DS2,DS3,DM2,DM3,C
      INTEGER II,COL,ROW,PNTCO,PNTCO1,PNTRW1
C------------Dimension of D-matrix
      DATA IMATD/3/
C---------------------------------------------------
      DATA ZERO/0.0D0/,DEUX/2.D0/,X05/0.5D0/,RADN/.572957795130823D2/ 
      DATA N0/0/
      DATA EPS/1.D-6/
      SQRT(X)=DSQRT(X)
      ATAN2(X,Y)=DATAN2(X,Y)
 
      NDIM0=NDIM
      
C----------ELASTIC PARAMETERS                  
      E  =VPREE(1)
      XNU=VPREE(2)
      RHO=VPREE(4)
      H  =VPREE(5)
C     width of the beam : w=1
      G =E/2./(1.+XNU)
      DB2 = E*(H**3)/12.
      DB3 = H**3/12
      DS2 = G*H/1.2
      DS3 = G*H/1.2
      DM2 = RHO
      DM3 = RHO
      
C----------TEST      
C      WRITE(*,*) 'ELEM08 : ICODE =',ICODE
      
C-------------CHOOSE FUNCTION TO BE EXECUTED 

      GO TO (100,200,300,400,500,600,700,800),ICODE
 
C================================================================================
C------------RETURN ELEMENT PARAMETERS IN COMMON /RGDT/ 
C================================================================================
100   IF (NDIM.EQ.2) THEN
         IDLE0=6
         INEL0=3
         IPG0=3
      ENDIF
      IF (NDIM.EQ.3) THEN
         IDLE0=27
         INEL0=9
         IPG0=9 
Cserendipity elt IDLE0=24
Cserendipity elt INEL0=8
      ENDIF    
      
      CFLAG=0

      RETURN

C================================================================================
C-----------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C-----------THEIR DERIVATIVES AT G.P.
C================================================================================
200   IF (NDIM0.EQ.2) THEN
          NDIM=1
          IPGKED(1)=2
      ENDIF    
      IF (NDIM0.EQ.3) THEN
          NDIM=2
          IPGKED(1)=3
          IPGKED(2)=3
      ENDIF    
      CALL GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
      IF(M.GE.2) THEN
         WRITE(MP,2000)IPG
C         I0=1
C         DO IG=1,IPG
C            I1=I0+NDIM-1
C            WRITE(MP,2010)VCPG(IG),(VKPG(I),I=I0,I1)
C            I0=I0+NDIM
C         ENDDO
      ENDIF
         
      CALL NI08(VKPG,VNI) 
Cserendipity elt CALL NI02(VKPG,VNI)
      IF (M.GE.2) THEN
         I1=3*INEL*IPG
         WRITE(MP,2020) (VNI(I),I=1,I1)
      ENDIF
      
2000  FORMAT(/I5,'  GAUSS POINTS'/10X,'VCPG',25X,'VKPG')      
2010  FORMAT(1X,F20.15,5X,3F20.15)
2020  FORMAT(/'FUNCTIONS N AND DERIVATIVES'/ (1X,8E12.5))
      
      NDIM=NDIM0      
      RETURN
      
C      
C================================================================================
C----------COMPUTE ELEMENT STIFFNESS MATRIX
C----------  .THE BENDING CONTRIBUTION USES THE FULL GAUSS INTEGRATION (3x3)
C----------  .THE SHEAR CONTRIBUTION USES THE REDUCED GAUSS INTEGRATION (2x2)
C================================================================================
C------ INITIALIZE VKE
300   DO I=1,IDLE*(IDLE+1)/2
        VKE(I)=ZERO
      ENDDO

C-----------------------------------        
C------ BENDING STIFFNESS MATRIX
C-----------------------------------
      IF (NDIM0.EQ.2) THEN
          NDIM=1
          IPGKED(1)=2
      ENDIF    
      IF (NDIM0.EQ.3) THEN
          NDIM=2
          IPGKED(1)=3
          IPGKED(2)=3
      ENDIF

      CALL GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
      IF(M.GE.2) THEN
         WRITE(MP,2000)IPG
c         I0=1
c         DO IG=1,IPG
c            I1=I0+NDIM-1
c            WRITE(MP,2010)VCPG(IG),(VKPG(I),I=I0,I1)
c            I0=I0+NDIM
c         ENDDO
      ENDIF

      CALL NI08(VKPG,VNI)
Cserendipity elt CALL NI02(VKPG,VNI)
            
C-------CHANGE NDIM COORDINATES TO NDIM-1 COORDINATES
      CALL COOR23D(VCORE,VCORE1)
      
C-------form elastic matrix D of the plate
      IF (NDIM0.EQ.3) THEN
        CALL D02(VPREE,VDE)
        IF (M.GE.2) WRITE(MP,2033) (VDE(I),I=1,9)
      ENDIF
2033  FORMAT(/' MATRIX D'/1X,9E12.5)

C-------LOOP OVER G.P
      I1=1+INEL
      
      DO 333 IG=1,IPG
      
C----------------EVALUATE THE JACOBIAN, ITS INVERSE AND ITS DETERMINANT      
      CALL JACOB(VNI(I1),VCORE1,NDIM,INEL,VJ,VJ1,DETJ)
      IF (DETJ.LT.EPS) WRITE(MP,2040)IEL,IG,DETJ
      IF (M.GE.2) WRITE(MP,2050)VJ,VJ1,DETJ
2040  FORMAT(' *** ELEM ',I5,' G.P. ',I3,' DET(J)=',E12.5)      
2050  FORMAT(/' JACOBIAN=',4E12.5 / ' J INVERS=',4E12.5/' DETJ=',E12.5)
        
C----------------PERFORM D*COEF
      IF (NDIM0.EQ.2) C=VCPG(IG)*DB2*DETJ
      IF (NDIM0.EQ.3) THEN
          C=VCPG(IG)*DB3*DETJ
          DO I=1,9
             VDE1(I)=VDE(I)*C
          ENDDO
      ENDIF
             
C----------------FORM DERIVATIVES OF Ni
      CALL DNIDX(VNI(I1),VJ1,NDIM,INEL,VNIX)
      IF (M.GE.2) WRITE(MP,2063)(VNIX(I),I=1,INEL*NDIM)      
2063  FORMAT(/' VNIX'/(1X,9E12.5))

C----------------FORM MATRIX VKE FOR BENDING
      IF (NDIM0.EQ.2) THEN
         VKE( 3)=VKE( 3) + VNIX(1)*VNIX(1)*C
         VKE( 8)=VKE( 8) + VNIX(1)*VNIX(2)*C
         VKE(10)=VKE(10) + VNIX(2)*VNIX(2)*C
         VKE(17)=VKE(17) + VNIX(1)*VNIX(3)*C
         VKE(19)=VKE(19) + VNIX(2)*VNIX(3)*C
         VKE(21)=VKE(21) + VNIX(3)*VNIX(3)*C
      ENDIF
      IF (NDIM0.EQ.3) THEN
          CALL B08(VNIX,INEL,VBBE)
          IF (M.GE.2) WRITE(MP,2073)(VBBE(I),I=1,81)
2073  FORMAT(/' MATRIX B bending'/(1X,10E12.5))
          CALL BTDB(VKE,VBBE,VDE1,IDLE,IMATD,N0)
      ENDIF
          
333   I1=I1+NDIM0*INEL
        
C-----------------------------------        
C------ SHEAR STIFFNESS MATRIX
C------(with reduced integration)
C------ added to bending stiffness matrix
C-----------------------------------        
C
C----------------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C----------------THEIR DERIVATIVES AT G.P.
C---------------- FOR reduced integration in the case of the plate !!!
      IF (NDIM0.EQ.2) THEN
          NDIM=1
          IPGKED(1)=2
      ENDIF    
      IF (NDIM0.EQ.3) THEN
          NDIM=2
          IPGKED(1)=2
          IPGKED(2)=2
      ENDIF
      CALL GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
      IF(M.GE.2) THEN
         WRITE(MP,2000)IPG
c         I0=1
c         DO IG=1,IPG
c            I1=I0+NDIM-1
c            WRITE(MP,2010)VCPG(IG),(VKPG(I),I=I0,I1)
c            I0=I0+NDIM
c         ENDDO
      ENDIF

      CALL NI08(VKPG,VNI)
Cserendipity CALL NI02(VKPG,VNI)
      
C2000  FORMAT(/I5,'  GAUSS POINTS'/10X,'VCPG',25X,'VKPG')      
C2010  FORMAT(1X,F20.15,5X,3F20.15)    
C----------------LOOP OVER G.P
      I1=1+INEL 
      I2=1
      DO 334 IG=1,IPG
      
C----------------EVALUATE THE JACOBIAN, ITS INVERSE AND ITS DETERMINANT
      CALL JACOB(VNI(I1),VCORE1,NDIM,INEL,VJ,VJ1,DETJ)
      IF (DETJ.LT.EPS) WRITE(MP,2040)IEL,IG,DETJ
      IF (M.GE.2) WRITE(MP,2050)VJ,VJ1,DETJ
      
C----------------PERFORM D*COEF
      IF (NDIM0.EQ.2) C=VCPG(IG)*DS2*DETJ
      IF (NDIM0.EQ.3) C=VCPG(IG)*DS3*DETJ
            
C----------------FORM DERIVATIVES OF Ni
      CALL DNIDX(VNI(I1),VJ1,NDIM,INEL,VNIX)
      IF (M.GE.2) WRITE(MP,2064)(VNIX(I),I=1,NDIM*INEL)
2064  FORMAT(/' VNIX'/(1X,9E12.5))

C----------------FORM MATRIX VKE
      IF (NDIM0.EQ.2) THEN
          VKE( 1)=VKE( 1) + C*VNIX(1)*VNIX(1)
          VKE( 2)=VKE( 2) - C*VNIX(1)*VNI(I2)
          VKE( 3)=VKE( 3) + C*VNI(I2)*VNI(I2)
          VKE( 4)=VKE( 4) + C*VNIX(1)*VNIX(2)
          VKE( 5)=VKE( 5) - C*VNIX(2)*VNI(I2)
          VKE( 6)=VKE( 6) + C*VNIX(2)*VNIX(2)
          VKE( 7)=VKE( 7) - C*VNIX(1)*VNI(I2+1)
          VKE( 8)=VKE( 8) + C*VNI(I2)*VNI(I2+1)
          VKE( 9)=VKE( 9) - C*VNIX(2)*VNI(I2+1)
          VKE(10)=VKE(10) + C*VNI(I2+1)*VNI(I2+1)
          VKE(11)=VKE(11) + C*VNIX(1)*VNIX(3)
          VKE(12)=VKE(12) - C*VNIX(3)*VNI(I2)          
          VKE(13)=VKE(13) + C*VNIX(2)*VNIX(3)
          VKE(14)=VKE(14) - C*VNIX(3)*VNI(I2+1)
          VKE(15)=VKE(15) + C*VNIX(3)*VNIX(3)
          VKE(16)=VKE(16) - C*VNIX(1)*VNI(I2+2)
          VKE(17)=VKE(17) + C*VNI(I2)*VNI(I2+2)
          VKE(18)=VKE(18) - C*VNIX(2)*VNI(I2+2)
          VKE(19)=VKE(19) + C*VNI(I2+1)*VNI(I2+2)
          VKE(20)=VKE(20) - C*VNIX(3)*VNI(I2+2)
          VKE(21)=VKE(21) + C*VNI(I2+2)*VNI(I2+2)
      ENDIF
      IF (NDIM0.EQ.3) THEN
          CALL BB08(VNI(I2),VNIX,INEL,VBSE)
          IF (M.GE.2) WRITE(MP,2074)(VBSE(I),I=1,54)
2074  FORMAT(/' MATRIX B shear'/(1X,10E12.5))
          CALL BTB08(VKE,VBSE,C,IDLE,NDIM,N0)
      ENDIF
      I1=I1+NDIM0*INEL
334   I2=I2+NDIM0*INEL


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

      NDIM=NDIM0
      RETURN 
      
C================================================================================
C----------COMPUTE ELEMENT MASS MATRIX
C================================================================================
C-------INITIALIZE VKE
500   DO I=1,IDLE*(IDLE+1)/2
        VKE(I)=ZERO
      ENDDO

      IF (NDIM0.EQ.2) THEN
          NDIM=1
          IPGKED(1)=3
      ENDIF    
      IF (NDIM0.EQ.3) THEN
          NDIM=2
          IPGKED(1)=3
          IPGKED(2)=3
      ENDIF

C-------CHANGE NDIM COORDINATES TO NDIM-1 COORDINATES
      CALL COOR23D(VCORE,VCORE1)
      
C-------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C-------THEIR DERIVATIVES AT G.P.
      CALL GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
      IF (M.GE.2) THEN
         WRITE(MP,2000)IPG
C         I0=1
C         DO IG=1,IPG
C            I1=I0+NDIM-1
C            WRITE(MP,2010)VCPG(IG),(VKPG(I),I=I0,I1)
C            I0=I0+NDIM
C         ENDDO
      ENDIF

      CALL NI08(VKPG,VNI)
Cserendipity CALL NI02(VKPG,VNI)

C2000  FORMAT(/I5,'  GAUSS POINTS'/10X,'VCPG',25X,'VKPG')      
C2010  FORMAT(1X,F20.15,5X,3F20.15)
     
C-------LOOP OVER THE G.P
      I1=1+INEL
      I2=0
      
      DO 553 IG=1,IPG
      
C----------------EVALUATE THE JACOBIAN, ITS INVERSE AND ITS DETERMINANT      
      CALL JACOB(VNI(I1),VCORE1,NDIM,INEL,VJ,VJ1,DETJ)
      IF (DETJ.LT.EPS) WRITE(MP,2040)IEL,IG,DETJ
      IF (M.GE.2) WRITE(MP,2050)VJ,VJ1,DETJ      
C2040  FORMAT(' *** ELEM ',I5,' G.P. ',I3,' DET(J)=',E12.5)      
C2050  FORMAT(/' JACOBIAN=',4E12.5 / ' J INVERS=',4E12.5/' DETJ=',E12.5)
        
C----------------PERFORM D*COEF
      IF (NDIM0.EQ.2) C=VCPG(IG)*DETJ*DM2
      IF (NDIM0.EQ.3) C=VCPG(IG)*DETJ*DM3
      
C----------------ACCUMULATE MASS TERMS
      IF (NDIM0.EQ.2) THEN
          VKE( 1)=VKE( 1) + VNI(I2+1)*VNI(I2+1)*C*H
          VKE( 3)=VKE( 3) + VNI(I2+1)*VNI(I2+1)*C*(H**3/12)
          VKE( 4)=VKE( 4) + VNI(I2+1)*VNI(I2+2)*C*H
          VKE( 6)=VKE( 6) + VNI(I2+2)*VNI(I2+2)*C*H
          VKE( 8)=VKE( 8) + VNI(I2+1)*VNI(I2+2)*C*(H**3/12)
          VKE(10)=VKE(10) + VNI(I2+2)*VNI(I2+2)*C*(H**3/12)
          VKE(11)=VKE(11) + VNI(I2+1)*VNI(I2+3)*C*(H)
          VKE(13)=VKE(13) + VNI(I2+2)*VNI(I2+3)*C*(H)
          VKE(15)=VKE(15) + VNI(I2+3)*VNI(I2+3)*C*(H)
          VKE(17)=VKE(17) + VNI(I2+1)*VNI(I2+3)*C*(H**3/12)
          VKE(19)=VKE(19) + VNI(I2+2)*VNI(I2+3)*C*(H**3/12)
          VKE(21)=VKE(21) + VNI(I2+3)*VNI(I2+3)*C*(H**3/12)
      ENDIF
      IF (NDIM0.EQ.3) THEN    
          IDL=0
          DO J=1,INEL
             JJ=I2+J
             J0=1+IDL*(IDL+1)/2
             DO I=1,J
                II=I2+I
                D=VNI(II)*VNI(JJ)*C
                VKE(J0)=VKE(J0)+D*(H)
                J1=J0+IDL+2
                VKE(J1)=VKE(J1)+D*(H**3/12)
                J1=J1+IDL+3
                VKE(J1)=VKE(J1)+D*(H**3/12)
                J0=J0+NDIM+1
             ENDDO
             IDL=IDL+NDIM+1
          ENDDO
      ENDIF
      
      I1=I1+NDIM0*INEL
553   I2=I2+NDIM0*INEL


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
      
      NDIM=NDIM0
      RETURN
            
C********************************************************************************************
C********************************************************************************************      
                         
C----------EVALUATE ELEMENT TANGENT MATRIX
400   NDIM=NDIM0
      RETURN
C-------------EVALUATE THE ELEMENT RESIDUAL
600   NDIM=NDIM0
      RETURN
C------------EVALUATE VOLUMIC FORCES, FX, FY PER UNIT VOLUME
700   NDIM=NDIM0
      RETURN
C-------------EVALUATE AND PRINT STRESSES AT G.P.
800   NDIM=NDIM0
      RETURN
      
      END
