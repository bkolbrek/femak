      SUBROUTINE ELEM95(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1      VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C********************************************************************C
C ELEM95.FOR                                                         C
C  16/3-1995     NTH                                                 C
C  Marc Dhainaut & Tonni Johansen                                    C
C********************************************************************C
C  Coupling element: plate/beam - porous material using Biot-theory  C                                                       C
C  Direction 1 and 2 of the porous material must coincide with       C
C  direction x and y of the plate. Only the third porous direction   C
C  can be + or - the direction z of the plate                        C
C                                                                    C
C  the method used is based on penalty function                      C
C--------------------------------------------------------------------C
C  NDIM = 2 : BEAM (3 nodes) + POROUS (3 nodes) = 6 NODES (18 DOF)   C
C  NDIM = 3 : PLATE (9 nodes) + POROUS (8 nodes) = 17 NODES (75 DOF) C
C--------------------------------------------------------------------C
C  Evaluate element information according to ICODE value             C
C     ICODE.EQ.1  return of parameters                               C
C     ICODE.EQ.2  evaluate interpolation functions and               C
C                 numerical integration coeffisients                 C
C     ICODE.EQ.3  stiffness matrix (VKE)                             C
C     ICODE.EQ.5  not used                                           C
C     ICODE.EQ.6  not used                                           C
C     ICODE.EQ.8  not used                                           C
C--------------------------------------------------------------------C
C  Element properties                                                C
C     VPREE(1) Young's modulus                                       C
C     VPREE(2) Poisson's coeffiesient                                C
C     VPREE(3) Plate thickness                                       C
C     VPREE(4) = 1 normal of beam/plate towards the  x direction     C
C              =-1 normal of beam/plate towards the -x direction     C
C              = 2 normal of beam/plate towards the  y direction     C
C              =-2 normal of beam/plate towards the -y direction     C
C              = 3 normal of beam/plate towards the  z direction     C
C              =-3 normal of beam/plate towards the -z direction     C
C     VPREE(5) = 1 or 2 or 3
C               positive bending angle 1 of the plate induce positiveC
C               displacement of the interface nodes in x or y or z   C
C              =-1 or -2 or -3                                       C
C               positive bending angle 1 of the plate induce negativeC
C               displacement of the interface nodes in x or y or z   C
C     VPREE(6) = 1 or 2 or 3                                         C
C               positive bending angle 2 of the plate induce positiveC
C               displacement of the interface nodes in x or y or z   C
C              =-1 or -2 or -3                                       C
C               positive bending angle 2 of the plate induce negativeC
C               displacement of the interface nodes in x or y or z   C     
C--------------------------------------------------------------------C
C  Symmetry of the element                                           C
c     The matrix Ke of the element is always symmetric               C
C********************************************************************C
C Theory and implementation of this element is described in FEMAK    C
C programmer's manual and user's manual                              C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
C $DECLARE      
      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VKE1(1),VFE(1)
      COMPLEX*16 VIE(1),VIE1(1),VNDSRC(1),VPRSG(1),VSRE(1),VDLEC(1)
      INTEGER KTYPRSG(1)
      

C-----------CHARACTERISTIC DIMENSIONS OF THE ELEMENT
C----------(VALID  UP TO 3 DIMENSIONS)
C     DIMENSION VCPG(IPG),VKPG(NDIM*IPG)
      REAL*8 VCPG(27),VKPG(81)   
C     DIMENSION VJ(NDIM*NDIM),VJ1(NDIM*NDIM)
      REAL*8 VJ(9),VJ1(9)
C     DIMENSION VNI((1+NDIM)*INEL*IPG),IPGKED(NDIM)
      REAL*8 VNI(216),VNIS(243)
      REAL*8 COEF(3)
      INTEGER IPGKED(3)

      SAVE VCPG,VKPG,VNI,VNIS


      INTEGER IDL,J0,LENR,IDLD,J1,I,I0,IG,I1,INI,IKE
      INTEGER COL,ROW,PNTCO,PNTCO1,PNTRW1
      REAL*8 DETJ,EPS,ZERO,R,C,ALPHA
      REAL*8 NR(3),NRLEN,GKSI(3),GETA(3)
      CHARACTER*2 AXE         
      REAL*8 H,H1,H2
      INTEGER IAXE,INORM,IT1,IT2
C------------Number of Gauss points in ksi,eta,dzeta direction
      DATA IPGKED/3,3,3/
C---------------------------------------------------
      DATA EPS/1.D-6/,ZERO/0.0/
     
      NDIM0=NDIM
      INEL0=INEL     
C
      H     = VPREE(3)
      IAXE  = VPREE(4)
      IT1   = VPREE(5)
      IT2   = VPREE(6)
      IF (ABS(IAXE).EQ.1) AXE='XX'
      IF (ABS(IAXE).EQ.2) AXE='YY'
      IF (ABS(IAXE).EQ.3) AXE='ZZ'
C   problem of sign at the interface :      
      INORM=IAXE/ABS(IAXE)
      H1   =(IT1/ABS(IT1))*H
      IF (NDIM.EQ.3) H2=(IT2/ABS(IT2))*H      

      IF(NSYM.EQ.0) THEN
        IKE=IDLE*(IDLE+1)/2
      ELSE
        IKE=IDLE*IDLE
      ENDIF

C-------------CHOOSE FUNCTION TO BE EXECUTED
      GO TO (100,200,300,400,500,600,700,800,900,1000,1100),ICODE

C
C================================================================================
C------------RETURN ELEMENT PARAMETERS IN COMMON /RGDT/
C================================================================================
100   IF (NDIM.EQ.2) THEN
        IDLE0=18
        INEL0=6
        IPG0=9
      ELSEIF (NDIM.EQ.3) THEN
        IDLE0=75
        INEL0=17
        IPG0=27
      ENDIF

      CFLAG=0
      
      RETURN

C
C================================================================================
C-----------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C-----------THEIR DERIVATIVES AT G.P.
C================================================================================
200   IF (NDIM0.EQ.2) THEN
        NDIM=1
        INEL=3
      ELSEIF (NDIM0.EQ.3) THEN
        NDIM=2
        INEL=8
      ENDIF
      
      CALL GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
              
      IF(M.GE.2) THEN
        WRITE(MP,2000)IPG
2000  FORMAT(/I5,'  GAUSS POINTS'/10X,'VCPG',25X,'VKPG')
        I0=1
        DO IG=1,IPG
          I1=I0+NDIM-1
          WRITE(MP,2010)VCPG(IG),(VKPG(I),I=I0,I1)
2010  FORMAT(1X,F20.15,5X,3F20.15)          
          I0=I0+NDIM
        ENDDO
      ENDIF

      CALL NI01(VKPG,VNI)
      IF (M.GE.2) THEN
        I1=3*INEL*IPG
        WRITE(MP,2020) (VNI(I),I=1,I1)
2020  FORMAT(/'FUNCTIONS N AND DERIVATIVES'/ (1X,8E12.5))        
      ENDIF
      
      IF (NDIM0.EQ.2) THEN
        DO I=1,(NDIM+1)*INEL*IPG
           VNIS(I)=VNI(I)
        ENDDO   
      ELSEIF (NDIM0.EQ.3) THEN
        INEL=9
        CALL NI08(VKPG,VNIS)
      ENDIF
      
      NDIM=NDIM0
      INEL=INEL0
      RETURN
 
C================================================================================
C----------COMPUTE ELEMENT STIFFNESS MATRIX
C================================================================================
C--------INITIALIZE VKE
300   DO I=1,IKE
        VKE(I)=ZERO
      ENDDO

C........The penalty number ALPHA
      ALPHA=VPREE(1)*1.D7

C..........Set to tricky number of nodes to get correct shape functions and Gauss points.
      IF (NDIM0.EQ.2) THEN
        INEL=3
        INELS=3
      ELSEIF (NDIM0.EQ.3) THEN
        INEL=8
        INELS=9
      ENDIF

C..........Loop over Gauss Points.
C...........if NDIM=3 : The determinant and therefore the normal to the plate will
C...........be calculated from the 8 porous nodes (we use VNI and INEL=8, not VNI9 and INEL=9)
      INI=1+INEL
      IN1=0
      IN2=0
            
      DO 310 IG=1,IPG
      
C----------------EVALUATE THE JACOBIAN (ALSO ITS INVERSE AND ITS DETERMINANT)
      CALL JACOB(VNI(INI),VCORE,NDIM,INEL,VJ,VJ1,DETJ)


c.........Find Jacobi determinant and normal to element

      IF (NDIM0.EQ.2) THEN        
C.........2D-case
C.........Calculate DETJ and the normal to the line
          DETJ=SQRT(VJ(3)**2+VJ(1)**2)
          IF (DETJ.LT.EPS) WRITE(MP,2042)IEL,IG,DETJ
2042  FORMAT(' *** ELEM ',I5,' G.P. ',I3,' DET(J)=',E12.5)

          NR(1)=VJ(3)/DETJ
          NR(2)=-VJ(1)/DETJ    
          IF (ABS(NR(1)).LT.EPS) THEN
            NR(1)=ZERO
          ENDIF
          IF (ABS(NR(2)).LT.EPS) THEN
            NR(2)=ZERO
          ENDIF

      ELSEIF (NDIM0.EQ.3) THEN
C.........3D-case
C.........Calculate DETJ and the normal to the surface, according to the
C.........porous cordinates (i.e : VNI) not to the plate coordinates (i.e : VNI9)
          DETJ=SQRT((VJ(4)*VJ(8)-VJ(7)*VJ(5))**2
     1             +(VJ(7)*VJ(2)-VJ(1)*VJ(8))**2
     2             +(VJ(1)*VJ(5)-VJ(4)*VJ(2))**2)
          IF (DETJ.LT.EPS) WRITE(MP,2043)IEL,IG,DETJ
2043  FORMAT(' *** ELEM ',I5,' G.P. ',I3,' DET(J)=',E12.5)
     
          NR(1) = (VJ(4)*VJ(8)-VJ(7)*VJ(5))/DETJ
          NR(2) = (VJ(7)*VJ(2)-VJ(1)*VJ(8))/DETJ
          NR(3) = (VJ(1)*VJ(5)-VJ(4)*VJ(2))/DETJ              
          IF (ABS(NR(1)).LT.EPS) THEN
            NR(1)=ZERO
          ENDIF
          IF (ABS(NR(2)).LT.EPS) THEN
            NR(2)=ZERO
          ENDIF
          IF (ABS(NR(3)).LT.EPS) THEN
            NR(3)=ZERO
          ENDIF

      ENDIF

C........Write results for Jacobian and normal
      IF (M.GE.2) THEN
        IF(NDIM.EQ.2) THEN
          WRITE(MP,2050)VJ,VJ1,DETJ
          WRITE(MP,2055)NR
        ELSEIF (NDIM.EQ.3) THEN
          WRITE(MP,2051)VJ,VJ1,DETJ
          WRITE(MP,2056)NR
        ENDIF
      ENDIF

2050  FORMAT(/' JACOBIAN=',4E12.5 / ' J INVERS=',
     1         4E12.5/' DETJ=',E12.5)
2051  FORMAT(/' JACOBIAN=',9E12.5 / ' J INVERS=',
     1         9E12.5/' DETJ=',E12.5)
2055  FORMAT(/' NORMAL=',E12.5,'i+',E12.5,'j')
2056  FORMAT(/' NORMAL=',E12.5,'i+',E12.5,'j+',E12.5,'k')


C------- penalty number * Gauss weight * det(J)
      C = ALPHA*VCPG(IG)*DETJ
      

C----------------FORM MATRIX VKE
C................the filling in of VKE procedes in that way :
C.......II = 1,INEL represents the lignes of the block matrix
C..........JJ =II,INEL reprensents the columns of the block matrix
C..............each block has size IDLPxIDLP (=4x4 for 2D, 6x6 for 3D)
C..............those blocks represent the cross product of porous dof.
C..........JJ = 1,INELS reprensents the columns of the block matrix
C..............each block has size IDLPxIDLS (=4x2 for 2D, 4x3 for 3D) 
C..............those blocks represent the coupling porous dof / plate dof
C.......II = 1,INELS represents the last lines of the block matrix
C..........JJ =II,INELS represents the last columns of the block matrix
C..............each block has size IDLPxIDLP
C..............those blocks represent the cross product of plate dof.
C.. IDLP = dof for a porous node
C...IDLS = dof for a structural node

      IF (NDIM0.EQ.2) THEN
          IDLP=4
          IDLS=2
      ELSEIF (NDIM0.EQ.3) THEN
          IDLP=6
          IDLS=3
      ENDIF      
      
      i=0
            
C        porous/porous in 2 or 3 dimension
      DO JJ=1,INEL             
          DO indic = 1,IDLP/2
             idl=i*(i+1)/2 + indic
             DO II=1,JJ
               VKE(idl) = VKE(idl)+ (VNI(IN1+II)*VNI(IN1+JJ))*C
               idl=idl+IDLP
             ENDDO
             i=i+1
          ENDDO
          DO indic = IDLP/2+1,IDLP
             idl=i*(i+1)/2 + indic
             DO II=1,JJ
               in=indic-IDLP/2
               VKE(idl)=VKE(idl)+
     1                 (VNI(IN1+II)*NR(in)*VNI(IN1+JJ)*NR(in))*C
               idl=idl+IDLP
             ENDDO
             i=i+1
          ENDDO
      ENDDO
      
C        porous/plate (ii=1,inel)  and plate/plate (ii=1,jj)
c................
c     NDIM = 2   
C................              
      IF (NDIM0.EQ.2) THEN
        IF (AXE.EQ.'XX') THEN      
         DO JJ = 1,INELS
             idl=i*(i+1)/2 + 1
             DO II=1,INEL
               VKE(idl)=VKE(idl) - VNI(IN1+II)*VNIS(IN2+JJ)*C*(INORM)
               idl = idl + 4
             ENDDO
             DO II=1,JJ
               VKE(idl)= VKE(idl) + (VNIS(IN2+II)*VNIS(IN2+JJ))*C
               idl = idl + 2
             ENDDO
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 1
               VKE(idl)=VKE(idl) - VNI(IN1+II)*VNIS(IN2+JJ)*C*(H1/2)
               idl = idl + 3
             ENDDO             
             DO II=1,JJ
               idl = idl + 1
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 1
             ENDDO
             i=i+1
         ENDDO
        ELSEIF (AXE.EQ.'YY') THEN         
         DO JJ = 1,INELS
             idl=i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 1
               VKE(idl)=VKE(idl) - VNI(IN1+II)*VNIS(IN2+JJ)*C*(INORM)
               idl = idl + 3
             ENDDO
             DO II=1,JJ
               VKE(idl)= VKE(idl) + (VNIS(IN2+II)*VNIS(IN2+JJ))*C
               idl = idl + 2
             ENDDO
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               VKE(idl)=VKE(idl) - VNI(IN1+II)*VNIS(IN2+JJ)*C*(H1/2)
               idl = idl + 4
             ENDDO             
             DO II=1,JJ
               idl = idl + 1
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 1
             ENDDO
             i=i+1
         ENDDO
        ENDIF 
c.....................
c         NDIM = 3    
c.....................
      ELSEIF (NDIM0.EQ.3) THEN
C      IF     ((AXE.EQ.'XX').AND.(ABS(IT1).EQ.2)) THEN 
       IF     ((AXE.EQ.'XX').AND.(ABS(IT1).EQ.2)) THEN
        DO JJ=1,INELS                
             idl=i*(i+1)/2 + 1
             DO II=1,INEL
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(INORM)
               idl = idl + 6
             ENDDO
             DO II=1,JJ
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C
               idl = idl + 3
             ENDDO   
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 1
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(H1/2)
               idl = idl + 5
             ENDDO
             DO II=1,JJ
               idl = idl + 1
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 2
             ENDDO  
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 2
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(H2/2)
               idl = idl + 4
             ENDDO
             DO II=1,JJ
               idl = idl + 2
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 1
             ENDDO  
             i=i+1
         ENDDO                       
C      ELSEIF ((AXE.EQ.'XX').AND.(ABS(IT1).EQ.3)) THEN
       ELSEIF ((AXE.EQ.'XX').AND.(ABS(IT1).EQ.3)) THEN
        DO JJ=1,INELS                
             idl=i*(i+1)/2 + 1
             DO II=1,INEL
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(INORM)
               idl = idl + 6
             ENDDO
             DO II=1,JJ
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C
               idl = idl + 3
             ENDDO   
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 2
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(H1/2)
               idl = idl + 4
             ENDDO
             DO II=1,JJ
               idl = idl + 1
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 2
             ENDDO  
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 1
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(H2/2)
               idl = idl + 5
             ENDDO
             DO II=1,JJ
               idl = idl + 2
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 1
             ENDDO  
             i=i+1
         ENDDO                       
C      ELSEIF ((AXE.EQ.'YY').AND.(ABS(IT1).EQ.1)) THEN
       ELSEIF ((AXE.EQ.'YY').AND.(ABS(IT1).EQ.1)) THEN
        DO JJ=1,INELS                
             idl=i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 1
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(INORM)
               idl = idl + 5
             ENDDO
             DO II=1,JJ
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C
               idl = idl + 3
             ENDDO   
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(H1/2)
               idl = idl + 6
             ENDDO
             DO II=1,JJ
               idl = idl + 1
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 2
             ENDDO  
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 2
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(H2/2)
               idl = idl + 4
             ENDDO
             DO II=1,JJ
               idl = idl + 2
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 1
             ENDDO  
             i=i+1
         ENDDO                
C      ELSEIF ((AXE.EQ.'YY').AND.(ABS(IT1).EQ.3)) THEN
       ELSEIF ((AXE.EQ.'YY').AND.(ABS(IT1).EQ.3)) THEN
        DO JJ=1,INELS                
             idl=i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 1
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(INORM)
               idl = idl + 5
             ENDDO
             DO II=1,JJ
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C
               idl = idl + 3
             ENDDO   
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 2
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(H1/2)
               idl = idl + 4
             ENDDO
             DO II=1,JJ
               idl = idl + 1
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 2
             ENDDO  
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(H2/2)
               idl = idl + 6
             ENDDO
             DO II=1,JJ
               idl = idl + 2
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 1
             ENDDO  
             i=i+1
         ENDDO               
C      ELSEIF ((AXE.EQ.'ZZ').AND.(ABS(IT1).EQ.1)) THEN
       ELSEIF ((AXE.EQ.'ZZ').AND.(ABS(IT1).EQ.1)) THEN
        DO JJ=1,INELS                
             idl=i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 2
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(INORM)
               idl = idl + 4
             ENDDO
             DO II=1,JJ
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C
               idl = idl + 3
             ENDDO   
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(H1/2)
               idl = idl + 6
             ENDDO
             DO II=1,JJ
               idl = idl + 1
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 2
             ENDDO  
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 1
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(H2/2)
               idl = idl + 5
             ENDDO
             DO II=1,JJ
               idl = idl + 2
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 1
             ENDDO  
             i=i+1
         ENDDO       
C      ELSEIF ((AXE.EQ.'ZZ').AND.(ABS(IT1).EQ.2)) THEN
       ELSEIF ((AXE.EQ.'ZZ').AND.(ABS(IT1).EQ.2)) THEN
        DO JJ=1,INELS                
             idl=i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 2
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(INORM)
               idl = idl + 4
             ENDDO
             DO II=1,JJ
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C
               idl = idl + 3
             ENDDO   
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               idl = idl + 1
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(H1/2)
               idl = idl + 5
             ENDDO
             DO II=1,JJ
               idl = idl + 1
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 2
             ENDDO  
             i=i+1
             idl = i*(i+1)/2 + 1
             DO II=1,INEL
               VKE(idl)=VKE(idl)- VNI(IN1+II)*VNIS(IN2+JJ)*C*(H2/2)
               idl = idl + 6
             ENDDO
             DO II=1,JJ
               idl = idl + 2
               VKE(idl)=VKE(idl)+ (VNIS(IN2+II)*VNIS(IN2+JJ))*C*(H/2)**2
               idl = idl + 1
             ENDDO  
             i=i+1
         ENDDO
       ENDIF        

      ENDIF            
      
      IN1=IN1+NDIM0*INEL
      IN2=IN2+NDIM0*INELS
      
  310 CONTINUE
 

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


      INEL=INEL0
      RETURN
      
C
C================================================================================
C----------MASS MATRIX
C================================================================================

500   DO I=1,IKE
        VKE(I)=ZERO
      ENDDO
      
      RETURN

C********************************************************************************
C********************************************************************************

C----------EVALUATE ELEMENT TANGENT MATRIX
400   RETURN

C    
C-------------EVALUATE THE ELEMENT RESIDUAL
600   RETURN

C
C------------EVALUATE FE
700   RETURN

C
C-------------EVALUATE AND PRINT GRADIENTS AT G.P.
800   RETURN
C
C ----------EVALUATE THE IMPEDANCE MATRIX
900   RETURN


C-----------EVALUATE THE SOURCE MATRIX
1000  RETURN


C----------EVALUTE THE LOSS MATRIX
1100  RETURN

      END

