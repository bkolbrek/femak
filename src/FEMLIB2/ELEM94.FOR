      SUBROUTINE ELEM94(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1      VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C********************************************************************C
C ELEM94.FOR                                                         C
C  16/3-1995     NTH                                                 C
C  Tonni F. Johansen                                                 C
C********************************************************************C
C  Coupling element: elastic solid - porous material using           C
C  Biot-theory                                                       C
C  Coupling always on one side of element.                           C
C  Speciality for FEMAK: separate nodes for elastic and porous domainC
C  Coupling using penalty function node by node                      C
C  Penalty number: VPREE(1)*1.D7                                     C
C                                                                    C
C  2, dimensions or axisymmetry                                      C
C     2 dimension: 6 nodes 3 -for solid, 3 for por.mat. (18 dofs)    C
C                  (3 noded isoparametric element)                   C
C     axisymmetry: 6 nodes 3 -for fluid, 3 for por.mat.  (18 dofs)   C
C                  (3 noded isoparametric element)                   C
C  Number of integration points: none                                C
C  Number of degrees of freedom per node:                            C
C                             elastic material : 2                   C
C                             porous material  : 4                   C
C  Element matrix or vector formed by this subbprogram               C
C  according to icode value :                                        C
C     ICODE.EQ.1  return of parameters                               C
C     ICODE.EQ.2  evaluate interpolation functions and               C
C                 numerical integration coeffisients                 C
C     ICODE.EQ.3  stiffness matrix (VKE)                             C
C     ICODE.EQ.5  not used                                           C
C     ICODE.EQ.6  not used                                           C
C     ICODE.EQ.8  not used                                           C
C  Element properties                                                C
C     VPREE(1) Young's modulus                                       C
C     VPREE(2) Poisson's coeffiesient                                C
C     VPREE(3)  specific mass                                        C
C********************************************************************C
C Theory and implementation of this element is described in FEMAK    C
C programmer's manual and user's manual                              C
C
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
      REAL*8 VNI(2160)
      INTEGER IPGKED(3)

      SAVE VCPG,VKPG,VNI


      INTEGER IDL,J0,LENR,IDLD,J1,I,I0,IG,I1,I2,INI,IKE
      INTEGER COL,ROW,PNTCO,PNTCO1,PNTRW1
      REAL*8 DETJ,EPS,ZERO,R,C
      REAL*8 NR(3),NRLEN,GKSI(3),GETA(3)
      
C---------NB!! Initiates VNDSRC somewhat tricky !!!      
C      EQUIVALENCE(VNDSRC(1),VA(LNDSRC))

C------------Number of Gauss points in ksi,eta,dzeta direction
      DATA IPGKED/3,3,3/
c      DATA IPGKED/2,2,2/
C---------------------------------------------------
      DATA EPS/1.D-6/,ZERO/0.0/
     
      
C
      IF(NSYM.EQ.0) THEN
        IKE=IDLE*(IDLE+1)/2
      ELSE
        IKE=IDLE*IDLE
      ENDIF

C-------------CHOOSE FUNCTION TO BE EXECUTED
      GO TO (100,200,300,400,500,600,700,800,900,1000,1100),ICODE

C
C------------RETURN ELEMENT PARAMETERS IN COMMON /RGDT/
100   IF (NDIM.EQ.2) THEN
        IDLE0=18
        INEL0=6
        IPG0=9
      ELSEIF (NDIM.EQ.3) THEN
        IDLE0=72
        INEL0=16
        IPG0=27
      ENDIF 


      CFLAG=0
      RETURN

C
C-----------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C-----------THEIR DERIVATIVES AT G.P.
c-----------NB! Do not really use G.P. but the nodes  NB!

200   CONTINUE
C..........Must introduce a "trick" to fool GAUSS to give a number of gauss-points according to the
C..........double set of nodes.
      IF (NDIM.EQ.2) THEN
         VKPG(1)=-1.0
         VKPG(2)=0.0
         VKPG(3)=1.0
         IPG=3                  
      ELSEIF (NDIM.EQ.3) THEN
         VKPG(1)=-1.0
         VKPG(2)=-1.0
         VKPG(3)=0.0
         VKPG(4)=-1.0
         VKPG(5)=1.0
         VKPG(6)=-1.0
         VKPG(7)=1.0
         VKPG(8)=0.0
         VKPG(9)=1.0
         VKPG(10)=1.0
         VKPG(11)=0.0
         VKPG(12)=1.0
         VKPG(13)=-1.0
         VKPG(14)=1.0
         VKPG(15)=-1.0
         VKPG(16)=0.0
         IPG=8
      ENDIF  

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


C..........Must introduce a "trick" to fool NI01 to give a number of shape functions according to the
C..........double set of nodes.
      IF (NDIM.EQ.2) THEN
        INEL=3
        NDIM=1
        CALL NI01(VKPG,VNI)
        INEL=6
        NDIM=2
      ELSEIF (NDIM.EQ.3) THEN
        INEL=8   
        NDIM=2
        CALL NI01(VKPG,VNI)
        INEL=16
        NDIM=3
      ENDIF  

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
        VKE(I)=ZERO
      ENDDO

C........The penalty number C
      C=VPREE(1)*1.D7

C..........Set to tricky number of nodes to get correct shape functions and Gauss points.
      IF (NDIM.EQ.2) THEN
        INEL=3
      ELSEIF (NDIM.EQ.3) THEN
        INEL=8
      ENDIF  

C--------Loop through all nodes.
      INI=1+INEL
      I2=0

C..........Note only loop through all nodes, therefore matrix pointer outside the 
C..........gauss point (node) loop.
      IDL=1

      DO  IG=1,IPG
C----------------EVALUATE THE JACOBIAN, ITS INVERSE AND ITS DETERMINANT
        CALL JACOB(VNI(INI),VCORE,NDIM,INEL,VJ,VJ1,DETJ)


c.........Find Jacobi determinant and normal to element
        IF (NDIM.EQ.2) THEN
C...........2D-case        

          DETJ=SQRT(VJ(3)*VJ(3)+VJ(1)*VJ(1))
          NR(1)=VJ(3)/DETJ    
          IF (ABS(NR(1)).LT.EPS) THEN
            NR(1)=ZERO
          ENDIF
          NR(2)=-VJ(1)/DETJ
          IF (ABS(NR(2)).LT.EPS) THEN
            NR(2)=ZERO
          ENDIF


        ELSEIF (NDIM.EQ.3) THEN

C.........3D-case

C.........Find dR/dksi and dR/deta, tangents to element in node
C.........GKSI(1) = dx/dksi, GKSI(2) = dy/dksi, GKSI(3) = dz/dksi
C.........GETA(1) = dx/deta, GETA(2) = dy/deta, GETA(3) = dz/deta
          GKSI(1)=0.0
          GKSI(2)=0.0
          GKSI(3)=0.0
          GETA(1)=0.0
          GETA(2)=0.0
          GETA(3)=0.0
          DO I=1,8
            GKSI(1)=GKSI(1)+VNI(INI+I-1)*VCORE(3*I-2)
            GKSI(2)=GKSI(2)+VNI(INI+I-1)*VCORE(3*I-1)
            GKSI(3)=GKSI(3)+VNI(INI+I-1)*VCORE(3*I)
            GETA(1)=GETA(1)+VNI(INI+I+7)*VCORE(3*I-2)
            GETA(2)=GETA(2)+VNI(INI+I+7)*VCORE(3*I-1)
            GETA(3)=GETA(3)+VNI(INI+I+7)*VCORE(3*I)
          ENDDO

C............Normal found from cross-product of tangents
          NR(1)=GKSI(2)*GETA(3)-GKSI(3)*GETA(2)
          NR(2)=GKSI(3)*GETA(1)-GKSI(1)*GETA(3)
          NR(3)=GKSI(1)*GETA(2)-GKSI(2)*GETA(1)
          NRLEN=SQRT(NR(1)*NR(1)+NR(2)*NR(2)+NR(3)*NR(3))
          NR(1)=NR(1)/NRLEN
          NR(2)=NR(2)/NRLEN
          NR(3)=NR(3)/NRLEN
          

        ENDIF

C........Write results for Jacobian and normal
        IF (M.GE.2) THEN
          IF(NDIM.EQ.2) THEN
            WRITE(MP,2050)VJ,VJ1,DETJ
            WRITE(MP,2055)NR
          ELSEIF (NDIM.EQ.3) THEN
C            WRITE(MP,2051)VJ,VJ1,DETJ
            WRITE(MP,2056)NR
          ENDIF
        ENDIF

2050  FORMAT(/' JACOBIAN=',4E12.5 / ' J INVERS=',
     1         4E12.5/' DETJ=',E12.5)
2051  FORMAT(/' JACOBIAN=',9E12.5 / ' J INVERS=',
     1         9E12.5/' DETJ=',E12.5)
2055  FORMAT(/' NORMAL=',E12.5,'i+',E12.5,'j')
2056  FORMAT(/' NORMAL=',E12.5,'i+',E12.5,'j+',E12.5,'k')



C----------------FORM MATRIX VKE, FOR ELASTIC FRAME
C.............................................
C.......J0 - points to ux - biot -material
C.......IDL - heigth of column ux -biot material
C.......IDLD - heigth of column ux -elastic material
C.......J1 - points to ux - elastic material
        J0=IDL*(IDL+1)/2

        IF (NDIM.EQ.2) THEN
C................Coupling elastic displacement and solid frame displacement
          VKE(J0) = C
          VKE(J0+IDL+1) = C
          LENR=2*NDIM*(INEL-IG+1)
          IDLD=NDIM*INEL*2+(IG-1)*NDIM+1
c          J1=LENR*(LENR+1)/2+IDL
          J1=IDLD*(IDLD-1)/2+IDL
          VKE(J1)=-C
          VKE(J1+IDLD+1)=-C
C................set relative displacement to zero
          J0=J0+NDIM*IDL+3              
          VKE(J0)=NR(1)*NR(1)*C
          J0=J0+IDL+2
          VKE(J0)=NR(1)*NR(2)*C
          VKE(J0+1)=NR(2)*NR(2)*C

C..............Elastic dofs' part
          J0=(IDLD)*(IDLD+1)/2
          VKE(J0)=C
          VKE(J0+IDLD+1)=C
        ELSEIF (NDIM.EQ.3) THEN
C................Coupling elastic displacement and solid frame displacement
          VKE(J0) = C
          VKE(J0+IDL+1) = C
          VKE(J0+2*IDL+3) = C
          LENR=2*NDIM*(INEL-IG+1)
          IDLD=NDIM*INEL*2+(IG-1)*NDIM+1
c          J1=LENR*(LENR+1)/2+IDL
          J1=IDLD*(IDLD-1)/2+IDL
          VKE(J1)=-C
          VKE(J1+IDLD+1)=-C
          VKE(J1+2*IDLD+3)=-C
C................set relative displacement to zero
          J0=J0+NDIM*IDL+6              
          VKE(J0)=NR(1)*NR(1)*C
          J0=J0+IDL+3
          VKE(J0)=NR(1)*NR(2)*C
          VKE(J0+1)=NR(2)*NR(2)*C
          J0=J0+IDL+4
          VKE(J0)=NR(1)*NR(3)*C
          VKE(J0+1)=NR(2)*NR(3)*C
          VKE(J0+2)=NR(3)*NR(3)*C


C..............Elastic dofs' part
          J0=(IDLD)*(IDLD+1)/2
          VKE(J0)=C
          VKE(J0+IDLD+1)=C
          VKE(J0+2*IDLD+3)=C
        ENDIF

C.............Point at next set of shapfunction and set height of next set of collomns in VKE.
        INI=INI+NDIM*INEL  
        IDL=IDL+2*NDIM
      ENDDO
      

C......Reset to correct number of nodes.        
      IF (NDIM.EQ.2) THEN
        INEL=6
      ELSEIF (NDIM.EQ.3) THEN
        INEL=16
      ENDIF  
      
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
400   RETURN

C
C----------MASS MATRIX


500   DO I=1,IKE
        VKE(I)=ZERO
      ENDDO
      
      RETURN


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

