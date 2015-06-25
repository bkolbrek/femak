      SUBROUTINE ELEM93(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1      VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C********************************************************************C
C ELEM93.FOR                                                         C
C  15/6-1996     NTH                                                  C
C  Tonni F. Johansen  & Marc Dhainaut                                C
C********************************************************************C
C  Coupling element: acoustic fluid - elastic solid                  C
C  Coupling always on one side of element.                           C
C  Speciality for FEMAK: separate nodes for fluid and porous domain. C
C                                                                    C
C  2, 3 dimensions or axisymmetry                                    C
C     2 dimension: 6 nodes 3 -for fluid, 3 for elastic solid (9 dofs)C
C                  (8 noded isoparametric element)                   C
C     axisymmetry: 6 nodes 3 -for fluid, 3 for elastic solid (9 dofs)C
C                  (8 noded isoparametric element)                   C
C     3 dimension: 16 nodes 8 -for fluid, 8 for elastic solid(32 dof)C
C                  (20 noded isoparametric element)                  C
C  Number of integration points: 3 in each direction                 C
C  Number of degrees of freedom per node:                            C
C                             fluid           : 1                    C
C                             elastic solid   : 2 or 4               C
C  Element matrix or vector formed by this subbprogram               C
C  according to icode value :                                        C
C     ICODE.EQ.1  return of parameters                               C
C     ICODE.EQ.2  evaluate interpolation functions and               C
C                 numerical integration coeffisients                 C
C     ICODE.EQ.3  stiffness matrix (VKE)                             C
C     ICODE.EQ.5  mass matrix (VKE)                                  C
C     ICODE.EQ.6  K.U. product (VFE)  .....not written.....NB!!      C
C     ICODE.EQ.8  print gradients     .....not written.....NB!!      C
C  Element properties                                                C
C     VPREE(1) Density of fluid                                      C
C********************************************************************C
C Theory and implementation of this element is described in FEMAK    C
C programmer's manual and user's manual                              C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
C   $DECLARE      
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


      INTEGER I,J,I0,IG,I1,I2,INI
      INTEGER IKE
      REAL*8 PI,DETJ,EPS,R,ZERO
      COMPLEX*16 C
      REAL*8 NR(3),NRLEN,GKSI(3),GETA(3)

C---------NB!! Initiates VNDSRC somewhat tricky !!!      
C      EQUIVALENCE(VNDSRC(1),VA(LNDSRC))

C------------Number of Gauss points in ksi,eta,dzeta direction
      DATA IPGKED/3,3,3/
c      DATA IPGKED/2,2,2/
C---------------------------------------------------
      DATA EPS/1.D-6/,ZERO/0.0/
      DATA PI/3.141592654/
     
      
C
      IKE=IDLE*IDLE

C-------------CHOOSE FUNCTION TO BE EXECUTED
      GO TO (100,200,300,400,500,600,700,800,900,1000,1100),ICODE

C
C------------RETURN ELEMENT PARAMETERS IN COMMON /RGDT/
100   IF (NDIM.EQ.2) THEN
        IDLE0=9
        INEL0=6
        IPG0=9
      ELSEIF (NDIM.EQ.3) THEN
        IDLE0=32
        INEL0=16
        IPG0=27
      ENDIF 


      CFLAG=0
      RETURN

C
C-----------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C-----------THEIR DERIVATIVES AT G.P.

200   CONTINUE
C..........Must introduce a "trick" to fool GAUSS to give a number of gauss-points according to the
C..........double set of nodes.
      IF (NDIM.EQ.2) THEN
        NDIM=1
        CALL GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
        NDIM=2
      ELSEIF (NDIM.EQ.3) THEN
        NDIM=2
        CALL GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
        NDIM=3
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

C..........Set to tricky number of nodes to get correct shape functions and Gauss points.
      IF (NDIM.EQ.2) THEN
        INEL=3
      ELSEIF (NDIM.EQ.3) THEN
        INEL=8
      ENDIF  

C--------Loop through all gauss-points.
      INI=1+INEL
      I2=0
      DO  IG=1,IPG
C----------------EVALUATE THE JACOBIAN, ITS INVERSE AND ITS DETERMINANT
        CALL JACOB(VNI(INI),VCORE,NDIM,INEL,VJ,VJ1,DETJ)

        IF (NDIM.EQ.2) THEN

C-----------If axisymmetry; find radius
          IF (IAX.NE.0) THEN
            R=0.0
            DO I=1,INEL
              R=R+VNI(6*(IG-1)+I)*VCORE((I-1)*NDIM+2)
            ENDDO
          ENDIF

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

C.........Find dR/dksi and dR/deta
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
          NR(1)=GKSI(2)*GETA(3)-GKSI(3)*GETA(2)
          NR(2)=GKSI(3)*GETA(1)-GKSI(1)*GETA(3)
          NR(3)=GKSI(1)*GETA(2)-GKSI(2)*GETA(1)
          NRLEN=SQRT(NR(1)*NR(1)+NR(2)*NR(2)+NR(3)*NR(3))
          NR(1)=NR(1)/NRLEN
          NR(2)=NR(2)/NRLEN
          NR(3)=NR(3)/NRLEN
          
          DETJ=(VJ(4)*VJ(8)-VJ(7)*VJ(5))**2
          DETJ=DETJ+(VJ(7)*VJ(2)-VJ(1)*VJ(8))**2
          DETJ=DETJ+(VJ(1)*VJ(5)-VJ(4)*VJ(2))**2
          DETJ=SQRT(DETJ)          

        ENDIF

        IF (DETJ.LT.EPS) WRITE(MP,2040)IEL,IG,DETJ
2040  FORMAT(' *** ELEM ',I5,' G.P. ',I3,' DET(J)=',E12.5)

C........Write results for Jacobian and normal
        IF (M.GE.2) THEN
          IF(NDIM.EQ.2) THEN
            WRITE(MP,2050)(VJ(I),I=1,4),(VJ1(I),I=1,4),DETJ
            WRITE(MP,2055)(NR(I),I=1,2)
          ELSEIF (NDIM.EQ.3) THEN
            WRITE(MP,2051)(VJ(I),I=1,9),(VJ1(I),I=1,9),DETJ
            WRITE(MP,2056)(NR(I),I=1,3)
          ENDIF
        ENDIF

2050  FORMAT(/' JACOBIAN=',4E12.5 / ' J INVERS=',
     1         4E12.5/' DETJ=',E12.5)
2051  FORMAT(/' JACOBIAN=',9E12.5 / ' J INVERS=',
     1         9E12.5/' DETJ=',E12.5)
2055  FORMAT(/' NORMAL=',E12.5,'i+',E12.5,'j')
2056  FORMAT(/' NORMAL=',E12.5,'i+',E12.5,'j+',E12.5,'k')


C----------------PERFORM COEF
        C=-VCPG(IG)*DETJ
        IF (IAX.NE.0) THEN
          C=C*R*2*PI
        ENDIF


C----------------FORM MATRIX VKE, FOR ELASTIC FRAME
        I0=IDLE*(IDLE-INEL)+1
        

C.........LOOP THROUGH ALL COLOUMNS (J)
        DO J=1,INEL

C.........LOOP THROUGH ALL ROWS (I)
          DO I=1,INEL
            IF (NDIM.EQ.2) THEN
              VKE(I0) = VKE(I0) + VNI(I2+I)*VNI(I2+J)*C*NR(1)
              VKE(I0+1) = VKE(I0+1)+VNI(I2+I)*VNI(I2+J)*C*NR(2)
            ELSEIF (NDIM.EQ.3) THEN
              VKE(I0) = VKE(I0) + VNI(I2+I)*VNI(I2+J)*C*NR(1)
              VKE(I0+1) = VKE(I0+1)+VNI(I2+I)*VNI(I2+J)*C*NR(2)
              VKE(I0+2) = VKE(I0+2)+VNI(I2+I)*VNI(I2+J)*C*NR(3)
            ENDIF
            I0=I0+NDIM
          ENDDO
          I0=I0+INEL
        ENDDO      
        I2=I2+NDIM*INEL
        INI=INI+NDIM*INEL
      ENDDO
      

C......Reset to correct number of nodes.        
      IF (NDIM.EQ.2) THEN
        INEL=6
      ELSEIF (NDIM.EQ.3) THEN
        INEL=16
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

C..........Set to tricky number of nodes to get correct shape functions and Gauss points.
      IF (NDIM.EQ.2) THEN
        INEL=3
      ELSEIF (NDIM.EQ.3) THEN
        INEL=8
      ENDIF  

C--------Loop through all gauss-points.
      INI=1+INEL
      I2=0
      DO  IG=1,IPG
C----------------EVALUATE THE JACOBIAN, ITS INVERSE AND ITS DETERMINANT
        CALL JACOB(VNI(INI),VCORE,NDIM,INEL,VJ,VJ1,DETJ)

        IF (NDIM.EQ.2) THEN
C-----------If axisymmetry; find radius
          IF (IAX.NE.0) THEN
            R=0.0
            DO I=1,INEL
              R=R+VNI(6*(IG-1)+I)*VCORE((I-1)*NDIM+2)
            ENDDO
          ENDIF
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

C.........Find dR/dksi and dR/deta
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
          NR(1)=GKSI(2)*GETA(3)-GKSI(3)*GETA(2)
          NR(2)=GKSI(3)*GETA(1)-GKSI(1)*GETA(3)
          NR(3)=GKSI(1)*GETA(2)-GKSI(2)*GETA(1)
          NRLEN=SQRT(NR(1)*NR(1)+NR(2)*NR(2)+NR(3)*NR(3))
          NR(1)=NR(1)/NRLEN
          NR(2)=NR(2)/NRLEN
          NR(3)=NR(3)/NRLEN
          
          DETJ=(VJ(4)*VJ(8)-VJ(7)*VJ(5))**2
          DETJ=DETJ+(VJ(7)*VJ(2)-VJ(1)*VJ(8))**2
          DETJ=DETJ+(VJ(1)*VJ(5)-VJ(4)*VJ(2))**2
          DETJ=SQRT(DETJ)          

        ENDIF

        IF (DETJ.LT.EPS) WRITE(MP,2040)IEL,IG,DETJ
C........Write results for Jacobian and normal
        IF (M.GE.2) THEN
          IF(NDIM.EQ.2) THEN
            WRITE(MP,2050)(VJ(I),I=1,4),(VJ1(I),I=1,4),DETJ
            WRITE(MP,2055)(NR(I),I=1,2)
          ELSEIF (NDIM.EQ.3) THEN
            WRITE(MP,2051)(VJ(I),I=1,9),(VJ1(I),I=1,9),DETJ
            WRITE(MP,2056)(NR(I),I=1,3)
          ENDIF
        ENDIF

C----------------PERFORM COEF
C        C=-W2*VCPG(IG)*DETJ*VPREE(1)
        C=VCPG(IG)*DETJ*VPREE(1)
        IF (IAX.NE.0) THEN
          C=C*R*2*PI
        ENDIF

C----------------FORM MATRIX VKE, FOR ELASTIC FRAME
        I0=IDLE-INEL+1
        

C.........LOOP THROUGH ALL COLOUMNS (J)
        DO J=1,INEL

C.........LOOP THROUGH ALL ROWS (I)
          DO I=1,INEL
            IF (NDIM.EQ.2) THEN
              VKE(I0) = VKE(I0) + VNI(I2+I)*VNI(I2+J)*C*NR(1)
              VKE(I0+IDLE) = VKE(I0+IDLE)+VNI(I2+I)*VNI(I2+J)*C*NR(2)
            ELSEIF (NDIM.EQ.3) THEN
              VKE(I0) = VKE(I0) + VNI(I2+I)*VNI(I2+J)*C*NR(1)
              VKE(I0+IDLE) = VKE(I0+IDLE)+VNI(I2+I)*VNI(I2+J)*C*NR(2)
              VKE(I0+2*IDLE) = VKE(I0+2*IDLE)+
     1                              VNI(I2+I)*VNI(I2+J)*C*NR(3)
            ENDIF
            I0=I0+1
          ENDDO
          I0=I0+NDIM*IDLE-INEL
        ENDDO      
        I2=I2+NDIM*INEL
        INI=INI+NDIM*INEL
      ENDDO
      

C......Reset to correct number of nodes.        
      IF (NDIM.EQ.2) THEN
        INEL=6
      ELSEIF (NDIM.EQ.3) THEN
        INEL=16
      ENDIF  


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

