      SUBROUTINE ELEM91(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1      VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
c      SUBROUTINE ELEM91(VCORE,VPRNE,VPREE,VDLE,VKE,VFE)
C**********************************************************************C
C ELEM91.FOR                                                           C
C  3/10-1994  Kompaniet Numerisk lyd                                   C
C  Marc Dhainaut                                                       C
C**********************************************************************C
C  NDIM = 2 : BEAM (3 nodes)/ACOUSTIC(3 nodes) COUPLING ELEMENT        C
C  NDIM = 3 : PLATE(9 nodes)/ACOUSTIC(8 nodes) COUPLING ELEMENT        C  
C----------------------------------------------------------------------C
C  Evaluate element information according to ICODE value               C
C     ICODE.EQ.1  return of parameters                                 C
C     ICODE.EQ.2  evaluate interpolation functions and                 C
C                 numerical integration coefficients                   C
C     ICODE.EQ.3  stiffnes matrix (VKE)                                C
C     ICODE.EQ.5  mass matrix (VKE)                                    C 
C----------------------------------------------------------------------C
C  Element properties                                                  C
C     VPREE(1) : fluid density                                         C
C     VPREE(2) : =0 the plate normal points towards the acoustic elt   C
C                =1 the plate normal points outside the acoustic elt   C
C----------------------------------------------------------------------C
C  Symmetry of the element                                             C
c     The matrcies (Ke, Me) of the element are NON symmetric           C
C**********************************************************************C

      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
c      COMMON/COOR/NDIM
c      COMMON/ASSE/NSYM,NKG,NKE,NDLE
c      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
c     &            IPG,ICODE,IDLE0,INEL0,IPG0
c      COMMON/ES/M,MR,MP

      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VKE1(1),VFE(1)
      COMPLEX*16 VIE(1),VIE1(1),VNDSRC(1),VPRSG(1),VSRE(1),VDLEC(1)
      INTEGER KTYPRSG(1)
c       DIMENSION VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VFE(1)
 
C-----------CHARACTERISTIC DIMENSIONS OF THE ELEMENT
C----------(VALID  UP TO 3 DIMENSIONS)
C     DIMENSION VCPG(IPG),VKPG(NDIM*IPG)
      DIMENSION VCPG(9),VKPG(18)
C     DIMENSION VJ(NDIM*NDIM),VJ1(NDIM*NDIM)
      DIMENSION VJ(4),VJ1(4)
C     DIMENSION VNI((1+NDIM)*INEL*IPG),IPGKED(NDIM)
      DIMENSION VNI(216),VNIS(243),IPGKED(2)
C     DIMENSION VCORE1(3*INEL)
      DIMENSION VCORE1(24)


      SAVE VCPG,VKPG,VNI,VNIS
C
      REAL*8 RHO,ISENS,C      
C------------ Number of Gauss points in ksi,eta direct
      DATA IPGKED/3,3/
C---------------------------------------------------
      DATA ZERO/0.0D0/,DEUX/2.D0/,X05/0.5D0/,RADN/.572957795130823D2/
      DATA EPS/1.D-6/

      SQRT(X)=DSQRT(X)
      ATAN2(X,Y)=DATAN2(X,Y)

      NDIM0=NDIM
      INEL0=INEL
            
C----------ELASTIC PARAMETERS         
      RHO  = VPREE(1)
      ISENS= VPREE(2)
      
C----------TEST
C      WRITE(*,*) 'ELEM91 : ICODE =',ICODE
      
C     
C-------------CHOOSE FUNCTION TO BE EXECUTED 

      GO TO (100,200,300,400,500),ICODE

C 
C================================================================================
C------------RETURN ELEMENT PARAMETERS IN COMMON /RGDT/ 
C================================================================================
100   IF (NDIM.EQ.2) THEN
         IDLE0=9
         INEL0=6
         IPG0=3
      ENDIF
      IF (NDIM.EQ.3) THEN
         IDLE0=35
         INEL0=17
         IPG0=9
      ENDIF

      CFLAG=0

      NDIM=NDIM0
      INEL=INEL0      
      RETURN

C
C================================================================================
C-----------EVALUATE COORDINATES, WEIGHTS, FUNCTIONS N AND
C-----------THEIR DERIVATIVES AT G.P.
C================================================================================
200   IF (NDIM0.EQ.2) THEN
          NDIM=1
          INEL=3
          INELS=3
      ENDIF    
      IF (NDIM0.EQ.3) THEN
          NDIM=2
          INEL=8
          INELS=9
      ENDIF    
        
      CALL GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
      IF (M.GE.2) THEN
         WRITE(MP,2000)IPG
         I0=1
         DO IG=1,IPG
            I1=I0+NDIM-1
            WRITE(MP,2010)VCPG(IG),(VKPG(I),I=I0,I1)
            I0=I0+NDIM
         ENDDO
      ENDIF
      
      CALL NI01(VKPG,VNI)
      IF (M.GE.2) THEN
         I1=NDIM0*INEL*IPG
         WRITE(MP,2020) (VNI(I),I=1,I1)
      ENDIF
      
      INEL=INELS
      CALL NI08(VKPG,VNIS)
      IF (M.GE.2) THEN
          I1=NDIM0*INEL*IPG
          WRITE(MP,2020) (VNIS(I),I=1,I1)
      ENDIF
      
2000  FORMAT(/I5,'  GAUSS POINTS'/10X,'VCPG',25X,'VKPG')
2010  FORMAT(1X,F20.15,5X,3F20.15)
2020  FORMAT(/'FUNCTIONS N AND DERIVATIVES'/ (1X,8E12.5))
      
      INEL=INEL0
      NDIM=NDIM0      
      RETURN
      
C      
C================================================================================
C----------COMPUTE ELEMENT "STIFFNESS" MATRIX
C================================================================================
C------ INITIALIZE VKE
300   DO I=1,IDLE*IDLE
        VKE(I)=ZERO
      ENDDO
      
C------ DIMENSION NDIM-1      
      IF (NDIM0.EQ.2) THEN
          NDIM=1
          INEL=3
          INELS=3
      ENDIF    
      IF (NDIM0.EQ.3) THEN
          NDIM=2
          INEL=8
          INELS=9
      ENDIF
          
C------ CHANGE NDIM COORDINATES TO NDIM-1 COORDINATES
      CALL COOR23D(VCORE,VCORE1)
              
C------ LOOP OVER G.P
      I1=1+INEL
      IA=1
      IS=1
      
      DO 330 IG=1,IPG 
      
C----------------EVALUATE THE JACOBIAN, ITS INVERSE AND ITS DETERMINANT
      CALL JACOB(VNI(I1),VCORE1,NDIM,INEL,VJ,VJ1,DETJ)
      IF (DETJ.LT.EPS) WRITE(MP,2040)IEL,IG,DETJ
2040  FORMAT(' *** ELEM ',I5,' G.P. ',I3,' DET(J)=',E12.5)
      IF (M.GE.2) WRITE(MP,2050)VJ,VJ1,DETJ
2050  FORMAT(/' JACOBIAN=',4E12.5 / ' J INVERS=',
     1        4E12.5/' DETJ=',E12.5)
     
C----------------PERFORM COEF
      C=-VCPG(IG)*DETJ
      IF (ISENS.EQ.0) C=-C
      
C----------------FORM MATRIX VKE
      I0=IDLE*(IDLE-INEL)+1
c or  io=idle*(  inels  )+1      
      DO J=1,INEL
         DO I=1,INELS
            VKE(I0) = VKE(I0) + VNIS(IS+I-1)*VNI(IA+J-1)*C
            I0=I0+NDIM0
         ENDDO
         I0=I0+INEL
      ENDDO
            
      I1=I1+NDIM0*INEL
      IA=IA+NDIM0*INEL
330   IS=IS+NDIM0*INELS

      INEL=INEL0
      NDIM=NDIM0
      RETURN 
       
C             
C================================================================================
C----------COMPUTE ELEMENT "MASS" MATRIX
C================================================================================
C------ INITIALIZE VKE
500   DO I=1,IDLE*IDLE
        VKE(I)=ZERO
      ENDDO

C------ DIMENSION NDIM-1      
      IF (NDIM0.EQ.2) THEN
          NDIM=1
          INEL=3
          INELS=3
      ENDIF    
      IF (NDIM0.EQ.3) THEN
          NDIM=2
          INEL=8
          INELS=9
      ENDIF
          
C------ CHANGE NDIM COORDINATES TO NDIM-1 COORDINATES
      CALL COOR23D(VCORE,VCORE1)
      
C------ LOOP OVER G.P
      I1=1+INEL
      IA=1
      IS=1
      
      DO 530 IG=1,IPG
      
C----------------EVALUATE THE JACOBIAN, ITS INVERSE AND ITS DETERMINANT
      CALL JACOB(VNI(I1),VCORE1,NDIM,INEL,VJ,VJ1,DETJ)
      IF (DETJ.LT.EPS) WRITE(MP,2060)IEL,IG,DETJ
2060  FORMAT(' *** ELEM ',I5,' G.P. ',I3,' DET(J)=',E12.5)
      IF (M.GE.2) WRITE(MP,2070)VJ,VJ1,DETJ
2070  FORMAT(/' JACOBIAN=',4E12.5 / ' J INVERS=',
     1          4E12.5/' DETJ=',E12.5)
     
C----------------PERFORM RHO*COEF
      C=RHO*VCPG(IG)*DETJ
      IF (ISENS.EQ.0) C=-C

C----------------FORM MATRIX VKE
      I0=IDLE-INEL+1
c or  i0=  inels  +1    
      DO J=1,INELS
         DO I=1,INEL
            VKE(I0) = VKE(I0) + VNI(IA+I-1)*VNIS(IS+J-1)*C
            I0=I0+1
         ENDDO
         I0=I0+NDIM0*IDLE-INEL
      ENDDO      

      I1=I1+NDIM0*INEL
      IA=IA+NDIM0*INEL
530   IS=IS+NDIM0*INELS

      INEL=INEL0
      NDIM=NDIM0
      RETURN      

C********************************************************************************************
C********************************************************************************************      
                         
C----------EVALUATE ELEMENT TANGENT MATRIX
400   NDIM=NDIM0
      INEL=INEL0
      RETURN
            
      END
