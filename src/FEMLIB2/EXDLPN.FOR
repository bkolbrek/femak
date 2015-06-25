      SUBROUTINE EXDLPN(KDLNC)
C**********************************************************************C
C EXDLPN.FOR                                                           C
C  28/9-1994   Kompaniet Numerisk lyd                                  C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  To execute  Block 'DLPN'                                            C
C  Read number of DOF per node                                         C
C  Programed after Dhatt & Touzot, p.384                               C
C**********************************************************************C
      INTEGER KDLNC(1)


      INTEGER K1(15)
      COMMON/TRVL/K1
      
      include 'femak.cmn'
      
      INTEGER J,I
      
      IF (M.GT.0) WRITE(MP,2000)
2000  FORMAT(//' GROUP OF D.O.F. '/)


C---------READ A GROUP CARD
C10    READ(M1,1000)IDLN,K1
10    READ(MR,1000)IDLN,K1
1000  FORMAT(16I5)

      IF(M.GT.0) WRITE(MP,2010)IDLN,K1
2010  FORMAT('>>>>>',16I5)    

      IF (IDLN.GT.0) THEN
        IF (IDLN.GT.NDLN) CALL ERREUR(21,IDLN,NDLN,1)
C---------STORE D.O.F. NUMBERS
20      DO I=1,15
          J=K1(I)
          IF(J.LE.0) GOTO 10
          IF (J.GT.NNT) CALL ERREUR(22,J,NNT,1)
          KDLNC(J+1)=IDLN
        ENDDO      
              
C        READ(M1,1010)K1
        READ(MR,1010)K1
1010    FORMAT(5X,15I5)
        IF (M.GT.0) WRITE(MP,2020) K1
2020    FORMAT('>>>>>',5X,15I5)    
        GO TO 20
        
      ENDIF


C---------TOTAL NUMBER OF D.O.F.
      NDLT=0
      J=NNT+1
      DO I=2,J        
        NDLT=NDLT+KDLNC(I)  
      ENDDO      
      

      RETURN
      END
