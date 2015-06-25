      SUBROUTINE EXPRND(VPRNG)
C**********************************************************************C
C EXPRND.FOR							       C
C  4/2-1992	Kompaniet Numerisk lyd				       C
C  Tonni F. Johansen						       C
C**********************************************************************C
C  To execute	Block 'PRND'					       C
C  To read nodal properties					       C
C  Programmert etter Dhatt & Touzot, s.392			       C
C**********************************************************************C

      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
      
      INTEGER KV(16),ICOD(10)
      REAL*8 V(10),H(20)
      COMMON/TRVL/KV,V,H,ICOD

      REAL*8 VPRNG(1)
      
      REAL*8 VTOT
      INTEGER I,L7,L10,L16,NODGRP,IHOW
      
      DATA L10/10/,L7/7/,L16/16/

      IF (M.GE.0) WRITE(MP,2010)
2010  FORMAT(// 'NODAL PROPERTIES CARDS '/)

C----------READ A B.C. GROUP CARD : 10 CODES + PRESCRIBED VALUE
      DO WHILE (.TRUE.)
        READ(MR,1000)NODGRP,(V(I),I=2,L10)
        V(1)=NODGRP
1000  FORMAT(I5,9F10.0)
        IF (M.GE.0) WRITE(MP,2020)NODGRP,(V(I),I=1,L7)
2020  FORMAT('>>>>>>',I5,9E12.5)

C----------CHECK FOR A BLANK CARD
        VTOT=NODGRP
        DO  I=2,L10
         VTOT=VTOT+V(I)
        ENDDO

C.........Stop execution of rutine here !!!
        IF (VTOT.EQ.0.0) RETURN

C---------READ CODIFICATION OF THE NODE NUMBER
        IHOW=1
        DO WHILE (IHOW.NE.0)
          READ(MR,1040)IHOW
1040    FORMAT(I1)    
        
          IF (IHOW.NE.0) THEN
C---------READ NODE CARDS


C-----------DIFFERENT BEHAVIOUR AS FUNCTION OF IHOW
          IF (IHOW.EQ.1) THEN
50          READ(MR,1020)(KV(IN),IN=1,L16)
            IF (M.GE.0) WRITE(MP,2040)(KV(IN),IN=1,L16)
1020  FORMAT(16I6)
2040  FORMAT('>>>>>>',10X,16I6)
C-----------FORM NEQ

            IN=1
            DO WHILE(IN.LT.17.AND.KV(IN).GT.0)

C-----------END OF GROUP OF B.C. OR END OF NODES OR ANALYSIS OF A NODE
              IF (KV(IN).GT.0)THEN
                IF (I2.GT.NNT) THEN
                  CALL ERREUR(32,I2,NNT,1)
                ELSE                
                  DO I=1,NPRN
                    VPRNG((KV(IN)-1)*NPRN+I)=V(I)
                  ENDDO
                ENDIF
              ENDIF
                
              IN=IN+1                
            ENDDO
            IF (IN.EQ.17) GOTO 50
          ELSEIF (IHOW.EQ.2) THEN

            READ(MR,1020)(KV(IN),IN=1,L16)
            IF (M.GE.0) WRITE(MP,2040)(KV(IN),IN=1,L16)

            IF (KV(3).EQ.0)KV(3)=1

            DO I2=KV(1),KV(2),KV(3)
              IF (I2.GT.0) THEN
                IF (I2.GT.NNT) THEN
                  CALL ERREUR(32,I2,NNT,1)
                ELSE                
                  DO I=1,NPRN
                    VPRNG((I2-1)*NPRN+I)=V(I)
                  ENDDO
                ENDIF
              ENDIF
            ENDDO
          ENDIF
C---------ADDITIONNAL INITIAL CONDITIONS
        ENDIF
       
        ENDDO

        

      ENDDO

c   Dhatt & Touzot version:
C------------------------------------------------------------------------
C----------TO READ PROPERTIES NODEWISE
C      I1=NNT*NPRN
C      READ(MR,1000) (VPRNG(I),I=1,I1)
C1000  FORMAT(8F10.0)
C      IF (M.GE.0) WRITE(MP,2000) (VPRNG(I),I=1,I1)
C2000  FORMAT(//' CARDS OF NODAL PROPERTIES '/ (' >>>>>',8E12.5))
C
      RETURN
      END
