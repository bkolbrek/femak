      SUBROUTINE EXCOOR(VCORG,KDLNC)
c****************************************************************************c 
C excoor.for                                                                 c
c   31/1-1992                                                                c
c    Kompaniet Numerisk lyd, Tonni F. Johansen                               c
c****************************************************************************c
c  Subroutine EXCOOR                                                         C
C   to execute block 'COOR'                                                  c
c   read nodal coordinates                                                   C
C     Programmed after Dhatt & Touzot p.385-386                              C
c****************************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

c      COMMON/COOR/NDIM,NNT,NDLN,NDLT,FAC(3)
c      COMMON/ES/M,MR,MP,M1
      COMMON/TRVL/X1(3),X2(3)

      include 'femak.cmn'
      
      DIMENSION VCORG(1),KDLNC(1)
      DATA SPECL/1.23456789D31/

C------------------------------------------------------------------------------
C------------Initialize coordinates
      I1=(NNT-1)*NDIM+1
      DO 10 I=1,I1,NDIM
10      VCORG(I)=SPECL

C----------READ NODAL DATA CARDS
      IF (M.GT.0) WRITE(MP,2000)
2000  FORMAT(//' NODAL DATA CARDS'/)

C20    READ(M1,1000) IN1,X1,IN2,X2,INCR,IDLN
20    READ(MR,1000) IN1,X1,IN2,X2,INCR,IDLN
1000  FORMAT(2(I5,3F10.0),2I5)
      IF (M.GT.0) WRITE(MP,2010) IN1,X1,IN2,X2,INCR,IDLN
2010  FORMAT('>>>>>',2(I5,3E12.5),2I5)
      IF(IN1.LE.0) GO TO 60

C-----------DECODE THE CARD
      IF (IN1.GT.NNT) CALL ERREUR(11,IN1,NNT,0)
      IF (IN2.GT.NNT) CALL ERREUR(12,IN2,NNT,0)
      IF (IN2.LE.0) IN2=IN1
      IF (IDLN.GT.NDLN) CALL ERREUR(13,IDLN,NDLN,0)
      IF (IDLN.LE.0) IDLN=NDLN
      IF (INCR.EQ.0) INCR=1

      I1=(IN2-IN1)/INCR
      I2=IN1+I1*INCR
      IF (I1.EQ.0) I1=1
      IF (IN2.NE.I2) CALL ERREUR(14,IN2,IN2,0)

C--------------GENERATE NODES BY INTERPOLATION
      DO 30 I=1,NDIM
        X1(I)=X1(I)*FAC(I)
        X2(I)=X2(I)*FAC(I)
30      X2(I)=(X2(I)-X1(I))/I1

      I1=0
      I2=(IN1-1)*NDIM+1
      I3=(INCR-1)*NDIM

      DO 50 IN=IN1,IN2,INCR
        KDLNC(IN+1)=IDLN
        IF (VCORG(I2).NE.SPECL) CALL ERREUR(15,IN,IN,0)
        DO 40 I=1,NDIM
          VCORG(I2)=X1(I)+X2(I)*I1
          I2=I2+1
40      CONTINUE
        I1=I1+1
        I2=I2+I3
50    CONTINUE

      GO TO  20

C--------------CHECK FOR MISSING NODES
60    I1=NNT*NDIM+1
      I2=0
      I3=NNT+1
      DO 90 I=1,NNT
        I1=I1-NDIM
        I3=I3-1
        IF (VCORG(I1)-SPECL) 70,80,70

70      IF (I2.EQ.0) I2=I3
        GO TO 90
80      IF (I2.EQ.0) CALL ERREUR(16,I3,I3,0)
        IF (I2.NE.0) CALL ERREUR(17,I3,I3,1)
90    CONTINUE
      IF (I2.NE.NNT) CALL ERREUR(18,NNT,I2,0)

C----------TOTAL NUMBERS OF D.O.F
      NDLT=0
      I1=NNT+1
      DO 100 I=2,I1
        NDLT=NDLT+KDLNC(I)
100   CONTINUE

C-----------OUTPUT
      IF (M.LT.2) GO TO 120
      WRITE(MP,2020)
2020  FORMAT(/10X,'NODE D.O.F.',5X,'X',11X,'Y',11X,'Z'/)
      I1=1
      I2=NDIM

      DO 110 IN=1,NNT
        WRITE(MP,2030)IN,KDLNC(IN+1),(VCORG(I),I=I1,I2)
2030    FORMAT(10X,2I5,3E12.5)
        I1=I1+NDIM
        I2=I2+NDIM
110   CONTINUE

120   RETURN
      END
