      SUBROUTINE PRPVTS(VKGD)
C************************************************************
C     PRPVTS.FOR
C************************************************************
C     21/2 1992
C     Ulf R. Kristiansen   NTH
C     Programmed after Dhatt and Touzot p.414
C************************************************************
C   TO EVALUATE AND TO PRINT THE PIVOTS AND THE DETERMINANT OF MATRIX KG
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)

      include 'femak.cmn'    

      DIMENSION VKGD(1)
      DATA UN/1.D0/,GROS/1.D38/
      ABS(X)=DABS(X)
C-----------------------------------------------------------------------
      X1=GROS
      X2=GROS
      DET=UN
      IDET=0
C------  PRINT PIVOTS OF MATRIX KG
      IF(M.GE.2) WRITE(MP,2000)(VKGD(I),I=1,NEQ)
2000  FORMAT(/' GLOBAL MATRIX PIVOTS'/(1X,10E12.5))
      DO 50 I=1,NEQ
C------  ABSOLUTE VALUE OF MINIMUM PIVOT
      X=ABS(VKGD(I))
      IF(X.GT.X1) GO TO 10
      X1=X
      I1=I
C------  ALGEBRAIC VALUE OF MINIMUM PIVOT
      X=VKGD(I)
10    IF(X.GT.X2) GO TO 20
      X2=X
      I2=I
C------  DETERMINANT (BOUNDS : EXPONENT + OR - 10)
20    DET=DET*VKGD(I)
30    DET1=ABS(DET)
      IF(DET1.LT.1.D10) GO TO 40
      DET=DET*1.D-10
      IDET=IDET+10
40    IF(DET1.GT.1.D-10) GO TO 50
      DET=DET*1.D10
      IDET=IDET-10
      GO TO 30
50    CONTINUE
C------  OUTPUT
      WRITE(MP,2010) X1,I1,X2,I2,DET,IDET
2010  FORMAT(/15X,'ABSOLUTE VALUE OF MINIMUM PIVOT   =',E12.5,' EQUATION
     1:',I5   /29X,              'ALGEBRAIC VALUE=',E12.5,' EQUATION:',
     2   I5   /29X,              'DETERMINANT      =',E12.5,' * 10 ** ',
     3   I5/)
      RETURN
      END
