      SUBROUTINE PRSOL(KDLNC,VCORG,VDIMP,KNEQ,VFG)
C**************************************************************
C     PRSOL.FOR
C**************************************************************
C     25/2/1992
C     Ulf R. Kristiansen NTH
C     Programmert etter Dhatt og Touzot s. 416
C**************************************************************
C     TO PRINT THE SOLUTION
C==============================================================
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRVL/V(10),FX(10)

      include 'femak.cmn'    

      DIMENSION VDIMP(1),KDLNC(1),VCORG(1),KNEQ(1),VFG(1)
      DATA RF/4H *  /,RL/4H    /,ZERO/0.D0/
C--------------------------------------------------------------
      X2=ZERO
      X3=ZERO
      WRITE(MP,2000)
2000  FORMAT(/' NODES',4X,'X',11X,'Y',11X,'Z',10X,'DEGREES OF FREEDOM (*
     1 = PRESCRIBED)'/)
      I2=0
      DO 50 IN=1,NNT
      I1=I2+1
      I2=I2+NDIM
      ID1=KDLNC(IN)+1
      ID2=KDLNC(IN+1)
      ID=ID2-ID1+1
      IF(ID2.LT.ID1) GO TO 50
      X1=VCORG(I1)
      IF(NDIM.GE.2) X2=VCORG(I1+1)
      IF(NDIM.GE.3) X3=VCORG(I1+2)
      J=ID1
      DO 40 I=1,ID
      JJ=KNEQ(J)
      IF(JJ) 10,20,30
10    V(I)=VDIMP(-JJ)
      FX(I)=RF
      GO TO 40
20    V(I)=ZERO
      FX(I)=RF
      GO TO 40
30    V(I)=VFG(JJ)
      FX(I)=RL
40    J=J+1
      WRITE(MP,2010)IN,X1,X2,X3,(V(II),FX(II),II=1,ID)
2010  FORMAT(1X,I5,3E12.5,5X,5(E12.5,A4)/47X,5(E12.5,A4))
50    CONTINUE
      RETURN
      END
