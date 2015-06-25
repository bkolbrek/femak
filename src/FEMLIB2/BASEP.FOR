      SUBROUTINE BASEP(VKSI,KEXP,KDER,VP)
C******************************************************************
C     BASEP.FOR
C********************************************************************
C     28/2 1992
C     Ulf R. Kristiansen  NTH
C     Programmed after Dhatt and Touzot p.56
C*********************************************************************     
C     TO EVALUATE THE POLYNOMIAL BASIS AND ITS DERIVATIVES AT POINT 
C     VKSI
C       INPUT    VKSI,KEXP,KDER,IDLE,IDEG,NDIM,M,MP
C       OUTPUT   VP
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)

      include 'femak.cmn'

      DIMENSION VKSI(1),KEXP(1),KDER(1),VP(1)
      DIMENSION PUISS(3,10)
      DATA ZERO/0.D0/,UN/1.D0/
C-----------------------------------------------------------------------
C------  FORM SUCCESSIVE POWERS OF KSI,ETA,DZETA
      DO 10 I=1,NDIM
      PUISS(I,1)=UN
      DO 10 ID=1,IDEG
10    PUISS(I,ID+1)=PUISS(I,ID)*VKSI(I)
C------  DERIVATIVES OF ORDER KDER WITH RESPECT TO KSI,ETA,DZETA
      DO 50 IDL=1,INEL
      C1=UN
      I0=(IDL-1)*NDIM
      DO 30 I=1,NDIM
      IDR=KDER(I)
      I0=I0+1
      IXP=KEXP(I0)+1
      J=IXP-IDR
      IF(J.LE.0) GO TO 40
      IF(IDR.LE.0) GO TO 30
      DO 20 ID=1,IDR
20    C1=C1*(IXP-ID)
30    C1=C1*PUISS(I,J)
      GO TO 50
40    C1=ZERO
50    VP(IDL)=C1
C------  PRINT POLYNOMIAL BASIS
      IF(M.LT.4) GO TO 60
      WRITE(MP,2000) (KDER(I),I=1,NDIM)
2000  FORMAT(/' POLYNOMIAL BASIS, DERIVATIVE OF ORDER ',3I2)
      WRITE(MP,2010) (VKSI(I),I=1,NDIM)
2010  FORMAT(19X,'AT POINT ',3E13.5)
      WRITE(MP,2020) (VP(I),I=1,INEL)
2020  FORMAT(/(1X,10E12.5))
60    RETURN
      END

