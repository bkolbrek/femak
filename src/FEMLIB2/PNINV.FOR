      SUBROUTINE PNINV(VKSI,KEXP,VP,K1,VPN)
C********************************************************************
C     PNINV.FOR
C********************************************************************
C     28/2/1992
C     Kompaniet Numerisk Lyd 
C     Ulf R. Kristiansen
C     Programmert etter Dhatt og Touzot s.54
C********************************************************************
C     EVALUATE THE PN-INVERSE MATRIX WHICH
C     CONTAINS THE COEFFICIENTS OF FUNCTIONS N
C       INPUT        VKSI,KEXP,INEL,IDLE,ITPE,M,MP
C       WORKSPACE    VP,K1
C       OUTPUT       VPN
C====================================================================
      IMPLICIT REAL*8(A-H,O-Z)
c      COMMON/COOR/NDIM
c      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,IPG
c      COMMON/ES/M,MR,MP

      include 'femak.cmn'    

      DIMENSION VKSI(1),KEXP(1),VP(1),K1(1),VPN(1),KDER(3)
      DATA ZERO/0.D0/,KDER/3*0/
C-----------------------------------------------------------------------
C
C......  FORM PN MATRIX (FOR ANY LAGRANGE TYPE ELEMENT)
C
C
      I0=1
      I1=1
      DO 20 IN=1,INEL
      CALL BASEP(VKSI(I1),KEXP,KDER,VP)
      I2=I0
      DO 10 IJ=1,INEL
      VPN(I2)=VP(IJ)
10    I2=I2+INEL
      I0=I0+1
20    I1=I1+NDIM
C
C......  END OF PN INFORMATION
C
C------ PRINT THE PN MATRIX
      IF(M.LT.4) GO TO 40
      WRITE(MP,2000)
2000  FORMAT(/' PN MATRIX'/)
      ID=(INEL-1)*INEL
      DO 30 I0=1,INEL
      I1=I0+ID
30    WRITE(MP,2010) (VPN(IJ),IJ=I0,I1,INEL)
2010  FORMAT(1X,10E13.5/(14X,9E13.5))
C------  INVERSE OF THE PN MATRIX
40    CALL INVERS(VPN,INEL,INEL,K1,DET)
      IF(DET.NE.ZERO) GO TO 50
      WRITE(MP,2020) ITPE
2020  FORMAT(' *** ERROR, PN SINGULAR, ELEMENT TYPE:',I3)
      STOP
C------  PRINT THE PN-INVERSE MATRIX
50    IF(M.LT.4) GO TO 70
      WRITE(MP,2030)
2030  FORMAT(/' PN INVERSE MATRIX'/)
      DO 60 I0=1,INEL
      I1=I0+ID
60    WRITE(MP,2010) (VPN(IJ),IJ=I0,I1,INEL)
70    RETURN
      END


