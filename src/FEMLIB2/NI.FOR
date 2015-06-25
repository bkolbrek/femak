      SUBROUTINE NI(VKSI,KEXP,KDER,VP,VPN,VNI)
C**********************************************************************
C     NI.FOR
C**********************************************************************
C     28/2 1992
C     Kompaniet Numerisk Lyd
C     Ulf R. Kristiansen
C     Programmert etter Dhat og Touzot s.55
C**********************************************************************
C     TO EVALUATE FUNCTIONS N AND THEIR DERIVATIVES 
C     AT POINT VKSI ON THE REFERENCE ELEMENT
C       INPUT    VKSI,KEXP,KDER,VP,VPN,IDLE,M,MP
C       OUTPUT   VNI
C======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
c      COMMON/COOR/NDIM
c      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,IPG
c      COMMON/ES/M,MR,MP

      include 'femak.cmn'    

      DIMENSION VKSI(1),KEXP(1),KDER(1),VP(1),VPN(1),VNI(1)
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C------  COMPUTE THE POLYNOMIAL BASIS AT OINT VKSI
      CALL BASEP(VKSI,KEXP,KDER,VP)
C------  P*(PN-INVERSE)  PRODUCT
      I0=1
      DO 20 IJ=1,INEL
      I1=I0
      C=ZERO
      DO 10 II=1,INEL
      C=C+VP(II)*VPN(I1)
10    I1=I1+1
      VNI(IJ)=C
20    I0=I0+INEL
C------  PRINT FUNCTIONS N
      IF(M.LT.3) GO TO 30
      WRITE(MP,2000) (KDER(I),I=1,NDIM)
2000  FORMAT(/' DERIVATIVE OF N WITH ORDER ',3I2)
      WRITE(MP,2010) (VKSI(I),I=1,NDIM)
2010  FORMAT(14X,'AT POINT ',3E13.5)
      WRITE(MP,2020) (VNI(I),I=1,INEL)
2020  FORMAT(/(1X,10E13.5))
30    RETURN
      END
