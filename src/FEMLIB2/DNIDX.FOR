      SUBROUTINE DNIDX(VNI,VJ1,NDIM,INEL,VNIX)
C********************************************************************C
C    DNIDX.FOR                                                       C
C     28/2-1992                                                      C
C   Tonni F. Johansen NTH                                            C
C********************************************************************C
C   Compute the derivatives of interpolation functions with respect  C
C   to x,y,z  (1, 2 and 3 dimensions)                                C                   C
C    input                                                           C 
C       VNI    derivatives of interpolation function                 C
c              W.R.T. ksi,eta,  dzeta                                C 
C       VJ1     invers of jacobian                                   C
C       NDIM    number of dimensions (1, 2 or 3)                     C 
C       INEL    number of nodes per element                          C
C     output                                                         C 
C        VNIX   x,y,z derivatives of interpolation function          C  
C********************************************************************C
C     Programmed after Dhatt & Touzot p.64                           C
C********************************************************************C

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION VNI(INEL,1),VJ1(NDIM,1),VNIX(INEL,1)
      DATA ZERO/0.D0/

C-------------------------------------------------------
      DO 20 I=1,NDIM
        DO 20 J=1,INEL
          C=ZERO
          DO 10 IJ=1,NDIM
            C=C+VJ1(I,IJ)*VNI(J,IJ)
10        CONTINUE
          VNIX(J,I)=C
20    CONTINUE

      RETURN
      END
