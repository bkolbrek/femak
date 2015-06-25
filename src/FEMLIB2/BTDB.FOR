      SUBROUTINE BTDB(VKE,VBE,VDE,IDLE,IMATD,NSYM)
C********************************************************************C
C BTDB.FOR                                                           C
C  19/3-1992                                                         C
C  Tonni F. Johansen  NTH                                            C
C********************************************************************C
C  To add the product B(T).D.B to VKE                                C
C    Input                                                           C
C      VKE   Element matrix nonsymmetrical (NSYM.EQ.1)               C
C                              symmetrical (NSYM.EQ.0)               C
C      VBE   Matrix B                                                C
C      VDE   Matrix D (Full)                                         C
C      IDLE  Total numbers of D.O.F. per element                     C
C      IMATD Dimension of matrix D (max.6)                           C
C    Output                                                          C
C      VKE                                                           C
C********************************************************************C
C  Programmert after Dhatt & Touzot, p.207                           C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION VKE(1),VBE(IMATD,1),VDE(IMATD,1),T(6)

      DATA ZERO/0.D0/

C-------------------------------------------------------
      IJ=1
      IMAX=IDLE
      DO 40 J=1,IDLE
        DO 20 I1=1,IMATD
          C=ZERO
          DO 10 J1=1,IMATD
10          C=C+VDE(I1,J1)*VBE(J1,J)
20        T(I1)=C
        IF (NSYM.EQ.0) IMAX=J
        DO 40 I=1,IMAX
          C=ZERO
          DO 30 J1=1,IMATD
30          C=C+VBE(J1,I)*T(J1)
          VKE(IJ)=VKE(IJ)+C
40      IJ=IJ+1

      RETURN
      END
