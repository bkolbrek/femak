      SUBROUTINE BAXI(VNIX,VNI,INEL,R,VBE)
C**********************************************************************C
C BAXI.FOR                                                             C
C  10/10-1994   NTH                                                    C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  To form matrix B for axisymmeteric material type elastic solid.     C
C  Now used for biot material.                                         C
C  Used with 8-noded quadrilateral.                                    C
C    Input                                                             C
C      VNIX Derivatives of int.pol.func w.r.t. x,y,z                   C
C      VNI  Shape functions, for this Gauss point                      C
C      INEL Number of interpolation functions                          C
C      R    R-coordinate, used for axisymmetric version                C
C    Output                                                            C
C      VBE  Matrix B                                                   C
C**********************************************************************C
C  Programmed after Dhatt & Touzot, s.206                              C
C**********************************************************************C
      
      INTEGER INEL
      REAL*8 VNIX(INEL,1),VBE(4,1),VNI(1)
      REAL*8 R
      INTEGER I,J
      REAL*8 C1,C2,ZERO
      DATA ZERO/0.D0/
      
      J=1
      DO I=1,INEL
        C1=VNIX(I,1)
        C2=VNIX(I,2)
        VBE(1,J)=C1
        VBE(1,J+1)=ZERO
        VBE(2,J)=ZERO
        VBE(2,J+1)=C2
        VBE(3,J)=ZERO
        VBE(3,J+1)=VNI(I)/R
c        VBE(3,J+1)=zero
        VBE(4,J)=C2
        VBE(4,J+1)=C1
        J=J+2
      ENDDO

      RETURN
      END
