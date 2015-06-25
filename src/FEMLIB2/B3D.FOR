      SUBROUTINE B3D(VNIX,INEL,VBE)
C**********************************************************************C
C B3D.FOR                                                              C
C  21/2-1995   NTH                                                     C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  To form matrix B (3 dimensional elasticity), used for biot type     C
C  material and elastic solid.                                         C
C    Input                                                             C
C      VNIX Derivatives of int.pol.func w.r.t. x,y,z                   C
C      INEL Number of interpolation functions                          C
C    Output                                                            C
C      VBE  Matrix B                                                   C
C**********************************************************************C
C  Programmed after Dhatt & Touzot, s.206                              C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER INEL
      REAL*8 VNIX(INEL,1),VBE(6,1)
      
      REAL*8 C1,C2,C3,ZERO
      INTEGER I,J
      DATA ZERO/0.D0/

      
      J=1
      DO I=1,INEL
        C1=VNIX(I,1)
        C2=VNIX(I,2)
        C3=VNIX(I,3)
        VBE(1,J)=C1
        VBE(1,J+1)=ZERO
        VBE(1,J+2)=ZERO

        VBE(2,J)=ZERO
        VBE(2,J+1)=C2
        VBE(2,J+2)=ZERO
          
        VBE(3,J)=ZERO
        VBE(3,J+1)=ZERO
        VBE(3,J+2)=C3
          
        VBE(4,J)=C2
        VBE(4,J+1)=C1
        VBE(4,J+2)=ZERO

        VBE(5,J)=ZERO
        VBE(5,J+1)=C3
        VBE(5,J+2)=C2

        VBE(6,J)=C3
        VBE(6,J+1)=ZERO
        VBE(6,J+2)=C1

        J=J+3
      ENDDO

      RETURN
      END
