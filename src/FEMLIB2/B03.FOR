      SUBROUTINE B03(VNIX,INEL,VBE)
C**********************************************************************C
C B03.FOR                                                              C
C  29/5-1996              NTH                                          C
C  Tonni F. Johansen  & Bruno Brouard                                  C
C**********************************************************************C
C  To form matrix B (3 dimensional elasticity)                         C
C    Input                                                             C
C      VNIX   Derivatives of int.pol.func w.r.t. x,y,z                 C
C      INEL   Number of interpolation functions                        C
C    Output                                                            C
C      VBE    Matrix B                                                 C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION VNIX(INEL,1),VBE(6,1)

      DATA ZERO/0.D0/


      J=1
      DO I=1,INEL
        C1=VNIX(I,1)
        C2=VNIX(I,2)
        C3=VNIX(I,2)
        VBE(1,J)=C1
        VBE(1,J+1)=ZERO
        VBE(1,J+2)=ZERO
        VBE(2,J)=ZERO
        VBE(2,J+1)=C2
        VBE(2,J+2)=ZERO
        VBE(3,J)=ZERO
        VBE(3,J+1)=ZERO
        VBE(3,J+2)=C3
        VBE(4,J)=C1
        VBE(4,J+1)=C2
        VBE(4,J+2)=ZERO
        VBE(5,J)=ZERO
        VBE(5,J+1)=C2
        VBE(5,J+2)=C3
        VBE(6,J)=C1
        VBE(6,J+1)=ZERO
        VBE(6,J+2)=C3
        J=J+3
      ENDDO

      RETURN
      END


