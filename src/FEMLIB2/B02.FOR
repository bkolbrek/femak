      SUBROUTINE B02(VNIX,INEL,VBE)
C**********************************************************************C
C B02.FOR							       C
C  19/3-1992	  Kompaniet Numerisk lyd			       C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  To form matrix B (2 dimensional elasticity)			       C
C    Input                                                             C
C	 VNIX	Derivatives of int.pol.func w.r.t. x,y,z	       C
C	 INEL	Number of interpolation functions		       C
C    Output                                                            C
C	 VBE	Matrix B					       C
C**********************************************************************C
C  Programmert etter Dhatt & Touzot, s.206			       C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION VNIX(INEL,1),VBE(3,1)

      DATA ZERO/0.D0/


      J=1
      DO I=1,INEL
        C1=VNIX(I,1)
        C2=VNIX(I,2)
        VBE(1,J)=C1
        VBE(1,J+1)=ZERO
        VBE(2,J)=ZERO
        VBE(2,J+1)=C2
        VBE(3,J)=C2
        VBE(3,J+1)=C1
        J=J+2
      ENDDO

      RETURN
      END
