      SUBROUTINE BB08(VNI,VNIX,INEL,VBSE)
C**********************************************************************C
C BB08.FOR                                                           C
C  23/1-1995  Kompaniet Numerisk lyd                                   C
C  Marc Dhainaut                                                       C
C**********************************************************************C
C  To form matrix B (2 dimensional elasticity)                         C
C    Input                                                             C
C      VNI  interpolation func.                                        C
C      VNIX Derivatives of int.pol.func w.r.t. x,y,z                   C
C      INEL Number of interpolation functions                          C
C    Output                                                            C
C      VBSE Matrix B                                                   C
C**********************************************************************C
C  Programmert etter Dhatt & Touzot, s.194-197                         C  
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION VNI(1),VNIX(INEL,1),VBSE(2,1)
      REAL*8 C0,C1,C2
      
      DATA ZERO/0.D0/


      J=1
      DO 10 I=1,INEL
      C0=VNI(I)
	C1=VNIX(I,1)
	C2=VNIX(I,2)
	VBSE(1,J)  = C1
	VBSE(1,J+1)= ZERO
	VBSE(1,J+2)= C0
	VBSE(2,J)  = C2
	VBSE(2,J+1)=-C0
	VBSE(2,J+2)= ZERO
10	J=J+3

      RETURN
      END
