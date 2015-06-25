      SUBROUTINE B08(VNIX,INEL,VBBE)
C**********************************************************************C
C B08.FOR                                                              C
C  23/1-1995  Kompaniet Numerisk lyd                                   C
C  Marc Dhainaut                                                       C                                              
C**********************************************************************C
C  To form matrix B (2 dimensional elasticity)                         C
C    Input                                                             C
C      VNIX Derivatives of int.pol.func w.r.t. x,y,z                   C
C      INEL Number of interpolation functions                          C
C    Output                                                            C
C      VBBE Matrix B                                                   C
C**********************************************************************C
C Programmert etter Dhatt & Touzot, s.194-197                          C 
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION VNIX(INEL,1),VBBE(3,1)
      REAL*8 C1,C2
      
      DATA ZERO/0.D0/


      J=1
      DO 10 I=1,INEL
	C1=VNIX(I,1)
	C2=VNIX(I,2)
	VBBE(1,J)  = ZERO
	VBBE(1,J+1)= ZERO
	VBBE(1,J+2)= C1
	VBBE(2,J)  = ZERO
	VBBE(2,J+1)=-C2
	VBBE(2,J+2)= ZERO
	VBBE(3,J)  = ZERO
	VBBE(3,J+1)=-C1
	VBBE(3,J+2)= C2
10	J=J+3

      RETURN
      END
