      SUBROUTINE BTB08(VKE,VBSE,COEF,IDLE,NDIM,NSYM)
C**********************************************************************C
C BTB08.FOR                                                            C
C  23/1-1995  Kompaniet Numerisk lyd                                   C
C  Marc Dhainaut                                                       C
C**********************************************************************C
C  To add the product B(T).B to VKE                                    C
C    Input                                                             C
C      VKE   Element matrix nonsymmetrical (NSYM.EQ.1)                 C
C                         symmetrical (NSYM.EQ.0)                      C
C      VBSE  Matrix B                                                  C
C      IDLE  Total numbers of D.O.F. per element                       C
C      IMATD Dimension of matrix D (max.6)                             C
C   Output                                                             C
C      VKE                                                             C
C**********************************************************************C
C  Programmert etter Dhatt & Touzot, s.194-197                         C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION VKE(1),VBSE(NDIM,1)
      REAL*8 C
      
      DATA ZERO/0.D0/

C-------------------------------------------------------
      IJ=1
      IMAX=IDLE
      DO 40 J=1,IDLE
	IF (NSYM.EQ.0) IMAX=J
	DO 40 I=1,IMAX
	  C=ZERO
	  DO 30 J1=1,NDIM
30	    C=C+VBSE(J1,I)*VBSE(J1,J)
	  VKE(IJ)=VKE(IJ)+COEF*C
40	IJ=IJ+1

      RETURN
      END
