      SUBROUTINE INVERSC(VP,N,IVP,K,DET)
C**********************************************************************C
C     INVERSC.FOR                                                      C
C**********************************************************************C
C     4/11/1992                                                        C
C     Kompaniet Numerisk Lyd                                           C
C     Tonni F. Johansen                                                C
C     Programmert etter Dhatt og Touzot s.57                           C
C**********************************************************************C
C     TO INVERT A NON-SYMMETRIC COMPLEX MATRIX WITH SEARCH OF A        C
C     NON-ZERO PIVOT IN A COLUMN                                       C
C       INPUT                                                          C
C          VP      MATRIX TO BE INVERTED                               C
C          N       ORDER OF THE MATRIX                                 C
C          IVP     DIMENSION OF THE MATRIX IN THE CALLING PROGRAM      C
C          K       INTEGER WORKING ARRAY WITH LENGTH N                 C
C       OUTPUT                                                         C
C          VP      INVERSE MATRIX                                      C
C          DET     DETERMINANT                                         C
C======================================================================C
C      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IVP,N
      INTEGER K(N)
      COMPLEX*16 VP(IVP,IVP),DET

      COMPLEX*16 PIV,C,UN,ZERO
      INTEGER I,II,I1,J,J1,JJ
      REAL*8 EPS


      DATA ZERO/(0.D0,0.D0)/,UN/(1.D0,0.D0)/,EPS/1.D-13/
c      ABS(X)=DABS(X)
C----------------------------------------------------------------------
      DET=UN
      DO I=1,N
	K(I)=I
      ENDDO

C------  START INVERSION
      DO II=1,N
C------  SEARCH FOR NON ZERO PIVOT IN COLUMN II
	DO I=II,N
	  PIV=VP(I,II)
	  IF(ABS(PIV).GT.EPS) GO TO 20
	ENDDO
	DET=ZERO
	RETURN
C------  EXCHANGE LINES II AND I
20      DET=DET*PIV
	IF(I.NE.II) THEN
	  I1=K(II)
	  K(II)=K(I)
	  K(I)=I1
	  DO J=1,N
	    C=VP(I,J)
	    VP(I,J)=VP(II,J)
	    VP(II,J)=C
	  ENDDO
	  DET=-DET
C------  NORMALIZE PIVOT LINE
	ELSE
	  C=UN/PIV
	  VP(II,II)=UN
	  DO J=1,N
	    VP(II,J)=VP(II,J)*C
	  ENDDO
	ENDIF
C------  ELIMINATION
	DO I=1,N
	  IF(I.NE.II) THEN
	    C=VP(I,II)
	    VP(I,II)=ZERO
	    DO J=1,N
	      VP(I,J)=VP(I,J)-C*VP(II,J)
	    ENDDO
	  ENDIF
	ENDDO


      ENDDO


C------  REORDER THE COLUMNS OF INVERSE MATRIX
      DO J=1,N
C------  FIND J1 SUCH THAT K(J1)=J
	DO J1=J,N
	  JJ=K(J1)
	  IF(JJ.EQ.J) GO TO 100
	ENDDO

100     IF(J.NE.J1) THEN
C------  EXCHANGE COLUMNS J AND J1
	  K(J1)=K(J)
	  DO I=1,N
	    C=VP(I,J)
	    VP(I,J)=VP(I,J1)
	    VP(I,J1)=C
	  ENDDO
	ENDIF
      ENDDO

      RETURN
      END
