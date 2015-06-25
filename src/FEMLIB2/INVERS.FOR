      SUBROUTINE INVERS(VP,N,IVP,K,DET)
C**********************************************************************C
C     INVERS.FOR						       C
C**********************************************************************C
C     28/2/1992 						       C
C     Kompaniet Numerisk Lyd					       C
C     Ulf R. Kristiansen					       C
C     Programmert etter Dhatt og Touzot s.57			       C
C**********************************************************************C
C     TO INVERT A NON-SYMMETRIC MATRIX WITH SEARCH OF A 	       C
C     NON-ZERO PIVOT IN A COLUMN				       C
C	INPUT							       C
C	   VP	   MATRIX TO BE INVERTED			       C
C	   N	   ORDER OF THE MATRIX				       C
C	   IVP	   DIMENSION OF THE MATRIX IN THE CALLING PROGRAM      C
C	   K	   INTEGER WORKING ARRAY WITH LENGTH N		       C
C	OUTPUT							       C
C	   VP	   INVERSE MATRIX				       C
C	   DET	   DETERMINANT					       C
C======================================================================C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VP(IVP,IVP),K(N)
      DATA ZERO/0.D0/,UN/1.D0/,EPS/1.D-13/
      ABS(X)=DABS(X)
C----------------------------------------------------------------------
      DET=UN
      DO 5 I=1,N
5     K(I)=I
C------  START INVERSION
      DO 80  II=1,N
C------  SEARCH FOR NON ZERO PIVOT IN COLUMN II
      DO 10 I=II,N
      PIV=VP(I,II)
      IF(ABS(PIV).GT.EPS) GO TO 20
10    CONTINUE
      DET=ZERO
      RETURN
C------  EXCHANGE LINES II AND I
20    DET=DET*PIV
      IF(I.EQ.II) GO TO 40
      I1=K(II)
      K(II)=K(I)
      K(I)=I1
      DO 30 J=1,N
      C=VP(I,J)
      VP(I,J)=VP(II,J)
30    VP(II,J)=C
      DET=-DET
C------  NORMALIZE PIVOT LINE
40    C=UN/PIV
      VP(II,II)=UN
      DO 50 J=1,N
50    VP(II,J)=VP(II,J)*C
C------  ELIMINATION
      DO 70 I=1,N
      IF(I.EQ.II) GO TO 70
      C=VP(I,II)
      VP(I,II)=ZERO
      DO 60 J=1,N
60    VP(I,J)=VP(I,J)-C*VP(II,J)
70    CONTINUE
80    CONTINUE
C------  REORDER THE COLUMNS OF INVERSE MATRIX
      DO 120 J=1,N
C------  FIND J1 SUCH THAT K(J1)=J
      DO 90 J1=J,N
      JJ=K(J1)
      IF(JJ.EQ.J) GO TO 100
90    CONTINUE
100   IF(J.EQ.J1) GO TO 120
C------  EXCHANGE COLUMNS J AND J1
      K(J1)=K(J)
      DO 110 I=1,N
      C=VP(I,J)
      VP(I,J)=VP(I,J1)
110   VP(I,J1)=C
120   CONTINUE
      RETURN
      END
