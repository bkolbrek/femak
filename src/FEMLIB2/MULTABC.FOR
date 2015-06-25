
      SUBROUTINE MULTABC(A,NA,MA,NRA,MRA,B,NB,MB,NRB,MRB,C,NRC,MRC,MPAR)
C**********************************************************************C
C     MULTABC.FOR						       C
C**********************************************************************C
C     5/11/1992							       C
C     Kompaniet Numerisk Lyd					       C
C     Tonni F. Johansen 					       C
C**********************************************************************C
C     multiplikasjon av to komplekse matriser			       C
C	 A    - matrix na x ma					       C
C	 B    - matrix nb x mb					       C
C	 C    - result of multiplication			       C
C	 NA,MA - dimension of matrix A				       C
C	 NRA,MRA - real dimension of matrix A			       C
C	 NB,MB - dimension of matrix B				       C
C	 NRB,MRB - real dimension of matrix B			       C
C	 NRC,MRC - real dimension of matrix C			       C
C	 MPAR  - Parameter beskriver multiplikasjon:		       C
C	   .EQ. 0  C=A*B					       C
C	   .EQ. 1  C=A'*B                  ('-transponert)	       C
C	   .EQ. 2  C=A*B'                                              C
C**********************************************************************C
      INTEGER NRA,MRA,NRB,MRB,NRC,MRC,NA,MA,NB,MB
      COMPLEX*16 A(NRA,MRA),B(NRB,MRB),C(NRC,MRC)

      INTEGER I,J,IJ

      COMPLEX*16 ZERO
      DATA ZERO/(0.D0,0.D0)/

      IF (MPAR.EQ.0) THEN
C..........C=A*B
	DO I=1,NA
	  DO J=1,MB
	    C(I,J)=ZERO
	    DO IJ=1,NB
	      C(I,J)=C(I,J)+A(I,IJ)*B(IJ,J)
	    ENDDO
	  ENDDO
	ENDDO
C..........Slutt C=A*B

      ELSEIF (MPAR.EQ.1) THEN
C..........C=A'*B

	DO I=1,MA
	  DO J=1,MB
	    C(I,J)=ZERO
	    DO IJ=1,NB
	      C(I,J)=C(I,J)+A(IJ,i)*B(IJ,j)
	    ENDDO
	  ENDDO
	ENDDO
C..........Slutt C=A'*B

      ELSEIF (MPAR.EQ.2) THEN
C..........C=A*B'

	DO I=1,NA
	  DO J=1,NB
	    C(I,J)=ZERO
	    DO IJ=1,MB
	      C(I,J)=C(I,J)+A(i,IJ)*B(J,IJ)
	    ENDDO
	  ENDDO
	ENDDO
C..........Slutt C=A*B'

      ENDIF

      RETURN
      END
