      SUBROUTINE ASSELC(IKG,IFG,IDLE,NSYM,KLOCE,KLD,VIE,VSRE,VIGS,
     1	VIGD,VIGI,VSRC)
C***************************************************************C
C     ASSELC.FOR						C
C***************************************************************C
C     29/7 1992 						C
C     Kompaniet Numerisk Lyd					C
C     Tonni F.Johansen						C
C     Programmert etter ASSEL (Dhatt og Touzot s.413)		C
C***************************************************************C
C     TO ASSEMBLE AN ELEMENT MATRIX AND/OR VECTOR		C
C     (MATRIX SYMMETRICAL OR NOT)				C
C	INPUT							C
C	   IKG	   IF IKG.EQ.1 ASSEMBLE ELEMENT MATRIX KE	C
C	   IFG	   IF IFG.EQ.1 ASSEMBLE ELEMENT VECTOR FE	C
C	   IDLE    ELEMENT NUMBER D.O.F.			C
C	   NSYM    0=SYMMETRIC PROBLEM, 1=UNSYMMETRIC PROBLEM	C
C	   KLOCE   ELEMENT LOCALIZATION VECTOR			C
C	   KLD	   CUMULATIVE COLUMN HEIGHTS OF KG		C
C	   VIE	   ELEMENT MATRIX AE (FULL OR UPPER TRIANGLE BY C
C		   DESCENDING COLUMNS)				C
C	   VSRE	   ELEMENT VECTOR SRE				C
C	OUTPUT							C
C	   VIGS,VIGD,VIGI    GLOBAL MATRIX (SKYLINES)		C
C		   (SYMMETRIC OR NOT)				C
C	   VSRC	   GLOBAL SOURCE VECTOR				C
C===============================================================C

      COMPLEX*16 VIE(1),VIGS(1),VIGD(1),VIGI(1),VSRE(1),VSRC(1)
      INTEGER KLOCE(1),KLD(1)
      INTEGER IKG,IFG,IDLE,NSYM

      INTEGER IQ,IEQ,IEQ0,IEQ1,JD,JL,I0,ID,IL,IJ,I
C---------------------------------------------------------------
C
C------  ASSEMBLE ELEMENT MATRIX
C
      IF(IKG.EQ.1) THEN

	IEQ0=IDLE
	IEQ1=1
C------  FOR EACH COLUMN OF KE
	DO 90 JD=1,IDLE
	  IF(NSYM.NE.1) IEQ0=JD
	  JL=KLOCE(JD)
	  IF(JL.GT.0) THEN
	    I0=KLD(JL+1)
	    IEQ=IEQ1
	    IQ=1
C------  FOR EACH ROW OF KE
	    DO 80 ID=1,IDLE
	      IL=KLOCE(ID)
	      IF(NSYM.NE.1) THEN
		IF(ID-JD.GE.0) THEN
		  IQ=ID
		ENDIF
	      ENDIF

	      IF(IL.GT.0) THEN
		IJ=JL-IL

		IF(IJ.EQ.0) THEN
C------  DIAGONAL TERMS OF KG
		  VIGD(IL)=VIGD(IL)+VIE(IEQ)
		ELSEIF (IJ.GT.0) THEN
C------  UPPER TRIANGLE TERMS OF KG
		  I=I0-IJ
		  VIGS(I)=VIGS(I)+VIE(IEQ)
		ELSEIF (IJ.LT.0) THEN
C------  LOWER TRIANGLE TERMS OF KG
		  IF(NSYM.EQ.1) THEN
		    I=KLD(IL+1)+IJ
		    VIGI(I)=VIGI(I)+VIE(IEQ)
		  ENDIF
		ENDIF
	      ENDIF

80	    IEQ=IEQ+IQ
	  ENDIF
90	IEQ1=IEQ1+IEQ0
C
C------  ASSEMBLE ELEMNT LOAD VECTOR
C
      ELSEIF(IFG.EQ.1) THEN
	DO 120 ID=1,IDLE
	  IL=KLOCE(ID)
	  IF(IL.GT.0) THEN
	    VSRC(IL)=VSRC(IL)+VSRE(ID)
	  ENDIF
120	CONTINUE
      ENDIF

      RETURN
      END  
