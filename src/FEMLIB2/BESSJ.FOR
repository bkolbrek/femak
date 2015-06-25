      DOUBLE PRECISION FUNCTION BESSJ0(X)
C********************************************************************C
C  BESSJ.FOR                                                         C
C  2/10-1992                                                         C
C  Tonni F. Johansen  NTH                                            C
C********************************************************************C
C  Returns besselfunction J0 for any X.                              C
C  Copied from Numerical recepies p.172-173                          C
C********************************************************************C

      REAL*8 X
      REAL*8 AX,Z,XX
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6
      REAL*8 S1,S2,S3,S4,S5,S6

      DATA R1,R2,R3,R4,R5,R6/57568490574.D0,-13362590354.D0,
     1	    651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0/
      DATA S1,S2,S3,S4,S5,S6/57568490411.D0,1029532985.D0,
     1	   9494680.718D0,59272.64853D0,267.8532712D0,1.D0/
      DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,
     1	  -.2073370639D-5,.2093887211D-6/
      DATA Q1,Q2,Q3,Q4,Q5/-.1562499995D-1,.1430488765D-3,-.6911147651D-5
     1	  ,.7621095161D-6,-.934945152D-7/



      IF (ABS(X).LT.8) THEN
	Y=X*X
	BESSJ0=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))/
     1	       (S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
	AX=ABS(X)
	Z=8./AX
	Y=Z*Z
	XX=AX-.785398164
	BESSJ0=SQRT(.636619772/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*
     1		    P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
      ENDIF


      RETURN

      END




      DOUBLE PRECISION FUNCTION BESSJ1(X)
C********************************************************************C
C  BESSJ.FOR                                                         C
C  2/10-1992                                                         C
C  Tonni F. Johansen NTH                                             C
C********************************************************************C
C  Returns besselfunction J1 for any X.                              C
C  Copied from Numerical recepies p.173-174                          C
C********************************************************************C

      REAL*8 X
      REAL*8 AX,Z,XX
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6
      REAL*8 S1,S2,S3,S4,S5,S6

      DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,242396853.1D0
     1	  ,-2972611.439D0,15704.48260D0,-30.16036606D0/
      DATA S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0,
     1	   18583304.74D0,99447.43394D0,376.9991397D0,1.D0/
      DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,.2457520174D-5
     1	  ,-.240337019D-6/,Q1,Q2,Q3,Q4,Q5/.04687499995D0,
     2	   -.2002690873D-3,.8449199096D-5,-.88228987D-6,.105787412D-6/



      IF (ABS(X).LT.8) THEN
	Y=X*X
	BESSJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))/
     1	       (S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
	AX=ABS(X)
	Z=8./AX
	Y=Z*Z
	XX=AX-2.3561194491
	BESSJ1=SQRT(.636619772/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*
     1		    P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
     2		    *SIGN(1.0D0,X)
      ENDIF

      RETURN
      END
