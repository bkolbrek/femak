      SUBROUTINE KAGINT(N,EMN,K,VI100KM)
C**********************************************************************C
C     KAGINT.FOR						       C
C     Tonni F. Johansen 					       C
C     Kompaniet Numerisk Lyd					       C
C     16/9/93							       C
C								       C
C**********************************************************************C
C     To assemble matrix E and T as described by Kagawa et.al.	       C
C     J.S.&V. vol.69  1980 p.207-228.				       C
C     Infinite integral method described in numerical recepies.        C
C**********************************************************************C
C   input:							       C
C      N  - Number of modes					       C
C      K  - wavenumber                                                 C
C   output:							       C
C      EMN  - Value of the integral				       C
C**********************************************************************C
      EXTERNAL INTGRR,INTGRI,BESSJ1,MIDINF,MIDSQL,MIDSQU
      DOUBLE PRECISION BESSJ1
      DOUBLE PRECISION INTGRR,INTGRI

      INTEGER N
      COMPLEX*16 EMN(N,N),VI100KM(1)
      REAL*8 K

      INTEGER IFL,IGIMP,IFREQ,IGSRC,IFL1,NI100NOD,NI100MOD,NI100FL
      REAL*8 I100RAD
      COMMON/IMPEL/IFL,IGIMP,IFREQ,IGSRC,IFL1,
     1		   NI100NOD,NI100MOD,NI100FL,I100RAD

      INTEGER M,MR,MP,M1,M2,M3
      COMMON/ES/M,MR,MP,M1,M2,M3

C.......Local COMMON-block for this file, improper integral
      REAL*8 KA
      INTEGER MODM,MODN
      COMMON/IMPINT/MODM,MODN,KA

      REAL*8 VREAL,VIMAG,DUMS

      REAL*8 ZERO,HIGHNUM
      DATA ZERO/0.D0/,HIGHNUM/1.D30/

C----------------Beginning of routine
      KA=K*I100RAD

      DO MODM=1,N
	DO MODN=MODM,N
          VREAL=ZERO
          VIMAG=ZERO

C.......Integrate from 0 to KA (real part)
          CALL QROMO(INTGRR,0.D0,KA,DUMS,MIDSQU)
          VREAL=DUMS
C.......Integrate from KA to 2KA (imaginary part)
          CALL QROMO(INTGRI,KA,2.*KA,DUMS,MIDSQL)
          VIMAG=DUMS
C.............Integrate from 2ka to infinity (imaginary part)
          CALL QROMO(INTGRI,2.*KA,HIGHNUM,DUMS,MIDINF)
          VIMAG=VIMAG+DUMS
C

	  EMN(MODM,MODN)=VI100KM(MODN)*2*I100RAD*DCMPLX(VREAL,VIMAG)
	  IF (MODM.NE.MODN) THEN
	    EMN(MODN,MODM)=VI100KM(MODM)*2*I100RAD*DCMPLX(VREAL,VIMAG)
	  ENDIF
        ENDDO
      ENDDO



      RETURN
      END

      SUBROUTINE QROMO(FUNC,A,B,SS,CHOOSE)
C**********************************************************************C
C     KAGINT.FOR						       C
C     Tonni F. Johansen 					       C
C     Kompaniet Numerisk Lyd					       C
C     16/9/93							       C
C								       C
C**********************************************************************C
C     To integrate improper integrals.                                 C
C     Infinite integral method described in numerical recepies.        C
C**********************************************************************C
C   input:							       C
C      FUNC - Name of function to evaluate integrand                   C
C      A    - Lower limit of integration inteval                       C
C      B    - Upper limit of integration inteval                       C
C      CHOOSE - Name of integration rule subroutine                    C
C**********************************************************************C
C   output:							       C
C      SS  - Value of the integral				       C
C**********************************************************************C

      EXTERNAL FUNC, CHOOSE
      DOUBLE PRECISION FUNC

      REAL*8 A,B,SS,DSS

      PARAMETER (EPS=1D-3,JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
      REAL*8 S(JMAXP),H(JMAXP)
C
      REAL*8 KA
      INTEGER MODM,MODN
      COMMON/IMPINT/MODM,MODN,KA

      H(1)=1
      DO 11 J=1,JMAX
        CALL CHOOSE(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL POLINT(H(J-KM),S(J-KM),K,0.D0,SS,DSS)
          IF (ABS(DSS).LT.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=H(J)/9.
11    CONTINUE
      PAUSE 'TOO MANY STEPS.'

      RETURN
      END



      SUBROUTINE MIDINF(FUNK,AA,BB,S,N)
C**********************************************************************C
C     KAGINT.FOR						       C
C     Tonni F. Johansen 					       C
C     Kompaniet Numerisk Lyd					       C
C     16/9/93							       C
C								       C
C**********************************************************************C
C     To integrate improper integrals, infinite domain                 C
C     Infinite integral method described in numerical recepies, p.118  C
C**********************************************************************C
C   input:							       C
C      FUNC - Name of function to evaluate integrand                   C
C      A    - Lower limit of integration inteval                       C
C      B    - Upper limit of integration inteval                       C
C      N    - Iteration number                                         C
C**********************************************************************C
C   output:							       C
C      S   - Value of the integral				       C
C**********************************************************************C

      EXTERNAL FUNK
      DOUBLE PRECISION FUNC,FUNK
C
      REAL*8 AA,BB,S
      INTEGER N
C
      REAL*8 A,B

      REAL*8 KA, X
      INTEGER MODM,MODN
      COMMON/IMPINT/MODM,MODN,KA


      FUNC(X)=FUNK(1.D0/X)/X**2
C
      B=1./AA
      A=1./BB
      IF (N.EQ.1) THEN
        S=(B-A)*FUNC(0.5D0*(A+B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/(3.*TNM)
        DDEL=DEL+DEL
        X=A+0.5*DEL
        SUM=0.0D0
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DDEL
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=(S+(B-A)*SUM/TNM)/3.
        IT=3*IT
      ENDIF
C
      RETURN
      END

      SUBROUTINE MIDSQL(FUNK,AA,BB,S,N)
C**********************************************************************C
C     KAGINT.FOR						       C
C     Tonni F. Johansen 					       C
C     Kompaniet Numerisk Lyd					       C
C     16/9/93							       C
C								       C
C**********************************************************************C
C     To integrate improper integrals, singularity at lower limit      C
C     Infinite integral method described in numerical recepies.p.120   C
C**********************************************************************C
C   input:							       C
C      FUNC - Name of function to evaluate integrand                   C
C      A    - Lower limit of integration inteval                       C
C      B    - Upper limit of integration inteval                       C
C      N    - Iteration number                                         C
C**********************************************************************C
C   output:							       C
C      S   - Value of the integral				       C
C**********************************************************************C

      EXTERNAL FUNK
      DOUBLE PRECISION FUNC,FUNK
C
      REAL*8 AA,BB,S
      INTEGER N
C
      REAL*8 A,B

      REAL*8 KA, X
      INTEGER MODM,MODN
      COMMON/IMPINT/MODM,MODN,KA


      FUNC(X)=2.*X*FUNK(AA+X**2)
      B=SQRT(BB-AA)
      A=0.
      IF (N.EQ.1) THEN
        S=(B-A)*FUNC(0.5D0*(A+B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/(3.*TNM)
        DDEL=DEL+DEL
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DDEL
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=(S+(B-A)*SUM/TNM)/3.
        IT=3*IT
      ENDIF
C
      RETURN
      END


      SUBROUTINE MIDSQU(FUNK,AA,BB,S,N)
C**********************************************************************C
C     KAGINT.FOR						       C
C     Tonni F. Johansen 					       C
C     Kompaniet Numerisk Lyd					       C
C     16/9/93							       C
C								       C
C**********************************************************************C
C     To integrate improper integrals, singularity at upper limit      C
C     Infinite integral method described in numerical recepies.p.120   C
C**********************************************************************C
C   input:							       C
C      FUNC - Name of function to evaluate integrand                   C
C      A    - Lower limit of integration inteval                       C
C      B    - Upper limit of integration inteval                       C
C      N    - Iteration number                                         C
C**********************************************************************C
C   output:							       C
C      S   - Value of the integral				       C
C**********************************************************************C
      EXTERNAL FUNK
      DOUBLE PRECISION FUNC,FUNK
C
      REAL*8 AA,BB,S
      INTEGER N
C
      REAL*8 A,B

      REAL*8 KA, X
      INTEGER MODM,MODN
      COMMON/IMPINT/MODM,MODN,KA


      FUNC(X)=2.*X*FUNK(BB-X**2)
      B=SQRT(BB-AA)
      A=0.0D0
      IF (N.EQ.1) THEN
        S=(B-A)*FUNC(0.5*(A+B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/(3.*TNM)
        DDEL=DEL+DEL
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DDEL
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=(S+(B-A)*SUM/TNM)/3.
        IT=3*IT
      ENDIF
C
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION INTGRR(X)
C**********************************************************************C
C     KAGINT.FOR						       C
C     Tonni F. Johansen 					       C
C     Kompaniet Numerisk Lyd					       C
C     16/9/93							       C
C								       C
C**********************************************************************C
C     Function used to evaluate real part of integrand as described    C
C     by Kagawa et.al. J.S.&V. vol.69  1980 p.207-228.                 C
C**********************************************************************C
C   input:							       C
C      X - Argument of the function                                    C
C**********************************************************************C
      REAL*8 X
C
      DOUBLE PRECISION BESSJ1
      REAL*8 KA
      INTEGER MODM,MODN
      COMMON/IMPINT/MODM,MODN,KA

      REAL*8 GAMMA(10)
C...........Zero in first order bessel function Gamma
      DATA GAMMA/0.0,3.8317059702,7.0155866698,10.1734681351,
     1	      13.3236919363,16.4706300509,19.6158585105,
     2	      22.7600843806,25.9036720876,29.0468285349/


      BJ=BESSJ1(X)
      INTGRR=X*X*X*BJ*BJ/(SQRT(KA*KA-X*X)*
     1    (X*X-GAMMA(MODM)*GAMMA(MODM))*(X*X-GAMMA(MODN)*GAMMA(MODN)))

      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION INTGRI(X)
C**********************************************************************C
C     KAGINT.FOR						       C
C     Tonni F. Johansen 					       C
C     Kompaniet Numerisk Lyd					       C
C     16/9/93							       C
C								       C
C**********************************************************************C
C     Function used to evaluate imaginary part of integrand as         C
C     describedby Kagawa et.al. J.S.&V. vol.69  1980 p.207-228.        C
C**********************************************************************C
C   input:							       C
C      X - Argument of the function                                    C
C**********************************************************************C
      REAL*8 X
      DOUBLE PRECISION BESSJ1
C
      REAL*8 KA
      INTEGER MODM,MODN
      COMMON/IMPINT/MODM,MODN,KA

      REAL*8 GAMMA(10)
C...........Zero in first order bessel function Gamma
      DATA GAMMA/0.0,3.8317059702,7.0155866698,10.1734681351,
     1	      13.3236919363,16.4706300509,19.6158585105,
     2	      22.7600843806,25.9036720876,29.0468285349/


      BJ=BESSJ1(X)
      INTGRI=X*X*X*BJ*BJ/(SQRT(X*X-KA*KA)*
     1    (X*X-GAMMA(MODM)*GAMMA(MODM))*(X*X-GAMMA(MODN)*GAMMA(MODN)))

      RETURN
      END
C
C
C



C
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
C**********************************************************************C
C     KAGINT.FOR						       C
C     Tonni F. Johansen 					       C
C     Kompaniet Numerisk Lyd					       C
C     16/9/93							       C
C								       C
C**********************************************************************C
C     Interpolation subroutine for polynomial interpolation            C
C     described in numerical recepies.p.82                             C
C**********************************************************************C
C   input:							       C
C      XA,YA - Arrays with x and y values                              C
C      N     - Length of XA,YA arrays                                  C
C      X     - X argument of estimate                                  C
C**********************************************************************C
C   output:							       C
C      Y   - Estimate                                                  C
C      DY  - estimate of error in estimate (!!)                        C
C**********************************************************************C

      PARAMETER (NMAX=10)
      INTEGER N
      REAL*8 XA(N),YA(N),C(NMAX),D(NMAX)
      REAL*8 DY,X,Y
      REAL*8 DIF,DIFT
C
      INTEGER NS,I
C
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1

      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF (DEN.EQ.0.) PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M) THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
