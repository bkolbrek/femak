C***********************************************************************
C CALCUL DE LA FONCTION DE BESSEL D'ORDRE ZERO D'UN ARGUMENT COMPLEXE
C
C***********************************************************************
      COMPLEX*16 FUNCTION BESS0(X)

      COMPLEX*16   X

      COMPLEX*16   E,MOINSJ,J0
      REAL*8       NF,PI,A
      INTEGER      N

      IF ( CDABS(X)-11. ) 2,3,3

C PETIT CABS(X) < 11.
 2    J0 = CMPLX(1.,0.)
      A = 1.

      DO 5 N=1,34
        CALL FAKULT(DBLE(N),NF)
        A = -A
        E = (A/(NF**2))*(X/2.)**(2.*N)
        J0 = J0 + E
        IF ( CDABS(E) - 0.00001 ) 4,4,5
 5    CONTINUE
      GOTO 4

C GRAND CABS(X) > 11.
 3    PI = 3.141592654
      MOINSJ = CMPLX(0.0,-1.0)
c      J0 = CSQRT(2./(PI*X)) * CCOS( X - PI/4.)
      J0 = 1. /SQRT(2.*PI) /CDSQRT(X) *
     1     (1 + CDEXP(MOINSJ *2*(X- PI/4.))) /CDEXP( MOINSJ *(X- PI/4.))


 4    CONTINUE
      BESS0 = J0

      END


C***********************************************************************
C CALCUL DE LA FONCTION DE BESSEL D'ORDRE UN D'UN ARGUMENT COMPLEXE
C
C***********************************************************************
      COMPLEX*16 FUNCTION BESS1(X)

      COMPLEX*16 X

      COMPLEX*16 E,MOINSJ,J1
      REAL*8     NF,PI,A,N1
      INTEGER    N



      IF (CDABS(X) - 13.) 2,3,3

C PETIT CABS(X) < 11.
 2    A=-1.
      J1=CMPLX(0.,0.)

      DO 5 N=1,34
        N1 = FLOAT(N) - 1.D00
        CALL FAKULT(N1,NF)
        A = -A
        E = A/(NF**2 *N) *(X/2.)**(2.*N1 + 1.)
        J1 = J1 + E
        IF (CDABS(E) - 0.00001) 4,4,5
 5    CONTINUE
      GOTO 4

C GRAND CABS(X) > 11.
 3    PI = 3.141592654
      MOINSJ = CMPLX(0.0,-1.0)
c      J1 = CDSQRT(2./(PI*X)) * CCOS( X - 3.*PI/4.)
      J1 = 1./SQRT(2.*PI)/CDSQRT(X)*(1 +CDEXP(MOINSJ *2*(X-3.*PI/4.)))/
     1                               CDEXP( MOINSJ *(X-3.*PI/4.))

 4    CONTINUE
      BESS1 = J1
      END

C***********************************************************************
C CALCUL DE LA FACTORIELLE DU NOMBRE NF = N!
C
C***********************************************************************
      SUBROUTINE FAKULT(N,NF)

      REAL*8 N,NF

      REAL*8 A,N1
      INTEGER L

      IF (N .LT. 0.D0) WRITE(*,*) 'FACTORIELLE D''UN NOMBRE NEGATIF !!!'
      IF (N - 1.) 8,8,9

C N = 1. OU BIEN N = 0.
 8    NF=1.
      GOTO 10

C N > 1.
 9    N1 = N - 1
      NF=N
      DO 11 L=1,N1
        A = N - FLOAT(L)
        NF = NF * A
 11   CONTINUE

 10   CONTINUE

      RETURN
      END
C***********************************************************************
