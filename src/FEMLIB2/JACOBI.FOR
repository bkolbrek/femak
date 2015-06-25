      SUBROUTINE JACOBI(VK,VM,N,NCYM,EPS,VALP0,VALP,VECT) 
C**********************************************************
C     JACOBI.FOR
C**********************************************************
C     13/3 1992
C     Kompaniet Numerisk Lyd
C     Programmert etter Dhatt og Touzot s.345
C     Tonni F. Johansen
C**********************************************************
C     To solve eigenproblems k-kambda, m by the generalized
C     jacobi method 
C     INPUT
C          VK       MATRIX K (UPPER TRIANGLE BY DESCENDING COLOUMNS)
C          VM       MATRIX M (UPPER TRIANGLE BY DESCENDING COLOUMNS)
C           N       ORDER OF MATRIX K AND M
C        NCYM       MAXIMUM ORDER OF SWEEPS ALLOWED  (15)
C         EPS       CONVERGENCE TOLERANCE    (1.D-12)
C     WORKSPACE
C       VALP0       WORKING VECTOR  (DIMENSION N)
C     OUTPUT
C        VALP       EIGENVALUES
C        VECT       EIGENVECTORS
C============================================================
      IMPLICIT REAL*8(A-H,O-Z)
c      COMMON/ES/M,MR,MP,M1,M2

      include 'femak.cmn'

      DIMENSION VK(1),VM(1),VALP0(N),VALP(N),VECT(N,N)

      DATA EPSD0/1.D-4/,ZERO/0.D0/,UN/1.D0/,DEUX/2.D0/,QUATR/4.D0/

      SQRT(X)=DSQRT(X)
      ABS(X)=DABS(X)
      EPS2=EPS*EPS
      ITR=0

C---------------VERIFY IF DIAGONAL TERMS ARE POSITIVE
      II=0
      DO 20 I=1,N
	II=II+I
        IF (VK(II).GT.ZERO.AND.VM(II).GT.ZERO) GO TO 10
        WRITE(MP,2000)I
2000    FORMAT('** ERROR, NEGATIVE DIAGONAL TERM IN JACOBI, ROW ',
     1          I5)
        STOP
10      VALP(I)=VK(II)/VM(II)
20      VALP0(I)=VALP(I)

C-------------INITIALIZE EIGENVECTORS
      
      DO 40 I=1,N
        DO 30 J=1,N
30        VECT(I,J)=ZERO
40      VECT(I,I)=UN

C-------------FOR EACH SWEEP
      DO 250 IC=1,NCYM
C------------DYNAMIC TOLERANCE
        EPSD=EPSD0**IC
C-------------SWEEP ROWWISE OVER UPPTER TRIANGLE
        IMAX=N-1
        II=0
        DO 180 I=1,IMAX
          I0=II+1
          II=II+I
          IP1=I+1
          IJ=II+I
          JJ=II
          DO 180 J=IP1,N
	    JP1=J+1
            JM1=J-1
            J0=JJ+1
            JJ=JJ+J
            J3=JJ-1

C----------------COMPUTE THE COUPLING FACTOR
            FK=(VK(IJ)*VK(IJ))/(VK(II)*VK(JJ))
            FM=(VM(IJ)*VM(IJ))/(VM(II)*VM(JJ))
            IF (FK.LT.EPSD.AND.FM.LT.EPSD) GO TO 180

C---------------COMPUTE THE TRANSFORMATION COEFFICIENTS  
            ITR=ITR+1
            C1=VK(II)*VM(IJ)-VM(II)*VK(IJ)
            C2=VK(JJ)*VM(IJ)-VM(JJ)*VK(IJ)
            C3=VK(II)*VM(JJ)-VM(II)*VK(JJ)
            DET=(C3*C3/QUATR)+(C1*C2)
            IF (DET.GE.ZERO) GO TO 50
            WRITE(MP,2005)I,J
2005        FORMAT('** ERROR, SINGULAR JACOBI TRANSFORMATION I=',I5,
     1              ' J=',I5)
            STOP

50          DET=SQRT(DET)
            D1=C3/DEUX+DET
            D2=C3/DEUX-DET
            D=D1
            IF (ABS(D2).GT.ABS(D1)) D=D2
            IF (D.EQ.ZERO) GO TO 60
            A=C2/D
            B=-C1/D
            GO TO 65
60          A=ZERO
            B=-VK(IJ)/VK(JJ)

C--------------------------MODIFY COLOUMNS OF K AND M
65          IF (I.EQ.1) GO TO 80
            IK=I0
            J1=IJ-1
            DO 70 JK=J0,J1
              C1=VK(IK)
              C2=VK(JK)
              VK(IK)=C1+B*C2
              VK(JK)=C2+A*C1
              C1=VM(IK)
              C2=VM(JK)
              VM(IK)=C1+B*C2
              VM(JK)=C2+A*C1
70            IK=IK+1

80          IF (I.EQ.JM1) GO TO 100
            IK=II+I
            J2=IJ+1
            IM=I
            DO 90 JK=J2,J3
              C1=VK(IK)
              C2=VK(JK)
              VK(IK)=C1+B*C2
              VK(JK)=C2+A*C1
              C1=VM(IK)
              C2=VM(JK)
              VM(IK)=C1+B*C2
              VM(JK)=C2+A*C1
              IM=IM+1
90            IK=IK+IM
100         IF (J.EQ.N) GO TO 120
            IK=IJ+J
            JK=JJ+J
            IM=J
            DO 110 JJK=JP1,N
              C1=VK(IK)
              C2=VK(JK)
              VK(IK)=C1+B*C2
              VK(JK)=C2+A*C1
              C1=VM(IK)
              C2=VM(JK)
              VM(IK)=C1+B*C2
              VM(JK)=C2+A*C1
              IM=IM+1
              IK=IK+IM
110           JK=JK+IM

120         C1=VK(II)
            C2=VK(IJ)
            C3=VK(JJ)
            B2=B*B
            BB=DEUX*B
            A2=A*A
            AA=DEUX*A
            VK(II)=C1+BB*C2+B2*C3
            VK(IJ)=ZERO
            VK(JJ)=C3+AA*C2+A2*C1
            C1=VM(II)
            C2=VM(IJ)
            C3=VM(JJ)
            VM(II)=C1+BB*C2+B2*C3
            VM(IJ)=ZERO
            VM(JJ)=C3+AA*C2+A2*C1

C------------------UPDATE EIGENVECTORS
            DO 170 IJ1=1,N
              C1=VECT(IJ1,I)
              C2=VECT(IJ1,J)
              VECT(IJ1,I)=C1+B*C2
170	      VECT(IJ1,J)=C2+A*C1
180         IJ=IJ+J

C-----------------UPDATE EIGENVALUES
            II=0
            DO 190 I=1,N
              II=II+I
              IF(VK(II).GT.ZERO.AND.VM(II).GT.ZERO) GO TO 190
              WRITE(MP,2000)I
              STOP
190           VALP(I)=VK(II)/VM(II)

            IF (M.GT.1) WRITE(MP,2010)IC,(VALP(I),I=1,N)
2010        FORMAT(/' EIGENVALUES, SWEEP ',I4/(1X,10E12.5))

C--------------------CHECK FOR CONVERGENCE OF EIGENVALUES
            DO 200 I=1,N
              IF (ABS(VALP(I)-VALP0(I)).GT.(EPS*VALP0(I))) GO TO 230  
200         CONTINUE

C--------------------CHECK FOR CONVERGENCE OF DIAGONAL TERMS
            JJ=1
            DO 210 J=2,N
              JJ=JJ+J
              JM1=J-1
              II=0
              DO 210 I=1,JM1
                II=II+I
                IJ=JJ-J+I
                FK=VK(IJ)*VK(IJ)/(VK(II)*VK(JJ))
                FM=VM(IJ)*VM(IJ)/(VM(II)*VM(JJ))
                IF (FK.GT.EPS2.OR.FM.GT.EPS2) GO TO 230
210           CONTINUE

C--------------------NORMALIZE EIGENVECTORS
            JJ=0
            DO 220 J=1,N
              JJ=JJ+J
              C1=SQRT(VM(JJ))
              DO 220 I=1,N
220             VECT(I,J)=VECT(I,J)/C1

C----------------------ACHIEVED CONVERGENCE
            IF (M.GT.0)WRITE(MP,2020)IC,ITR
2020        FORMAT(15X,'CONVERGENCE IN ',I4,' SWEEPS AND ',I5,
     1                 ' TRANSFORMATIONS IN JACOBI')
            RETURN

C-----------------TRANSFER VALP TO VALP0
230         DO 240 I=1,N
240           VALP0(I)=VALP(I)
250       CONTINUE

C---------FAIL TO CONVERGE
      WRITE(MP,2030)NCYM
2030  FORMAT('** ERROR, CONVERGENCE FAILURE IN JACOBI IN ',I4,' SWEEPS')
      STOP
 
              
      END
     
