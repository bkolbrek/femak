      SUBROUTINE MODFG(IDLE,NSYM,KLOCE,VDIMP,VKE,VFG)
C***********************************************************
C     MODFG.FOR
C***********************************************************
C     24/2 1992
C     Kompaniet Numerisk Lyd
C     Ulf R. Kristiansen
C     programmert etter Dhatt og Touzot s. 414
C***********************************************************
C     TO MODIFY VECTOR FG TO TAKE INTO ACCOUNT OF PRESCRIBED NON ZERO
C     D.O.F.  FOR A GIVEN ELEMENT
C       INPUT
C          IDLE     ELEMENT NUMBER OF D.O.F.
C          NSYM     0=SYMMETRIC PROBLEM, 1=NON SYMMETRIC PROBLEM
C          KLOCE    ELEMENT LOCALIZATION VECTOR
C          VDIMP    VALUES OF PRESCRIBED D.O.F.
C          VKE      ELEMENT MATRIX (FULL OR UPPER TRIANGLE 
C                   BY DESCENDING COLUMNS)
C       OUTPUT
C          VFG      GLOBAL LOAD VECTOR
C====================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION KLOCE(1),VDIMP(1),VKE(1),VFG(1)
      DATA ZERO/0.D0/
C--------------------------------------------------------------------
      IEQ0=IDLE
      IEQ1=1
C------  FOR EACH ROW OF ELEMENT MATRIX
      DO 50 JD=1,IDLE
      IF(NSYM.NE.1) IEQ0=JD
      IEQ=IEQ1
      JL=KLOCE(JD)
      IQ=1
      IF(JL) 10,50,50
10    JL=-JL
      DIMP=VDIMP(JL)
      IF(DIMP.EQ.ZERO) GO TO 50
C------  FOR EACH COLUMN OF ELEMENT MATRIX
      DO 40 ID=1,IDLE
      IL=KLOCE(ID)
      IF(NSYM.EQ.1) GO TO 30
      IF(ID-JD) 30,20,20
20    IQ=ID
30    IF(IL.GT.0) VFG(IL)=VFG(IL)-VKE(IEQ)*DIMP
40    IEQ=IEQ+IQ
50    IEQ1=IEQ1+IEQ0
      RETURN
      END
