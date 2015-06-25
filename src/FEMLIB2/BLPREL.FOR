      SUBROUTINE BLPREL
C**********************************************************************C
C BLPREL.FOR							       C
C  18/1-1992   Kompaniet Numerisk lyd				       C
C  Tonni F. Johansen						       C
C**********************************************************************C
C  Call and execute  Block 'PREL'				       C
C  To read element properties					       C
C  Programmert etter Dhatt & Touzot, s.393			       C
C**********************************************************************C

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 TBL

c      COMMON/PREL/NGPE,NPRE
c      COMMON/ES/M,MR,MP,M1
c      COMMON/LOC/LXX(5),LPREG
      COMMON/TRVL/IN(2)

      include 'femak.cmn'

      COMMON VA(1)

      DIMENSION TBL(2)
      DATA TBL/4HPREG,4HV   /

C----------------------------------------------------------------------
C     IF (M1.EQ.0) M1=MR
C--------READ NUMBER OF GROUPS AND PROPERTIES PER GROUP
CC     READ(M1,1000) IN
      READ(MR,1000) IN
1000  FORMAT(2I5)

      IF (IN(1).GT.0) NGPE=IN(1)
      IF (IN(2).GT.0) NPRE=IN(2)
      WRITE(MP,2000) M,NGPE,NPRE
2000  FORMAT(//' INPUT OF ELEMENT PROPERTIES (M=',I2,')'/' ',
     1	 35('=')/15X,'NUMBER OF GROUPS OF PROPERTIES   (NGPE)=',I5/
     2	 15X,'MUMBER OF PROPERTIES PER GROUP    (NPRE)=',I5)

      IF (LPREG.EQ.1) CALL ESPACE(NGPE*NPRE,1,TBL(1),LPREG)
      CALL ESPACE(NPRE,1,TBL(2),L1)
      CALL EXPREL(VA(LPREG),VA(L1))
      CALL VIDE(L1,1,TBL(2))

      RETURN
      END
