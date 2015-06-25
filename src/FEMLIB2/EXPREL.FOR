      SUBROUTINE EXPREL(VPREG,V1)
C**********************************************************************C
C EXPREL.FOR							       C
C  4/2-1992	Kompaniet Numerisk lyd				       C
C  Tonni F. Johansen						       C
C**********************************************************************C
C  Call and execute  Block 'PREL'				       C
C  Read element properties					       C
C  Programmert etter Dhatt & Touzot, s.393			       C
C**********************************************************************C

      IMPLICIT REAL*8 (A-H,O-Z)

c      COMMON/PREL/NGPE,NPRE
c      COMMON/ES/M,MR,MP,M1

      include 'femak.cmn'

      DIMENSION VPREG(1),V1(1)

C----------------------------------------------------------------------
      IF (M.GE.0) WRITE(MP,2000)
2000  FORMAT(//' CARDS OF ELEMENT PROPERTIES '/)

C--------READ A GROUP
      I1=MIN0(7,NPRE)
      J=1

C10    READ(M1,1000) IGPE,(V1(I),I=1,I1)
10    READ(MR,1000) IGPE,(V1(I),I=1,I1)
C000  FORMAT(I5,7F10.0)
1000  FORMAT(I5,7e12.5)
      IF (M.GE.0) WRITE(MP,2010) IGPE,(V1(I),I=1,I1)
2010  FORMAT(' >>>>> ',I5,7E12.5)


      IF (IGPE.LE.0) GO TO 40
      IF (IGPE.GT.NGPE) CALL ERREUR(41,IGPE,NGPE,1)
      IF (NPRE.LE.7) GO TO 20

C-------------READ THE PROPERTIES
C      READ(M1,1010)(V1(I),I=8,NPRE)
      READ(MR,1010)(V1(I),I=8,NPRE)
c1010  FORMAT(5X,7F10.0)
1010  FORMAT(7e12.5)
      IF (M.GE.0) WRITE(MP,2020)(V1(I),I=8,NPRE)
2020  FORMAT(' >>>>>  ',5X,7E12.5)

20    DO 30 I=1,NPRE
	VPREG(J)=V1(I)
	J=J+1
30    CONTINUE

      GO TO 10


40    RETURN
      END
