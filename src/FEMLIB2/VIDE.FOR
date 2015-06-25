      SUBROUTINE VIDE(IDEB,IREEL,TBL)
C****************************************************************************C
C   vide.for								     C
C   1/2-1992								     C
C   Kompaniet Numerisk lyd, Tonni F. Johansen                                C
C****************************************************************************C
C  To delete a table from VA, followd by compacting			     C
C    Input                                                                   C
C	IDEB		First position of table to be deleted		     C
C       IREEL           Table type                                           C
C                          .EQ.0    Integer                                  C
C                          .EQ.1    Real                                     C
C       TBL             Name of the table (A4)                               C
C                                                                            C
C	 Programmed after Dhatt & Touzot p.368-369			     C
C****************************************************************************C

      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 TBL
      
c      COMMON/ES/M,MR,MP
c      COMMON/ALLOC/NVA,IVA,IVAMAX,NREEL,NTBL
c      COMMON/LOC/LXX(25)

      INTEGER LXX(53)
      EQUIVALENCE(LXX(1),LCORG)
      
      include 'femak.cmn'    

      COMMON VA(1)

C---------------------------------------------------------------------------
C----------------SEARCH FOR THE FIRST POSITION OF NEXT TABLE
      I1=IVA+1
      DO 10 I=1,NTBL
	IF (LXX(I).LE.IDEB) GO TO 10
	IF (LXX(I).LT.I1) I1=LXX(I)
10    CONTINUE

C----------SHIFT ALL TABLES AFTER THIS
      ID=I1-IDEB
      IF (I1.EQ.IVA-1) GO TO 40

      DO 20 I=1,NTBL
	IF (LXX(I).GT.IDEB) LXX(I)=LXX(I)-ID
20    CONTINUE

      DO 30 I=I1,IVA
	J=I-ID
	VA(J)=VA(I)
30    CONTINUE

C----------PRINT
40    IVA=IVA-ID
      IF (M.GT.0) WRITE(MP,2000)TBL,ID,IDEB
2000  FORMAT(60X,'DELETED TABLE ',A4,' COMPACTING ',I7,
     1	' REAL WORDS AFTER VA(',I7,')')

      RETURN
      END
