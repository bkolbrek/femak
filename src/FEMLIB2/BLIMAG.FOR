      SUBROUTINE BLIMAG
C**********************************************************************C
C BLIMAG.FOR							       C
C  18/1-1992   Kompaniet Numerisk lyd				       C
C  Tonni F. Johansen						       C
C**********************************************************************C
C  Call and execute  Block 'IMAG'				       C
C  To print out the image of data cards 			       C
C  Programmert etter Dhatt & Touzot, s.383			       C
C**********************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)
c      COMMON/ES/M,MR,MP,M1
c      COMMON/TRVL/CART(20)

      include 'femak.cmn'

      DATA ICARTM/40/

C-----------------------------------------------------------------------
C      IF (M1.EQ.0)  M1=MR
      WRITE(MP,2000)
2000  FORMAT(///,1X,'IMAGE OF DATA CARDS'/1X,28('='),/)
      WRITE(MP,2005)
2005  FORMAT(/
     &	50X,'C O L U M N   N U M B E R',/,13X,'CARD',6X,
     &	10X,'1',9X,'2',9X,'3',9X,'4',9X,'5',9X,'6',9X,'7',9X,'8',/,
     &	12X,'NUMBER',6X,8('1234567890'),/,12X,8('-'),6X,80('-'))
      ICART=0
      ICART1=0
C10    READ(M1,1000,END=30) CART
10    READ(MR,1000,END=30) CART
1000  FORMAT(20A4)
      ICART=ICART+1
      ICART1=ICART1+1
      IF(ICART1.LE.ICARTM) GO TO 20
      WRITE(MP,2010)
2010  FORMAT(12X,8(1H-),6X,80(1H-),/,13X,'CARD',7X,8('1234567890'),/,
     &	12X,'NUMBER',6X,9X,'1',9X,'2',9X,'3',9X,'4',9X,'5',9X,'6',
     &	9X,'7',9X,'8',/,50X,'C O L U M N   N U M B E R')
      WRITE(MP,2015)
2015  FORMAT(1H1,//)
      WRITE(MP,2005)
      ICART1=0
20    WRITE(MP,2020)ICART,CART
2020  FORMAT(10X,I10,6X,20A4)
      GO TO 10
30    WRITE(MP,2010)
      WRITE(MP,2030)
2030  FORMAT(///,51X,'E N D   O F   D A T A',/,1H1)
C      REWIND M1
C      READ(M1,1000) CART
      REWIND MR
      READ(MR,1000) CART
      RETURN
      END
