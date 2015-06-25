      SUBROUTINE BLCOMT
C**********************************************************************C
C BLCOMT.FOR							       C
C  18/1-1992   Kompaniet Numerisk lyd				       C
C  Tonni F. Johansen						       C
C**********************************************************************C
C  Call and execute  Block 'COMT'				       C
C  Programmert etter Dhatt & Touzot, s.384			       C
C**********************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 BLANC,CART
c      COMMON/ES/M,MR,MP
      COMMON/TRVL/CART(20)

      include 'femak.cmn'

      DATA BLANC/4H    /

C-----------------------------------------------------------------------
      WRITE(MP,2000)
2000  FORMAT(///'COMMENTS'/' ',10('=')/)
c-----------Read a comment card
10    READ(MR,1000)CART
1000  FORMAT(20A4)
C-----------Search for a wholly blanc card
      DO 20 I=1,20
	IF(CART(I).NE.BLANC) GO TO 30
20    CONTINUE
      RETURN
30    WRITE(MP,2010) CART
2010  FORMAT(1X,20A4)
      GO TO 10

      END
