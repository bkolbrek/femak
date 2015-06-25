      SUBROUTINE BLSOLC
C**********************************************************************C
C BLSOLC.FOR                                                           C
C  2/3-1992   Kompaniet Numerisk lyd                                   C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  Call and execute  Block 'SOLC'                                      C
C  To read concentrated loads                                          C 
C     Programmert etter Dhatt & Touzot, s.400                          C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*4 TBL
c      COMMON /RESO/NEQ
c      COMMON /ES/M,MR,MP,M1
c      COMMON /LOC/LCORG,LDLNC,LNEQ,LXX(15),LFG

      include 'femak.cmn'

      COMMON VA(1)

      DATA TBL/4HFG  /

C------------------------------------------------------------------------
C      IF (M1.EQ.0)M1=MR
      WRITE(MP,2000)M
2000  FORMAT(//' INPUT OF CONCENTRATED LOADS (M=',I2,')'/' ',
     1      39('='))
      IF (LFG.EQ.1) CALL ESPACE(NEQ,1,TBL,LFG)
      CALL EXSOLC(VA(LFG),VA(LDLNC),VA(LNEQ))

      RETURN
      END
