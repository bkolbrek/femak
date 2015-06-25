      SUBROUTINE BLLINM
C***********************************************************************
C     BLELEM.FOR
C***********************************************************************
C     21/2-1992
C     Kompaniet Numerisk Lyd
C     ULF R. KRISTIANSEN
C     Programmert etter Dhatt og Touzot s.408
C***********************************************************************
C     TO CALL BLOCK "LINM"
C     ASSEMBLE AND SOLVE A LINEAR PROBLEM IN CORE
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include 'femak.cmn'

      COMMON VA(1)
C   $DECLARE
      INTEGER IN

      REAL*4 TBL(8)
      DATA TBL/4HKGS ,4HKGD ,4HKGI ,4HFG  ,4HKE  ,4HFE  ,4HRES ,4HDLE /
C-----------------------------------------------------------------------
      READ(MR,1000) IN
1000  FORMAT(1I5)
      IF(IN.NE.0) NRES=1
      WRITE(MP,2000) M,NRES
2000  FORMAT(//' ASSEMBLING AND LINEAR SOLUTION (M=',I2,')'/' ',30('=')/
     1  15X,' INDEX FOR RESIDUAL COMPUTATION    (NRES)=',I5)
      IF(LKGS.EQ.1) CALL ESPACE(NKG,1,TBL(1),LKGS)
      IF(LKGD.EQ.1) CALL ESPACE(NEQ,1,TBL(2),LKGD)
      IF(NSYM.EQ.1.AND.LKGI.EQ.1) CALL ESPACE(NKG,1,TBL(3),LKGI)
      IF(LFG.EQ.1) CALL ESPACE(NEQ,1,TBL(4),LFG) 
      IF(LKE.EQ.1) CALL ESPACE(NKE,1,TBL(5),LKE)
      IF(LFE.EQ.1) CALL ESPACE(NDLE,1,TBL(6),LFE)
      IF(LRES.EQ.1) CALL ESPACE(NDLT,1,TBL(7),LRES)
      IF(LDLE.EQ.1) CALL ESPACE(NDLE,1,TBL(8),LDLE)
      CALL EXLINM(VA(LLD),VA(LDIMP),VA(LLOCE),VA(LCORE),VA(LPRNE),
     1      VA(LPREE),VA(LNE),VA(LKE),VA(LKE1),VA(LIE),VA(LIE1),VA(LFE),
     2      VA(LSRE),VA(LKGS),VA(LKGD),VA(LKGI),VA(LFG),VA(LCORG),
     3      VA(LDLNC),VA(LNEQ),VA(LRES),VA(LDLE),VA(LNDSRC),VA(LTYPRSG),
     4      VA(LPRSG))
      RETURN
      END
