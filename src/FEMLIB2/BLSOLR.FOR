      SUBROUTINE BLSOLR
C**********************************************************************C
C BLSOLR.FOR							       C
C  20/3-1992	Kompaniet Numerisk lyd				       C
C  Tonni F. Johansen						       C
C**********************************************************************C
C  Call and execute  Block 'SOLR'				       C
C  To assemble distributed loads (element function 7)		       C
C  Programmert etter Dhatt & Touzot, s.404			       C
C**********************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)

      REAL*4 TBL

      include 'femak.cmn'

      COMMON VA(1)

      DIMENSION TBL(8)
      DATA TBL/4HFG  ,4HKE  ,4HFE  ,4HDLE ,4HKGS ,4HKGD ,4HKGI ,
     1       4HRES /

C---------------------------------------------------
C      IF (M1.EQ.0) M1=MR
C      IF (M2.EQ.0) M2=ME
      WRITE(MP,2000)M
2000  FORMAT(//'ASSEMBLING OF DISTIBUTED LOADS (M=',I2,')'/
     1	     1X,40('=')/)

      IF (LFG.EQ.1) CALL ESPACE(NEQ,1,TBL(1),LFG)
      IF (LKE.EQ.1) CALL ESPACE(NKE,1,TBL(2),LKE)
      IF (LFE.EQ.1) CALL ESPACE(NDLE,1,TBL(3),LFE)
      IF (LDLE.EQ.1) CALL ESPACE(NDLE,1,TBL(4),LDLE)
      IF (LKGS.EQ.1) CALL ESPACE(NKG,1,TBL(5),LKGS)
      IF (LKGD.EQ.1) CALL ESPACE(NEQ,1,TBL(6),LKGD)
      IF (NSYM.EQ.1.AND.LKGI.EQ.1) CALL ESPACE(NKG,1,TBL(7),LKGI)
      IF (LRES.EQ.1) CALL ESPACE(NDLT,1,TBL(8),LRES)

      CALL EXSOLR(VA(LLD),VA(LDIMP),VA(LLOCE),VA(LCORE),VA(LPRNE),
     1  VA(LPREE),VA(LNE),VA(LKE),VA(LKE1),VA(LIE),VA(LIE1),VA(LFE),
     2  VA(LKGS),VA(LKGD),VA(LKGI),VA(LFG),VA(LCORG),VA(LDLNC),VA(LNEQ),
     3  VA(LRES),VA(LDLE),VA(LNDSRC),VA(LTYPRSG),VA(LPRSG))

      RETURN
      END
