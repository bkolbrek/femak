      SUBROUTINE BLTEMP
C**********************************************************************C
C  BLTEMP.FOR                                                          C
C  02/12-1994   Kompaniet Numerisk lyd                                 C
C  Bruno Brouard                                                       C
C**********************************************************************C
C  Call and execute  Block 'TEMP'                                      C
C  Programmed after Dhatt & Touzot, p.431                              C
C**********************************************************************C
C     TO CALL BLOCK 'TEMP'
C     TO SOLVE AN UNSTEADY PROBLEM (LINEAR OR NOT)
C**********************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)
      include 'femak.cmn'

      REAL*4 TBL
      COMMON VA(1)
      DIMENSION TBL(13),IN(2),XIN(3)
      DATA TBL/4HKGS ,4HKGD ,4HKGI ,4HFG  ,4HKE  ,
     1         4HFE  ,4HRES ,4HDLE ,4HDLG ,4HME  ,4HDLE0,4HDLG0,4HFG0 /
C----------------------------------------------------------------------- 
C      IF (M1.EQ.0) M1=MR
C      IF (M2.EQ.0) M2=ME
      WRITE(MP,2000) M
 2000 FORMAT(//' UNSTEADY SOLUTION (M=',I2,')'/1X,23('='))
C------- TO ALLOCATE SPACE
      IF(LKGS.EQ.1) CALL ESPACE(NKG,1,TBL(1),LKGS)
      IF(LKGD.EQ.1) CALL ESPACE(NEQ,1,TBL(2),LKGD)
      IF(NSYM.EQ.1.AND.LKGI.EQ.1) CALL ESPACE(NKG,1,TBL(3),LKGI)
      IF(LFG.EQ.1) CALL ESPACE(NEQ,1,TBL(4),LFG)
      IF(LKE.EQ.1) CALL ESPACE(NKE,1,TBL(5),LKE)
      IF(LFE.EQ.1) CALL ESPACE(NDLE,1,TBL(6),LFE)
      IF(LRES.EQ.1) CALL ESPACE(NEQ,1,TBL(7),LRES)
      IF(LDLE.EQ.1) CALL ESPACE(NDLE,1,TBL(8),LDLE)
      IF(LDLG.EQ.1) CALL ESPACE(NEQ,1,TBL(9),LDLG)
      IF(LME.EQ.1) CALL ESPACE(NKE,1,TBL(10),LME)
      IF(LDLE0.EQ.1) CALL ESPACE(NDLE,1,TBL(11),LDLE0)
      IF(LDLG0.EQ.1) CALL ESPACE(NEQ,1,TBL(12),LDLG0)
      IF(LFG0.EQ.1) CALL ESPACE(NEQ,1,TBL(13),LFG0)
C------- TO EXECUTE THE BLOCK
      CALL EXTEMP(VA(LCORG),VA(LDLNC),VA(LDIMP),VA(LNEQ),VA(LLD),
     1  VA(LLOCE),VA(LCORE),VA(LPRNE),VA(LPREE),VA(LNE),VA(LKE),VA(LME),
     2  VA(LFE),VA(LDLE),VA(LKGS),VA(LKGD),VA(LKGI),VA(LFG),VA(LRES),
     3  VA(LDLG),VA(LDLE0),VA(LDLG0),VA(LFG0))
      RETURN
      END
      