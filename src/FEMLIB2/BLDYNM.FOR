      SUBROUTINE BLDYNM
C***********************************************************************C
C     BLDYNM.FOR                                                        C
C***********************************************************************C
C     3/1-1995                                                          C
C     NTH                                                               C
C     Tonni F. Johansen                                                 C
C***********************************************************************C
C     TO CALL BLOCK "DYNM"                                              C
C     ASSEMBLE AND SOLVE A DYNAMIC PROBLEM IN CORE                      C
C     Primarely made for acoustic problems, however extension to plates C
C     and related type systems is wanted.                               C
C***********************************************************************C
C  Some comments:                                                       C
C   We builds this block on the structure of LINM. It is changed in     C
C   dynamical problems. This means complex number, loops for multiple   C
C   frequencies, and simple updating for each frequency.                C
C***********************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)
      include 'femak.cmn'

      COMMON VA(1)
C $DECLARE
      REAL*4 TBL(21)
      DATA TBL/4HKGS ,4HKGD ,4HKGI ,4HMGS ,4HMGD ,4HMGI ,4HIGS ,4HIGD ,
     1 4HIGI ,4HFG  ,4HKE  ,4HFE  ,4HRES ,4HDLE ,4HSRC ,4HIE  ,4HSRE ,
     2 4HFREQ,4HCOXT,4HKE1 ,4HIE1 /

C-----------------------------------------------------------------------
      READ(MR,1010)NFREQ,PRPAR,VELPAR,INTPAR,NPEXT
1010  FORMAT(5I5)

      WRITE(MP,2000) M,NRES,NFREQ,PRPAR,VELPAR,INTPAR,NPEXT

2000  FORMAT(//' ASSEMBLING AND DYNAMIC SOLUTION (M=',I2,')'/' ',
     1 30('=')/ 15X,' INDEX FOR RESIDUAL COMPUTATION    (NRES)=',I5,/
     2 15X,' NUMBER OF FREQUENCIES                     (NFREQ)=',I5,/
     3 15X,' PRINT PARAMETER                           (PRPAR)=',I5,/
     4 15X,' VELOCITY PARAMETER                       (VELPAR)=',I5,/
     5 15X,' INTENCITY PARAMETER                      (INTPAR)=',I5,/
     6 15X,' NO. OF EXTERNAL CALC. POINTS              (NPEXT)=',I5)

C...........Allocate space for calculation of external field.
      IF(NPEXT.GT.0) THEN
	IF(LCOEXT.EQ.1) CALL ESPACE(NDIM*NPEXT,1,TBL(19),LCOEXT)
	CALL RDEXTFLD(VA(LCOEXT))
      ENDIF


C.........Allocate space for system matrices
      IF(LKGS.EQ.1) CALL ESPACE(NKG,1,TBL(1),LKGS)
      IF(LKGD.EQ.1) CALL ESPACE(NEQ,1,TBL(2),LKGD)
      IF(NSYM.EQ.1.AND.LKGI.EQ.1) CALL ESPACE(NKG,1,TBL(3),LKGI)
      IF(LMGS.EQ.1) CALL ESPACE(NKG,1,TBL(4),LMGS)
      IF(LMGD.EQ.1) CALL ESPACE(NEQ,1,TBL(5),LMGD)
      IF(NSYM.EQ.1.AND.LMGI.EQ.1) CALL ESPACE(NKG,1,TBL(6),LMGI)
      IF(LIGS.EQ.1) CALL ESPACE(NKG,2,TBL(7),LIGS)
      IF(LIGD.EQ.1) CALL ESPACE(NEQ,2,TBL(8),LIGD)
      IF(NSYM.EQ.1.AND.LIGI.EQ.1) CALL ESPACE(NKG,2,TBL(9),LIGI)


      IF(LFG.EQ.1) CALL ESPACE(NEQ,1,TBL(10),LFG)
      IF(LSRC.EQ.1) CALL ESPACE(NEQ,2,TBL(15),LSRC)

      IF(LKE.EQ.1) CALL ESPACE(NKE,1,TBL(11),LKE)
      IF((LKE1.EQ.1).AND.(NSYM.EQ.1)) THEN
        CALL ESPACE(NDLE*(NDLE+1)/2,2,TBL(21),LKE1)
      ENDIF
      IF(LFE.EQ.1) CALL ESPACE(NDLE,1,TBL(12),LFE)
      IF(LRES.EQ.1) CALL ESPACE(NDLT,1,TBL(13),LRES)
      IF(LDLE.EQ.1) CALL ESPACE(NDLE,1,TBL(14),LDLE)
      IF(LIE.EQ.1) CALL ESPACE(NKE,2,TBL(16),LIE)

      IF(LSRE.EQ.1) CALL ESPACE(NDLE,2,TBL(17),LSRE)
      IF((LIE1.EQ.1).AND.(NSYM.EQ.1)) THEN
        CALL ESPACE(NDLE*(NDLE+1)/2,2,TBL(21),LIE1)
      ENDIF
      IF(LFREQ.EQ.1) CALL ESPACE(NFREQ,1,TBL(18),LFREQ)


      CALL EXDYNM(VA(LLD),VA(LDIMP),VA(LLOCE),VA(LCORE),VA(LPRNE),
     1        VA(LPREE),VA(LNE),VA(LKE),VA(LKE1),VA(LFE),VA(LKGS),
     2        VA(LKGD),VA(LKGI),VA(LMGS),VA(LMGD),VA(LMGI),VA(LIGS),
     3        VA(LIGD),VA(LIGI),VA(LFG),VA(LSRC),VA(LCORG),
     4        VA(LDLNC),VA(LNEQ),VA(LRES),VA(LDLE),VA(LIE),VA(LIE1),     
     5        VA(LSRE),VA(LFREQ),VA(LTYPRIG),VA(LPRIG),VA(LNDSRC),
     6        VA(LTYPRSG),VA(LPRSG),VA(LPRNG))


      RETURN
      END
