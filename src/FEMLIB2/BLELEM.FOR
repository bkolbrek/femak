      SUBROUTINE BLELEM
C**********************************************************************C
C BLELEM.FOR							       C
C  14/2-1992	Kompaniet Numerisk lyd				       C
C  Tonni F. Johansen						       C
C**********************************************************************C
C  Call and execute  Block 'ELEM'				       C
C  To read element data 					       C
C  Programmert etter Dhatt & Touzot, s.395			       C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 TBL

c      COMMON/COOR/NDIM,NNT,NDLN
c      COMMON/PRND/NPRN
c      COMMON/PREL/NGPE,NPRE
c      COMMON/ELEM/NELT,NELTR,NELTC,NNEL,NTPE,NGRE,ME,NIDENT,NPG
c      COMMON/ASSE/NSYM,NKG
c      COMMON/RESO/NEQ
c      COMMON/ES/M,MR,MP,M1,M2
c      COMMON/LOC/LCORG,LDLNC,LNEQ,LDIMP,LPRNG,LPREG,LLD,LLOCE,LCORE,LNE,
c     &		 LPRNE,LPREE,LDLE,LKE,LFE,LKGS,LKGD,LKGI,LFG,LRES,LDLG,
c     &		 LME,LDLE0,LDLG0,LFG0

c      INTEGER NGIMP,NGSRC,NFREQ,NSIMP,NSRC,PRPAR,VELPAR,INTPAR,NPEXT
c      LOGICAL UPDTFRQ,CMAT
c      COMMON/AKUS/NGIMP,NGSRC,NFREQ,NSIMP,NSRC,PRPAR,VELPAR,INTPAR,
c     +          NPEXT,UPDTFRQ,CMAT

      include 'femak.cmn'

      COMMON VA(1)

      DIMENSION TBL(6)
      INTEGER IN(8)
      DATA TBL/4HLD  ,4HLOCE,4HCORE,4HNE  ,4HPRNE,4HPREE/

C----------------------------------------------------------------------
c      IF (M1.EQ.0) M1=MR
c      IF (M2.EQ.0) M2=ME
c      READ(M1,1000)IN
      READ(MR,1000)IN
1000  FORMAT(8I5)

      IF (IN(1).GT.0) NELT=IN(1)
      IF (IN(2).GT.0) NELTR=IN(2)
      IF (IN(3).GT.0) NELTC=IN(3)
      IF (IN(4).GT.0) NNEL=IN(4)
      IF (IN(5).GT.0) NTPE=IN(5)
      IF (IN(6).GT.0) NGRE=IN(6)
      IF (IN(7).NE.0) NSYM=1
      IF (IN(8).NE.0) NIDENT=1

      IF (NELTC.GT.0) CMAT=.TRUE.

      WRITE(MP,2000)M,NELT,NELTR,NELTC,NNEL,NTPE,NGRE,NSYM,NIDENT
2000  FORMAT(//' INPUT OF ELEMENTS (M=',I2,')'/' ',20('=')/
     1	     15X,'MAX. NUMBER OF ELEMENTS              (NELT)=',I5/
     2	     15X,'MAX. NO. ELEMENTS, REAL COEFF.      (NELTR)=',I5/
     3	     15X,'MAX. NO. ELEMENTS, COMPLEX COEFF.   (NELTC)=',I5/
     4	     15X,'MAX. NUMBER OF NODES PER ELEMENTS    (NNEL)=',I5/
     5	     15X,'DEFAULT ELEMENT TYPE                 (NTPE)=',I5/
     6	     15X,'NUMBER OF GROUPS OF ELEMENTS         (NGRE)=',I5/
     7	     15X,'INDEX FOR NON SYMMETRIC PROBLEM      (NSYM)=',I5/
     8	     15X,'INDEX FOR IDENTICAL ELEMENTS       (NIDENT)=',I5/)

      IF (LLD.EQ.1) CALL ESPACE(NEQ+1,0,TBL(1),LLD)
      IF (LLOCE.EQ.1) CALL ESPACE(NNEL*NDLN,0,TBL(2),LLOCE)
      IF (LCORE.EQ.1) CALL ESPACE(NNEL*NDIM,1,TBL(3),LCORE)
      IF (LNE.EQ.1) CALL ESPACE(NNEL,0,TBL(4),LNE)
      IF (NPRN.GT.0.AND.LPRNE.EQ.1)
     1		  CALL ESPACE(NNEL*NPRN,1,TBL(5),LPRNE)
      IF (NPRE.GT.0.AND.LPREE.EQ.1) CALL ESPACE(NPRE,1,TBL(6),LPREE)

      CALL EXELEM(VA(LCORG),VA(LDLNC),VA(LPRNG),VA(LPREG),VA(LLOCE),
     1	     VA(LCORE),VA(LNE),VA(LPRNE),VA(LPREE),VA(LNEQ),VA(LLD))


      WRITE(MP,2010) NKG,NPG
2010  FORMAT(15X,'LENGTH OF A TRIANGLE IN KG           (NKG) =',I10/
     1	     15X,'NUMBER OF INTEGRATION POINTS         (NPG) =',I10/)


      RETURN
      END
