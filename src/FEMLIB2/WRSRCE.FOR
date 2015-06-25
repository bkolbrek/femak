      SUBROUTINE WRSRCE(MCH,KLOCE,VCORE,VPRNE,VPREE,KNE,VNDSRC,KTYPRSG)
C********************************************************************C
C WRSRCE.FOR                                                         C
C  2/1-1995                                                          C
C  Tonni F. Johansen  NTH                                            C
C********************************************************************C
C  write information about sources to a file MCH                     C
C   Based on  Dhatt & Touzot, WRELEM p.399                           C
C********************************************************************C

      include 'femak.cmn'    

      INTEGER MCH
      COMPLEX*16 VNDSRC(1)
      INTEGER KLOCE(1),KNE(1),KTYPRSG(1)
      REAL*8 VCORE(1),VPRNE(1),VPREE(1)

      INTEGER IPRNE1,IPREE1,I
C----------------------------------------------------------------------
      IPRNE1=IPRNE
      IF (IPRNE1.EQ.0) IPRNE1=1
      IPREE1=IPREE
      IF (IPREE1.EQ.0) IPREE1=1

      WRITE(MCH) IEL,ITPE,IGRE,IDLE,ICE,IPRNE1,IPREE1,INEL,CFLAG,
     1      (KLOCE(I),I=1,IDLE),(VCORE(I),I=1,ICE),
     2      (VPRNE(I),I=1,IPRNE1),(VPREE(I),I=1,IPREE1),
     3      (KNE(I),I=1,INEL),
     4       IFL,IGSRC

C      IF(KTYPRSG(IGSRC).EQ.20) THEN
        IF (NDIM.EQ.2) THEN
          WRITE(MCH)(VNDSRC(I),I=1,NFREQ*4)
        ELSEIF (NDIM.EQ.3) THEN
          WRITE(MCH)(VNDSRC(I),I=1,NFREQ*10)
        ENDIF
C      ENDIF

      RETURN
      END
