      SUBROUTINE RDSRCE(MCH,KLOCE,VCORE,VPRNE,VPREE,KNE,VNDSRC,KTYPRSG)
C********************************************************************C
C RDSRCE.FOR                                                         C
C  30/8-1992                                                         C
C  Tonni F. Johansen   NTH                                           C
C********************************************************************C
C  write information about sources to a file MCH.                    C
C   Based on  Dhatt & Touzot, WRELEM p.399                           C
C********************************************************************C

      include 'femak.cmn'    

      INTEGER MCH
      COMPLEX*16 VNDSRC(1)
      INTEGER KLOCE(1),KNE(1),KTYPRSG(1)
      REAL*8 VCORE(1),VPRNE(1),VPREE(1)

      INTEGER I
C----------------------------------------------------------------------
      READ(MCH) IEL,ITPE,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,CFLAG,
     1      (KLOCE(I),I=1,IDLE),(VCORE(I),I=1,ICE),
     2      (VPRNE(I),I=1,IPRNE),(VPREE(I),I=1,IPREE),
     3      (KNE(I),I=1,INEL),
     4      IFL,IGSRC

c      IF(KTYPRSG(IGSRC).EQ.20) THEN
        IF (NDIM.EQ.2) THEN
          READ(MCH)(VNDSRC(I),I=1,NFREQ*4)
        ELSEIF (NDIM.EQ.3) THEN
          READ(MCH)(VNDSRC(I),I=1,NFREQ*10)
        ENDIF
C........End srctype = 20
c      ENDIF

      RETURN
      END
