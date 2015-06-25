      SUBROUTINE WRELEM(MCH,KLOCE,VCORE,VPRNE,VPREE,KNE)
C********************************************************************C
C WRELEM.FOR                                                         C
C  14/2-1992                                                         C
C  Tonni F. Johansen   NTH                                           C
C********************************************************************C
C  write element properties of file MCH                              C
C  Programmed after Dhatt & Touzot, p.399                            C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)


      include 'femak.cmn'    

      DIMENSION KLOCE(1),VCORE(1),VPRNE(1),VPREE(1),KNE(1)

C----------------------------------------------------------------------
      IPRNE1=IPRNE
      IF (IPRNE1.EQ.0) IPRNE1=1
      IPREE1=IPREE
      IF (IPREE1.EQ.0) IPREE1=1

      WRITE(MCH) IEL,ITPE,IGRE,IDLE,ICE,IPRNE1,IPREE1,INEL,CFLAG,
     1      (KLOCE(I),I=1,IDLE),(VCORE(I),I=1,ICE),
     2      (VPRNE(I),I=1,IPRNE1),(VPREE(I),I=1,IPREE1),
     3      (KNE(I),I=1,INEL)


      RETURN
      END
