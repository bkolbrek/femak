      SUBROUTINE RDELEM(MCH,KLOCE,VCORE,VPRNE,VPREE,KNE)
C********************************************************************C
C RDELEM.FOR                                                         C
C  14/2-1992                                                         C
C  Tonni F. Johansen   NTH                                           C
C********************************************************************C
C  read element properties from file MCH                             C
C  Programmed after Dhatt & Touzot, p.399                            C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)


      include 'femak.cmn'    

      DIMENSION KLOCE(1),VCORE(1),VPRNE(1),VPREE(1),KNE(1)

C----------------------------------------------------------------------
      READ(MCH,END=999) IEL,ITPE,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,CFLAG,
     1      (KLOCE(I),I=1,IDLE),(VCORE(I),I=1,ICE),
     2      (VPRNE(I),I=1,IPRNE),(VPREE(I),I=1,IPREE),
     3      (KNE(I),I=1,INEL)


999   RETURN
      END
