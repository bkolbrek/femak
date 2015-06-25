      SUBROUTINE RDCIMP(MCH,KLOCE,VCORE,VPRNE,VPREE,KNE,IGRP)
C********************************************************************C
C RDCIMP.FOR                                                         C
C  11/8-1992                                                         C
C  Tonni F. Johansen  NTH                                            C
C********************************************************************C
C  write information about impedance surface on an element to a file C
C  MCH.                                                              C
C   Based on  Dhatt & Touzot, WRELEM s.399                           C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)


      include 'femak.cmn'    

      DIMENSION KLOCE(1),VCORE(1),VPRNE(1),VPREE(1),KNE(1)

C----------------------------------------------------------------------
      READ(MCH) IEL,ITPE,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,CFLAG,
     1      (KLOCE(I),I=1,IDLE),(VCORE(I),I=1,ICE),
     2      (VPRNE(I),I=1,IPRNE),(VPREE(I),I=1,IPREE),
     3      (KNE(I),I=1,INEL),
     4      IFL,IGRP


      RETURN
      END
