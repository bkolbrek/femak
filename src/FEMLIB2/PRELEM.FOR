      SUBROUTINE PRELEM(KLOCE,VCORE,VPRNE,VPREE,KNE)
C********************************************************************C
C PRELEM.FOR                                                         C
C  14/2-1992                                                         C
C  Tonni F. Johansen  NTH                                            C
C********************************************************************C
C  print data defining an element                                    C
C  Programmed after Dhatt & Touzot, p.399                            C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)


      include 'femak.cmn'    

      DIMENSION KLOCE(1),VCORE(1),VPRNE(1),VPREE(1),KNE(1)

C----------------------------------------------------------------------
      IF (M.GE.0) WRITE(MP,2000) IEL,ITPE,INEL,IDLE,IPRNE,IPREE,IGRE
2000  FORMAT(10X,'ELEMENT:',I5,' TYPE:',I2,' N.P.:',I2,'D.O.F.:',
     1	      I3,' N. PROP.:',I3,' EL. PROP:',I3,' GROUP:',I3)
      IF (M.GE.0) WRITE(MP,2010)(KNE(I),I=1,INEL)
2010  FORMAT(15X,'CONNECTIVITY (NE)',20I5/(32X,20I5))
      IF (M.LT.1) GO TO 10
      WRITE(MP,2020)(KLOCE(I),I=1,IDLE)
2020  FORMAT(15X,'LOCALIZATION (LOCE)',20I5/(32X,20I5))
      WRITE(MP,2030)(VCORE(I),I=1,ICE)
2030  FORMAT(15X,'COORDINATES (CORE)',8E12.5/(32X,8E12.5))
      IF (NPRN.GT.0) WRITE(MP,2040)(VPRNE(I),I=1,IPRNE)
2040  FORMAT(15X,'NOD. PROP. (PRNE)',8E12.5/(32X,8E12.5))
      IF (IPREE.GT.0) WRITE(MP,2050)(VPREE(I),I=1,IPREE)
2050  FORMAT(15X,'ELEM. PROP. (PREE)',8E12.5/(32X,8E12.5))

10    RETURN
      END
