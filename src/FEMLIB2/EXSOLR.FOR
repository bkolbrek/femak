      SUBROUTINE EXSOLR(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1    VIE,VIE1,VFE,VKGS,VKGD,VKGI,VFG,VCORG,KDLNC,KNEQ,VRES,VDLE,
     2    VNDSRC,KTYPRSG,VPRSG)
C**********************************************************************C
C EXSOLR.FOR                                                           C
C  3/1-1995          NTH                                               C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  To execute  Block 'SOLR'                                            C
C  To assemble distributed loads (element function 7)                  C
C  Programmed after Dhatt & Touzot, p.405                              C
C**********************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)

      include 'femak.cmn'
C$DECLARE

      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VKE1(1),VFE(1)
      COMPLEX*16 VIE(1),VIE1(1),VNDSRC(1),VPRSG(1)
      INTEGER KTYPRSG(1)

      INTEGER KLD(1),KLOCE(1),KNE(1)
      REAL*8 VDIMP(1),VKGS(1),VKGD(1),VKGI(1),VFG(1),VRES(1)

      INTEGER I
      
C---------------------------------------------------
C------------ASSEMBLE FG
C      CALL ASFG(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VFE,VKGS,
C     1		VKGD,VKGI,VFG,VDLE,VRES)
      CALL ASFG(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1  VIE,VIE1,VFE,VKGS,VKGD,VKGI,VFG,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG)

C-----------------OUTPUT
      IF (M.GE.1) WRITE(MP,2000) (VFG(I),I=1,NEQ)
2000  FORMAT(/'GLOBAL LOAD VECTOR   (FG)'/(1X,10E12.5))



      RETURN
      END
