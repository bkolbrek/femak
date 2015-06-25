      SUBROUTINE ASGRAD(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1  VIE,VIE1,VFE,VSRE,VKGS,VKGD,VKGI,VFG,VDLE,VRES,VNDSRC,
     2  KTYPRSG,VPRSG)
C********************************************************************** 
C     ASGRAD.FOR
C**********************************************************************
C     3/1 1995
C     NTH
C     Ulf R. Kristiansen
C     Programmed after Dhatt og Touzot s.410
C**********************************************************************  
C     TO EVALUATE AND PRINT GRADIENTS ( STRESSES ) AT ELEMENT G.P.
C       ( ELEMENT FUNCTION 8)
C======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      include 'femak.cmn'
C $DECLARE

      COMPLEX*16 VIE(1),VIE1(1),VSRE(1)
      REAL*8 VDIMP(1),VCORE(1),VPRNE(1),VPREE(1),VKE(1),VFE(1),VDLE(1)
      REAL*8 VRES(1),VKE1(1),VKGS(1),VKGD(1),VKGI(1)
      INTEGER KLD(1),KLOCE(1),KNE(1)

      COMPLEX*16 VNDSRC,VPRSG(1)
      INTEGER KTYPRSG(1)
      
      INTEGER IE
      
C----------------------------------------------------------------------
C---------------  REWIND ELEMENTS FILE (M1)
      REWIND M1
      ITPE1=0
C--------------   LOOP OVER THE ELEMENTS
      DO IE=1,NELT
C-------------    READ THE ELEMENT
          CALL RDELEM(M1,KLOCE,VCORE,VPRNE,VPREE,KNE)
C-------------    EVALUATE INTERPOLATION FUNCTION IF REQUIRED
        IF(ITPE.EQ.ITPE1) GO TO 10
        ICODE=2
        CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VSRE,VKE,VKE1,VIE,VIE1,VFE,
     1            VSRE,VNDSRC,KTYPRSG,VPRSG)
C------------     FIND ELEMENT D.O.F.
10      CALL DLELM(KLOCE,VFG,VDIMP,VDLE)
C------------     COMPUTE AND PRINT STRESSES OR GRADIENTS
        ICODE=8
        CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VSRE,VKE,VKE1,VIE,VIE1,VFE,
     1            VSRE,VNDSRC,KTYPRSG,VPRSG)
        ITPE1=ITPE
      ENDDO
      
      RETURN
      END

