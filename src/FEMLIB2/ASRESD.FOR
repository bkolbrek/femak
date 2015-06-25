      SUBROUTINE ASRESD(IRESD,IREAC,KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,
     1  KNE,VKE,VKE1,VIE,VIE1,VFE,VSRE,VKGS,VKGD,VKGI,VFG,VDLE,VRES,
     2  VNDSRC,KTYPRSG,VPRSG,VREAC)
C*********************************************************************
C     ASRESD.FOR
C*********************************************************************
C     3/1 1995
C     NTH
C     Ulf R. Kristiansen
C     Programmed after Dhatt og Touzot s. 411
C*********************************************************************
C     TO ASSEMBLE INTERNAL RESIDUALS IN VRES (IF IRESD .EQ. 1)
C     AND EXTERNAL REACTIONS IN VREAC (IF IREAC.EQ.1)
C=====================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      include 'femak.cmn'
C  $DECLARE

      INTEGER IRESD,IREAC
      COMPLEX*16 VIE(1),VIE1(1),VSRE(1)
      REAL*8 VDIMP(1),VCORE(1),VPRNE(1),VPREE(1),VKE(1),VFE(1),VDLE(1)
      REAL*8 VRES(1),VKE1(1),VKGS(1),VKGD(1),VKGI(1),VREAC(1)
      INTEGER KLD(1),KLOCE(1),KNE(1)

      COMPLEX*16 VNDSRC(1),VPRSG(1)
      INTEGER KTYPRSG(1)
      
      INTEGER I,ID,IE
C-----------------------------------------------------------------------
C------------  REWIND ELEMENT FILE (M1)
      REWIND M1
      ITPE1=0
C------------  LOOP OVER THE ELEMENTS
      DO IE=1,NELT
C------------  READ AN ELEMENT ON FILE M2
        CALL RDELEM(M1,KLOCE,VCORE,VPRNE,VPREE,KNE)
C------------  EVALUATE INTERPOLATION FUNCTION IF REQUIRED
        IF(ITPE.NE.ITPE1) THEN
          ICODE=2
          CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)                       
        ENDIF
C-----------   FIND ELEMENT D.O.F.
        CALL DLELM(KLOCE,VFG,VDIMP,VDLE)
C----------    EVALUATE ELEMENT REACTIONS
        ICODE=6
        CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C----------    PRINT ELEMENT REACTIONS
        IF(M.GE.2) WRITE(MP,2000) IEL,(VFE(I),I=1,IDLE)
2000    FORMAT(/' REACTIONS (FE) , ELEMENT:',I5,/(10X,10E12.5))
        IF(IRESD.EQ.1) THEN
C----------  ASSEMBLE INTERNAL RESIDUALS
          CALL ASSEL(0,1,IDLE,NSYM,KLOCE,KLD,VKE,VFE,VKGS,VKGD,VKGI,
     1        VRES)
        ENDIF
        IF(IREAC.EQ.1) THEN
C----------  ASSEMBLE EXTERNAL REACTIONS
C          MODIFY TERMS IN CLOCE SUCH THAT PRESCRIBED D.O.F. ARE THE ONLY
C          ASSEMBLED ONES
          DO ID=1,IDLE
            IF(KLOCE(ID).LT.0) THEN
              KLOCE(ID)=-KLOCE(ID)
            ELSEIF(KLOCE(ID).GT.0) THEN
              KLOCE(ID)=0      
            ENDIF
          ENDDO
          CALL ASSEL(0,1,IDLE,NSYM,KLOCE,KLD,VKE,VFE,VKGS,VKGD,VKGI,
     1         VREAC)
        ENDIF

        ITPE1=ITPE
      ENDDO
      
      RETURN
      END

