      SUBROUTINE ASFG(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1  VIE,VIE1,VFE,VKGS,VKGD,VKGI,VFG,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG)
C*******************************************************************
C     ASFG.FOR
C*******************************************************************
C     2/1-1995
C     Tonni F. Johansen,   NTH
C     Programmert etter Dhatt og Touzot s.405
C*******************************************************************
C     TO ASSEMBLE DISTRIBUTE LOADS IN FG
C*******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include 'femak.cmn'

      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VKE1(1),VFE(1)
      COMPLEX*16 VIE(1),VIE1(1),VSRE(1),VNDSRC(1),VPRSG(1)
      INTEGER KTYPRSG(1)

      INTEGER KLD(1),KLOCE(1),KNE(1)
      REAL*8 VDIMP(1),VKGS(1),VKGD(1),VKGI(1),VFG(1),VRES(1)
      
      COMPLEX*16 VDLEC(1)
      
C-----------------------------------------------------------------------
C------ REWIND ELEMENT FILE FOR ELEMENTS WITH REAL COEFFICIENTS (M1)
      REWIND M1
      ITPE1=0
C------ LOOP OVER THE ELEMENTS
      DO IE=1,NELT

C------ READ AN ELEMENT ON FILE M1
        CALL RDELEM(M1,KLOCE,VCORE,VPRNE,VPREE,KNE)

C------ EVALUATE INTERPOLATION FUNCTIONS IF REQUIRED
        IF(ITPE.EQ.ITPE1) GO TO 10
        ICODE=2

        CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,VFE,
     1   VSRE,VNDSRC,KTYPRSG,VPRSG)

C------ EVALUATE ELEMENT VECTOR
10      ICODE=7
        CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1   VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)

C------ PRINT ELEMENT VECTOR VFE
        IF(M.GE.2) WRITE(MP,2000) IEL,(VFE(I),I=1,IDLE)
2000  FORMAT(/' VECTOR (FE) , ELEMENT:',I5/(10X,10E12.5))

C------ ASSEMBLE
        CALL ASSEL(0,1,IDLE,NSYM,KLOCE,KLD,VKE,VFE,VKGS,VKGD,VKGI,VFG)
        ITPE1=ITPE
      ENDDO

      RETURN
      END
