      SUBROUTINE ASMGC(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1  VIE,VIE1,VFE,VSRE,VKGS,VKGD,VKGI,VIGS,VIGD,VIGI,VFG,VSRC,VDLE,
     2   VRES,VNDSRC,KTYPRSG,VPRSG)
C*******************************************************************
C     ASMGC.FOR
C*******************************************************************
C     2/1-1995
C     Ulf R. Kristiansen & Tonni F. Johansen  NTH
C     Complex version of ASMG as described in Dhatt og Touzot s. 439
C*******************************************************************
C     TO ASSEMBLE THE GLOBAL MASS MATRIX (ELEMENT FUNCTION 5)
C===================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      include 'femak.cmn'

C $DECLARE

      INTEGER KLD(1),KLOCE(1),KNE(1),KEB(1)
      REAL*8 VDIMP(1),VCORE(1),VPRNE(1),VPREE(1),VKE(1),VKE1(1),VFE(1)
      REAL*8 VKGS(1),VKGD(1),VKGI(1),VFG(1),VDLE(1),VRES(1)
      COMPLEX*16 VIGS(1),VIGD(1),VIGI(1),VIE(1),VIE1(1),VSRC(1),vsre(1)

      COMPLEX*16 VNDSRC(1),VPRSG(1)
      INTEGER KTYPRSG(1)
      
      INTEGER I,IE,IKE

C-----------------------------------------------------------------------
C------ REWIND ELEMENT FILE FOR ELEMENTS WITH COMPLEX COEFFICIENTS (M2)
      REWIND M2
      ITPE1=0
C------  LOOP OVER THE ELEMENTS
      DO IE=1,NELTC
C------ READ AN ELEMENT ON FILE M2
        CALL RDELEM(M2,KLOCE,VCORE,VPRNE,VPREE,KNE)

C--------Elementtype .ge. 100 means special coupling element that will
C--------not contribute to mass or stiffness matrices.
        IF (ITPE.LT.100) THEN

C------ SKIP COMPUTATION IF IDENTICAL ELEMENTS ENCOUNTERED
          IF(NIDENT.EQ.1.AND.IE.GT.1) GO TO 20

          IF(ITPE.EQ.ITPE1) GO TO 10
          ICODE=2
          CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
10        ICODE=5

cccccc!!!!! midlertidig !!!!!!!!

c        IF ((NSYM.EQ.0).OR.(ITPE.EQ.11).OR.(ITPE.EQ.92).OR.
c     1                                       (ITPE.EQ.93)) THEN
            CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
          
c          ELSE
c            CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VKE,VKE1,VIE,VIE1,VFE,
c     1            VSRE,VNDSRC,KTYPRSG,VPRSG)
c            PNTCO=-IDLE+1
c            PNTCO1=1
c            DO COL=1,IDLE
c              PNTCO1=PNTCO1+COL-1
c              PNTCO=PNTCO+IDLE
c              PNTRW1=0
c              DO ROW=1,IDLE
c                IF (COL.GE.ROW) THEN
c                  VIE(PNTCO+ROW-1)=VIE1(PNTCO1+ROW-1)
c                  PNTRW1=PNTRW1+1
c                ELSE
c                  PNTRW1=PNTRW1+ROW-1
c                  VIE(PNTCO+ROW-1)=VIE1(PNTCO1+PNTRW1-1)         
c                END IF
c              END DO
c            END DO

c          ENDIF
C------ PRINT ELEMENT MATRIX
          IF(M.LT.2) GO TO 20
          IF(NSYM.EQ.0) IKE=IDLE*(IDLE+1)/2
          IF(NSYM.EQ.1) IKE=IDLE*IDLE
          WRITE(MP,2000) IEL,(VIE(I),I=1,IKE)
2000      FORMAT(/' MATRIX (IE - MASS) , ELEMENT:',I5/(10X,10E12.5))
C------ ASSEMBLE
20        CALL ASSELC(1,0,IDLE,NSYM,KLOCE,KLD,VIE,VSRE,VIGS,VIGD,VIGI,
     1        VSRC)

        ENDIF
C----------Finish test (ITPE.LT.100)

        ITPE1=ITPE
      ENDDO


      RETURN
      END
