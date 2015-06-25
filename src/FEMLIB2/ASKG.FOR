      SUBROUTINE ASKG(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1   VIE,VIE1,VFE,VKGS,VKGD,VKGI,VFG,VDLE,VRES,
     2   VNDSRC,KTYPRSG,VPRSG)
C*******************************************************************
C     ASKG.FOR
C*******************************************************************
C     2/1-1995
C     NTH
C     Ulf R. Kristiansen & Tonni F.Johansen
C     Programmert etter Dhatt og Touzot s.410 
C*******************************************************************
C     TO ASSEMBLE GLOBAL MATRIX KG (ELEMENT FUNCTION 3)
C     TAKING INTO ACCOUNT OF NON ZERO PRESCRIBED D.O.F.
C*******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include 'femak.cmn'

      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VDLE(1),VKE(1),VKE1(1),VFE(1)
      COMPLEX*16 VIE(1),VIE1(1),VSRE(1),VNDSRC(1),VPRSG(1)
      INTEGER KTYPRSG(1)

      INTEGER KLD(1),KLOCE(1),KNE(1)
      REAL*8 VDIMP(1),VKGS(1),VKGD(1),VKGI(1),VFG(1),VRES(1)


      
      ITPE1=0
C-----------------------------------------------------------------------
C------ REWIND ELEMENT FILE FOR ELEMENTS WITH REAL COEFFICIENTS (M1)
      REWIND M1
C------ LOOP OVER THE ELEMENTS
      DO IE=1,NELTR

C------ READ AN ELEMENT ON FILE M1
        CALL RDELEM(M1,KLOCE,VCORE,VPRNE,VPREE,KNE)

C--------Elementtype .ge. 100 means special coupling element that will
C--------not contribute to mass or stiffness matrices.
        IF (ITPE.LT.100) THEN

C------ SKIP COMPUTATION IF IDENTICAL ELEMENTS ENCOUNTERED
          IF(NIDENT.EQ.1.AND.IE.GT.1) GO TO 20
C------ EVALUATE INTERPOLATION FUNCTIONS IF REQUIRED
          IF(ITPE.EQ.ITPE1) GO TO 10
          ICODE=2

          CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C------ FORM ELEMENT MATRIX
10        ICODE=3
C          IF ((NSYM.EQ.0).OR.(ITPE.EQ.2)) THEN
            CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C C         ELSE
C            CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VKE,VKE1,VIE,VIE1,VFE,
C     1            VSRE,VNDSRC,KTYPRSG,VPRSG)

C            PNTCO=-IDLE+1
C            PNTCO1=1
C            DO COL=1,IDLE
C              PNTCO1=PNTCO1+COL-1
C              PNTCO=PNTCO+IDLE
C              PNTRW1=0
C              DO ROW=1,IDLE
C                IF (COL.GE.ROW) THEN
C                  VKE(PNTCO+ROW-1)=VKE1(PNTCO1+ROW-1)
C                  PNTRW1=PNTRW1+1
C                ELSE
C                  PNTRW1=PNTRW1+ROW-1
C                  VKE(PNTCO+ROW-1)=VKE1(PNTCO1+PNTRW1-1)         
C                END IF
C              END DO
C            END DO
C          ENDIF

C------ PRINT ELEMENT MATRIX
        IF(M.LT.2) GO TO 20
        IF(NSYM.EQ.0) IKE=IDLE*(IDLE+1)/2
        IF(NSYM.EQ.1) IKE=IDLE*IDLE
        WRITE(MP,2000) IEL,(VKE(I),I=1,IKE)
2000    FORMAT(/' MATRIX (KE) , ELEMENT:',I5/(10X,10E12.5))
C------  MODIFY FG FOR NON ZERO PRESCRIBED D.O.F.
20      IF(NCLNZ.NE.0) CALL MODFG(IDLE,NSYM,KLOCE,VDIMP,VKE,VFG)
C------ ASSEMBLE
        CALL ASSEL(1,0,IDLE,NSYM,KLOCE,KLD,VKE,VFE,VKGS,VKGD,VKGI,VFG)

C----------M
        ENDIF
C----------Finish test (ITPE.LT.100)

        ITPE1=ITPE
      ENDDO 
      
      RETURN
      END
