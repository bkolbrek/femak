      SUBROUTINE ASLGC(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1  VIE,VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VSRC,VDLE,VRES,VNDSRC,
     2  KTYPRSG,VPRSG,VFREQ)
C*******************************************************************
C     ASLGC.FOR
C*******************************************************************
C     3/1-1995
C     Tonni F. Johansen   NTH
C*******************************************************************
C     TO ASSEMBLE GLOBAL LOSS-MATRIX  (ELEMENT FUNCTION 11)
C     TAKING INTO ACCOUNT OF NON ZERO PRESCRIBED D.O.F.
C*******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include 'femak.cmn'
C$DECLARE

      INTEGER KLD(1),KLOCE(1),KNE(1),KEB(1)
      REAL*8 VDIMP(1),VCORE(1),VPRNE(1),VPREE(1),VKE(1),VKE1(1),VFE(1)
      REAL*8 VKGS(1),VKGD(1),VKGI(1),VFG(1),VDLE(1),VRES(1)
      COMPLEX*16 VIGS(1),VIGD(1),VIGI(1),VIE(1),VIE1(1),VSRC(1)

      COMPLEX*16 VNDSRC(1),VPRSG(1)
      INTEGER KTYPRSG(1)

      INTEGER IE,IKE,I



C-----------------------------------------------------------------------
C------ REWIND ELEMENT FILE (M2)
      REWIND M2
      ITPE1=0
C------ LOOP OVER THE ELEMENTS
      DO IE=1,NELTC

C------ READ AN ELEMENT ON FILE M2
        CALL RDELEM(M2,KLOCE,VCORE,VPRNE,VPREE,KNE)

C--------Loss matrix only for ELEM05
        IF (ITPE.EQ.5) THEN

C------ SKIP COMPUTATION IF IDENTICAL ELEMENTS ENCOUNTERED
          IF(NIDENT.EQ.1.AND.IE.GT.1) GO TO 20
C------ EVALUATE INTERPOLATION FUNCTIONS IF REQUIRED
          IF(ITPE.EQ.ITPE1) GO TO 10
          ICODE=2

          CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C------ FORM ELEMENT MATRIX
10        ICODE=11
C          IF ((NSYM.EQ.0).OR.(ITPE.EQ.11).OR.(ITPE.EQ.92)) THEN
            CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C          ELSE
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
C                  VIE(PNTCO+ROW-1)=VIE1(PNTCO1+ROW-1)
C                  PNTRW1=PNTRW1+1
C                ELSE
C                  PNTRW1=PNTRW1+ROW-1
C                  VIE(PNTCO+ROW-1)=VIE1(PNTCO1+PNTRW1-1)         
C                END IF
C              END DO
C            END DO

C          ENDIF
C------ PRINT ELEMENT MATRIX
          IF(M.LT.2) GO TO 20
          IF(NSYM.EQ.0) IKE=IDLE*(IDLE+1)/2
          IF(NSYM.EQ.1) IKE=IDLE*IDLE
          WRITE(MP,2000) IEL,(VIE(I),I=1,IKE)
2000      FORMAT(/' MATRIX (IE - LOSS) , ELEMENT:',I5/(10X,10E12.5))
C------  MODIFY FG FOR NON ZERO PRESCRIBED D.O.F.
20        IF(NCLNZ.NE.0) CALL MODFG(IDLE,NSYM,KLOCE,VDIMP,VKE,VFG)
C------ ASSEMBLE
          CALL ASSELC(1,0,IDLE,NSYM,KLOCE,KLD,VIE,VSRE,VIGS,VIGD,VIGI,
     1   VSRC)
        ENDIF
C----------Finish test (ITPE.LT.100)

        ITPE1=ITPE
      ENDDO
      
      RETURN
      END
