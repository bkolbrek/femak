      SUBROUTINE ASEULR(IKT,VCORG,KDLNC,VDIMP,KNEQ,KLD,KLOCE,VCORE,
     1  VPRNE,VPREE,KNE,VKE,VME,VFE,VDLE,VKGS,VKGD,VKGI,VFG,VRES,
     2  VDLG,VDLE0,VDLG0,VFG0)
C**********************************************************************C
C     ASEULR.FOR                                                       C
C  19/1-1994   NTH                                                     C
C  Bruno Brouard                                                       C
C**********************************************************************C
C     TO ASSEMBLE THE RESIDUALS AND THE GLOBAL MATRIX (IF IKT.EQ.1)    C
C     WHILE LOOPING OVER THE ELEMENTS (FOR EULER METHOD)               C
C**********************************************************************C
C  Programmed after Dhatt & Touzot, p.433-434                          C
C**********************************************************************C
C  Changed for new version of FEMAK. Note this version is not clear    C
C  for nonsymmetrical version yet. (VKE1 is not among the input        C
C  parameters.) Neither are the following variables (which are only    C
C  referenced as simple variables, (we don't need them (we think)))    C
c      VDLEC,VIE,VIE1,VSRE,VNDSRC,KTYPRSG,VPRSG                        C
C**********************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)

      INCLUDE 'FEMAK.CMN'


      DIMENSION VCORG(1),KDLNC(1),VDIMP(1),KNEQ(1),KLD(1),KLOCE(1),
     1  VCORE(1),VPRNE(1),VPREE(1),KNE(1),VKE(1),VME(1),VFE(1),VDLE(1),
     2  VKGS(1),VKGD(1),VKGI(1),VFG(1),VRES(1),VDLG(1),VDLE0(1),
     3  VDLG0(1),VFG0(1)

      COMPLEX*16 VDLEC,VIE,VIE1,VSRE,VNDSRC,VPRSG
      INTEGER KTYPRSG

      DATA UN/1.D0/
C-----------------------------------------------------------------------
      CC=DPAS*OMEGA
      IFE=0
      IF(ITER.GT.1) IFE=1
C--------   REWIND ELEMENT FILE (M1)
      REWIND M1
      ITPE1=0
C--------    LOOP OVER THE ELEMENTS
      DO 90 IE=1,NELT
C--------    READ AN ELEMENT
      CALL RDELEM(M1,KLOCE,VCORE,VPRNE,VPREE,KNE)
C--------    EVALUATE INTERPOLATION FUNCTIONS IF REQUIRED
      IF(ITPE.EQ.ITPE1) GO TO 10
      ICODE=2
      CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
c      CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VKE,VFE)
C---------   FIND ELEMENT D.O.F. FROM VFG
 10   CALL DLELM(KLOCE,VDLG,VDIMP,VDLE)
C---------  COMPUTE THE RESIDUAL K.U.
      ICODE=6
      CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
c      CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VKE,VFE)
C---------   COMPUTE MATRIX M
      ICODE=5
      CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VME,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
c      CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VME,VFE)
C---------   COMPUTE MATRIX K IF REQUIRED
      IF(IKT.EQ.0) GO TO 15
      ICODE=3
      CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
c      CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VKE,VFE)
C---------    RESIDUALS OF THE FIRST ITERATION IN EACH STEP (LINEAR)
 15   IF(ITER.GT.1) GO TO 20
      CALL ASSEL(0,1,IDLE,NSYM,KLOCE,KLD,VKE,VFE,VKGS,VKGD,VKGI,VFG)
      GO TO 60
C--------      RESIDUALS AFTER FIRST ITERATION
 20   CALL DLELM(KLOCE,VDLG0,VDIMP,VDLE0)
      DO 30 I=1,IDLE
      VDLE(I)=(VDLE0(I)-VDLE(I))/DPAS
 30   VFE(I)=-OMEGA*VFE(I)
C--------       PRODUCT  M . U
      VFE(1)=VFE(1)+VME(1)*VDLE(1)
      II=1
      DO 50 J=2,IDLE
      J1=J-1
      DO 40 I=1,J1
      II=II+1
      VFE(I)=VFE(I)+VME(II)*VDLE(J)
 40   VFE(J)=VFE(J)+VME(II)*VDLE(I)
      II=II+1
 50   VFE(J)=VFE(J)+VME(II)*VDLE(J)
C---------     MATRIX    M + DPAS.OMEGA.K
 60   IF(IKT.EQ.0) GO TO 80
      II=0
      DO 70 I=1,IDLE
      DO 70 J=I,IDLE
      II=II+1
 70   VKE(II)=VKE(II)*CC+VME(II)
C--------     ASSEMBLE THE RESIDUAL AND GLOBAL MATRIX
 80   CALL ASSEL(IKT,IFE,IDLE,NSYM,KLOCE,KLD,VKE,VFE,VKGS,VKGD,VKGI,
     1  VRES)
 90   ITPE1=ITPE
      RETURN
      END
