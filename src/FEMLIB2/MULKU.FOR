      SUBROUTINE MULKU(VKGS,VKGD,VKGI,KLD,VFG,NEQ,NSYM,VRES)
C************************************************************
C     MULKU.FOR
C************************************************************
C     25/2/1992
C     Kompaniet Numerisk Lyd
C     Ulf R. Kristiansen
C     Programmert etter Dhatt og Touzot s.417
C************************************************************
C     SUBPROGRAM :       MULK   3
C     TO ADD VECTOR RES TO THE PRODUCT OF MATRIX KG AND THE VECTOR FG
C       INPUT
C          VKGS,VKGD,VKGI  MATRIX KG STORED BY SKYLINE    MULK   6
C                          (SYM. OR NON SYM.)
C          KLD     ARRAY OF ADRESS OF COLUMN TOP TERMS IN KG
C          VFG     VECTOR FG
C          NEQ     ORDER OF VECTORS FG AND RES
C          NSYM    .EQ.1 IF NON SYMMETRIC PROBLEM
C          VRES    VECTOR RES
C       OUTPUT
C          VRES    VECTOR RES
C====================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VKGS(1),VKGD(1),VKGI(1),KLD(1),VFG(1),VRES(1)
C--------------------------------------------------------------------
C------  FOR EACH COLUMN OF MATRIX KG
      DO 20 IK=1,NEQ
      JHK=KLD(IK)
      JHK1=KLD(IK+1)
      LHK=JHK1-JHK
C------  DIAGONAL TERMS
      C=VKGD(IK)*VFG(IK)
      IF(LHK.LE.0) GO TO 20
      I0=IK-LHK
C------  ROW TWRMS
      IF(NSYM.NE.1) C=C+SCAL(VKGS(JHK),VFG(I0),LHK)
      IF(NSYM.EQ.1) C=C+SCAL(VKGI(JHK),VFG(I0),LHK)
C------  COLUMN TERMS
      J=JHK
      I1=IK-1
      DO 10 IJ=I0,I1
      VRES(IJ)=VRES(IJ)+VKGS(J)*VFG(IK)
10    J=J+1
20    VRES(IK)=VRES(IK)+C
      RETURN
      END
