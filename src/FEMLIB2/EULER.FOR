      SUBROUTINE EULER(VCORG,KDLNC,VDIMP,KNEQ,KLD,KLOCE,VCORE,VPRNE,   
     1  VPREE,KNE,VKE,VME,VFE,VDLE,VKGS,VKGD,VKGI,VFG,VRES,VDLG,
     2  VDLE0,VDLG0,VFG0)
C**********************************************************************C
C EULER.FOR                                                            C
C  02/12-1994   NTH                                                    C
C  Bruno Brouard                                                       C
C**********************************************************************C
C     ALGORITHM FOR EULER TYPE METHODS (IMPLICIT, EXPLICIT OR BOTH     C
C     ACCORDING TO OMEGA) FOR LINEAR OR NON LINEAR PROBLEMS.           C
C     THE NON LINEAR PROBLEM IS SOLVED BY A NEWTON-RAPHSON             C
C     METHOD                                                           C
C       IMETH.EQ.1   STANDARD NEWTON-RAPHSON                           C
C       IMETH.EQ.2   K IS CONSTANT                                     C
C       IMETH.EQ.3   K IS RECOMPUTED AT THE BEGINING OF EACH STEP      C
C**********************************************************************C
C  Programmed after Dhatt & Touzot, p.433                              C
C**********************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ASSE/NSYM,NKG
      COMMON/RESO/NEQ
      COMMON/NLIN/EPSDL,XNORM,OMEGA,XPAS,DPAS,DPAS0,NPAS,IPAS,NITER,
     1  ITER,IMETH
      COMMON/ES/M,MR,MP
      DIMENSION VCORG(1),KDLNC(1),VDIMP(1),KNEQ(1),KLD(1),KLOCE(1),
     1  VCORE(1),VPRNE(1),VPREE(1),KNE(1),VKE(1),VME(1),VFE(1),
     2  VDLE(1),VKGS(1),VKGD(1),VKGI(1),VFG(1),VRES(1),VDLG(1),
     3  VDLE0(1),VDLG0(1),VFG0(1)
      DATA ZERO/0.D0/,UN/1.D0/
C----------------------------------------------------------------------- 
C--------    DECIDE IF GLOBAL MATRIX IS TO BE REASSEMBLED
      IKT=0
      IF(IMETH.EQ.1) GO TO 10
      IF (DPAS.NE.DPAS0.AND.ITER.EQ.1) GO TO 10
      IF (IMETH.EQ.3.AND.ITER.EQ.1) GO TO 10
      GO TO 20
 10   IKT=1
C---------   INITIALIZE GLOBAL MATRIX TO ZERO IF NECESSARY
 20   IF(IKT.EQ.0) GO TO 30
      CALL INIT(ZERO,NKG,VKGS)
      CALL INIT(ZERO,NEQ,VKGD)
      IF(NSYM.EQ.1) CALL INIT(ZERO,NKG,VKGI)
C--------    ASSEMBLE RESIDUALS AND GLOBAL MATRIX IF REQUIRED
 30   CALL MAJ(UN,ZERO,NEQ,VFG0,VRES)
      CALL ASEULR(IKT,VCORG,KDLNC,VDIMP,KNEQ,KLD,KLOCE,VCORE,VPRNE,   
     1  VPREE,KNE,VKE,VME,VFE,VDLE,VKGS,VKGD,VKGI,VFG,VRES,VDLG,
     2  VDLE0,VDLG0,VFG0)
      C1=UN
      IF(ITER.GT.1) C1=C1-OMEGA
      DO 40 I=1,NEQ
 40   VRES(I)=DPAS*(VRES(I)-C1*VFG(I))      
C--------      SOLVE
      CALL SOL(VKGS,VKGD,VKGI,VRES,KLD,NEQ,MP,IKT,1,NSYM,ENERG)
C--------      UPDATE THE SOLUTION       
      CALL MAJ(UN,UN,NEQ,VRES,VDLG)
      RETURN
      END
      