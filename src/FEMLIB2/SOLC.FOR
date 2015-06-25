      Subroutine SOLC(VKGS,VKGD,VKGI,VFG,KLD,NEQ,MP,IFAC,ISOL,NSYM,
     1         ENERG)
C***********************************************************************
C     SOL.FOR
C***********************************************************************
C     25/2 1992
C     Kompaniet Numerisk Lyd
C     Ulf R. Kristiansen
C     Programmert etter  Dhatt og Touzot s.284
C***********************************************************************
C     TO SOLVE A LINEAR SYSTEM (SYMMETRICAL OR NOT).
C     THE MATRIX IS STORED IN CORE BY SKYLINES IN ARRAYS
C     VKGS,VKGD,VKGI
C       INPUT
C          VKGS,VKGD,VKGI   SYSTEM MATRIX : UPPER, DIAGONAL AND
C                           LOWER PARTS
C          VFG              SECOND MEMBER
C          KLD              ADRESSES OF COLUMN TOP TERMS
C          NEQ              NUMBER OF EQUATIONS
C          MP               OUTPUT DEVICE NUMBER
C          IFAC             IF IFAC.EQ.1 TRIANGULARIZE THE MATRIX
C          ISOL             IF ISOL.EQ.1 COMPUTE THE SOLUTION FROM
C                           TRIANGULARIZED MATRIX
C          NSYM             INDEX FOR NONSYMMETRIC PROBLEM
C       OUTPUT
C          VKGS,VKGD,VKGI   TRIANGULARIZED MATRIX (IF IFAC.EQ.1)
C          VFG              SOLUTION (IF ISOL.EQ.1)
C          ENERG            SYSTEM ENERGY (IF NSYM.EQ.0)
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)

      INTEGER KLD(1)
      COMPLEX*16 SCALC
      COMPLEX*16 VKGS(1),VKGD(1),VKGI(1),VFG(1)
      COMPLEX*16 ZERO,C1,C2,CDIAG,ENERG
      DATA ZERO/0.0D0/
C-----------------------------------------------------------------------
      IK=1
c........prøver med ny størrelse på sjekk for "PIVOT"      
c      IF(cabs(VKGD(1)).gt.0.0000001) GO TO 10
      IF(abs(VKGD(1)).gt.0.0000001) GO TO 10
      WRITE(MP,2000) IK  
      STOP
10    ENERG=ZERO
C
C------  FOR EACH COLUMN IK TO BE MODIFIED
C
      JHK=1
      DO 100 IK=2,NEQ
C------  ADRESS OF THE NEXT COLUMN TOP TERM IK+1
      JHK1=KLD(IK+1)
C------  HEIGHT OF COLUMN IK (INCLUDE UPPER AND DIAGONAL TERMS)
      LHK=JHK1-JHK
      LHK1=LHK-1
C------  ROW OF FIRST TERM TO BE MODIFIED IN COLUMN IK
      IMIN=IK-LHK1
      IMIN1=IMIN-1
C------  ROW OF LAST TERM TO BE MODIFIED IN COLUMN IK
      IMAX=IK-1
      IF(LHK1.LT.0) GO TO 100
      IF(IFAC.NE.1) GO TO 90
      IF(NSYM.EQ.1) VKGI(JHK)=VKGI(JHK)/VKGD(IMIN1)
      IF(LHK1.EQ.0) GO TO 40
C
C------  MODIFY NON-DIAGONAL TERM IN COLUMN IK
C
      JCK=JHK+1
      JHJ=KLD(IMIN)
C------  FOR EACH TERM LOCATED AT JCK AND CORRESPONDING TO COLUMN IJ
      DO 30 IJ=IMIN,IMAX
      JHJ1=KLD(IJ+1)
C------  NUMBER OF MODIFICATIVE TERMS FOR COEFFICIENT LOCATED AT JCK
      IC=MIN0(JCK-JHK,JHJ1-JHJ)
      IF(IC.LE.0.AND.NSYM.EQ.0) GO TO 20
      C1=ZERO
      IF(IC.LE.0) GO TO 17
      J1=JHJ1-IC
      J2=JCK-IC
      IF(NSYM.EQ.1) GO TO 15
      VKGS(JCK)=VKGS(JCK)-SCALC(VKGS(J1),VKGS(J2),IC)
      GO TO 20
 15   VKGS(JCK)=VKGS(JCK)-SCALC(VKGI(J1),VKGS(J2),IC)
      C1=SCALC(VKGS(J1),VKGI(J2),IC)
17    VKGI(JCK)=(VKGI(JCK)-C1)/VKGD(IJ)
20    JCK=JCK+1
30    JHJ=JHJ1
C
C------  MODIFY DIAGONAL TERM 
C
40    JCK=JHK
      CDIAG=ZERO
      DO 70 IJ=IMIN1,IMAX
      C1=VKGS(JCK)
      IF(NSYM.EQ.1) GO TO 50
      C2=C1/VKGD(IJ)
      VKGS(JCK)=C2
      GO TO 60
50    C2=VKGI(JCK)
60    CDIAG=CDIAG+C1*C2
70    JCK=JCK+1
      VKGD(IK)=VKGD(IK)-CDIAG
      IF(REAL(VKGD(IK))) 90,80,90
80    WRITE(MP,2000) IK
2000  FORMAT(' *** ERROR,ZERO PIVOT EQUATION ',I5)
      STOP
C
C------  SOLVE LOWER TRIANGULAR SYSTEM
C
90    IF(ISOL.NE.1) GO TO 100
      IF(NSYM.NE.1) VFG(IK)=VFG(IK)-SCALC(VKGS(JHK),VFG(IMIN1),LHK)
      IF(NSYM.EQ.1) VFG(IK)=VFG(IK)-SCALC(VKGI(JHK),VFG(IMIN1),LHK)
100   JHK=JHK1
      IF(ISOL.NE.1) RETURN
C
C------  SOLVE DIAGONAL SYSTEM
C
      IF(NSYM.EQ.1) GO TO 120
      DO 110 IK=1,NEQ
      C1=VKGD(IK)
      C2=VFG(IK)/C1
      VFG(IK)=C2
110   ENERG=ENERG+C1*C2*C2
C
C------  SOLVE DIAGONAL SYSTEM
C
120   IK=NEQ+1
      JHK1=KLD(IK)
130   IK=IK-1
      IF(NSYM.EQ.1) VFG(IK)=VFG(IK)/VKGD(IK)
      IF(IK.EQ.1) RETURN
      C1=VFG(IK)
      JHK=KLD(IK)
      JBK=JHK1-1
      IF(JHK.GT.JBK) GO TO 150
      IJ=IK-JBK+JHK-1
      DO 140 JCK=JHK,JBK
      VFG(IJ)=VFG(IJ)-VKGS(JCK)*C1
140   IJ=IJ+1
150   JHK1=JHK
      GO TO 130
      END
      
      COMPLEX*16 FUNCTION SCALC(X,Y,N)
C===============================================================
C     INNER PRODUCT OF VECTORS X AND Y OF LENGTH N
C       (FUNCTION TO BE WRITTEN EVENTUALLY IN ASSEMBLER)
C===============================================================
      COMPLEX*16 X(1),Y(1)
      COMPLEX*16 ZERO
C---------------------------------------------------------------
      ZERO=DCMPLX(0.D0,0.D0)
      SCALC=ZERO
      DO 10 I=1,N
10    SCALC=SCALC+X(I)*Y(I)
      RETURN
      END

      
