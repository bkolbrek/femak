      SUBROUTINE ASSEL(IKG,IFG,IDLE,NSYM,KLOCE,KLD,VKE,VFE,VKGS,
     1  VKGD,VKGI,VFG)
C***************************************************************
C     ASSEL.FOR
C***************************************************************
C     24/2 1992
C     Ulf R. Kristiansen, NTH
C     Programmed after Dhatt og Touzot p.413
C*************************************************************** 
C     TO ASSEMBLE AN ELEMENT MATRIX AND/OR VECTOR
C     (MATRIX SYMMETRICAL OR NOT)
C       INPUT
C          IKG     IF IKG.EQ.1 ASSEMBLE ELEMENT MATRIX KE
C          IFG     IF IFG.EQ.1 ASSEMBLE ELEMENT VECTOR FE
C          IDLE    ELEMENT NUMBER D.O.F.
C          NSYM    0=SYMMETRIC PROBLEM, 1=UNSYMMETRIC PROBLEM
C          KLOCE   ELEMENT LOCALIZATION VECTOR
C          KLD     CUMULATIVE COLUMN HEIGHTS OF KG
C          VKE     ELEMENT MATRIX KE (FULL OR UPPER TRIANGLE BY
C                  DESCENDING COLUMNS)
C          VFE     ELEMENT VECTOR FE
C       OUTPUT
C          VKGS,VKGD,VKGI    GLOBAL MATRIX (SKYLINES)
C                  (SYMMETRIC OR NOT)
C          VFG     GLOBAL LOAD VECTOR
C===============================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION KLOCE(1),KLD(1),VKE(1),VFE(1),VKGS(1),VKGD(1),
     1  VKGI(1),VFG(1)
C---------------------------------------------------------------
C
C------  ASSEMBLE ELEMENT MATRIX
C
      IF(IKG.NE.1) GO TO 100
      IEQ0=IDLE
      IEQ1=1
C------  FOR EACH COLUMN OF KE
      DO 90 JD=1,IDLE
      IF(NSYM.NE.1) IEQ0=JD
      JL=KLOCE(JD)
      IF(JL) 90,90,10
10    I0=KLD(JL+1)
      IEQ=IEQ1
      IQ=1
C------  FOR EACH ROW OF KE
      DO 80 ID=1,IDLE
      IL=KLOCE(ID)
      IF(NSYM.EQ.1) GO TO 30
      IF(ID-JD) 30,20,20
20    IQ=ID
30    IF(IL) 80,80,40
40    IJ=JL-IL
      IF(IJ) 70,50,60
C------  DIAGONAL TERMS OF KG
50    VKGD(IL)=VKGD(IL)+VKE(IEQ)
      GO TO 80
C------  UPPER TRIANGLE TERMS OF KG
60    I=I0-IJ
      VKGS(I)=VKGS(I)+VKE(IEQ)
      GO TO 80
C------  LOWER TRIANGLE TERMS OF KG
70    IF(NSYM.NE.1) GO TO 80 
      I=KLD(IL+1)+IJ
      VKGI(I)=VKGI(I)+VKE(IEQ)
80    IEQ=IEQ+IQ
90    IEQ1=IEQ1+IEQ0
C
C------  ASSEMBLE ELEMNT LOAD VECTOR
C
100   IF(IFG.NE.1) GO TO 130
      DO 120 ID=1,IDLE
      IL=KLOCE(ID)
      IF(IL) 120,120,110
110   VFG(IL)=VFG(IL)+VFE(ID)
120   CONTINUE
130   RETURN
      END  
