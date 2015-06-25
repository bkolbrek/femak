      SUBROUTINE PRRESD(VKGS,VKGD,VKGI,VFG,KLD,VRES)
C******************************************************************
C     PRRESD.FOR
C******************************************************************
C     25/2 1992
C     Ulf R. Kristiansen   NTH
C     programmed after Dhatt and Touzot p.415
C******************************************************************
C     TO COMPUTE AND PRINT THE RESIDUAL VECTOR K.U-F
C==================================================================
      IMPLICIT REAL*8(A-H,O-Z)

      include 'femak.cmn'    

      DIMENSION VKGS(1),VKGD(1),VKGI(1),VFG(1),KLD(1),VRES(1)
      DATA ZERO/0.D0/
      ABS(X)=DABS(X)
C------------------------------------------------------------------
      REWIND M3
C------  SKIP VECTOR FG NON MODIFIED BY B.C. ON FILE M3
      READ(M3) (VRES(I),I=1,NEQ)
C------  READ VECTOR FG MODIFIED BY B.C. AND MATRIX KG
      READ(M3) (VRES(I),I=1,NEQ)
      READ(M3) (VKGS(I),I=1,NKG),(VKGD(I),I=1,NEQ)
      IF(NSYM.EQ.1) READ(M3) (VKGI(I),I=1,NKG)
C------  EVALUATE THE RESIDUAL VECTOR
      DO 10 I=1,NEQ
10    VRES(I)=-VRES(I)
      CALL MULKU(VKGS,VKGD,VKGI,KLD,VFG,NEQ,NSYM,VRES)
      DO 20 I=1,NEQ
20    VRES(I)=-VRES(I)
      X1=ZERO
      DO 30 I=1,NEQ
      X=ABS(VRES(I))
      IF(X1.GE.X) GO TO 30
      X1=X
      I1=I
30    CONTINUE
      IF(M.GE.2) WRITE(MP,2000) (VRES(I),I=1,NEQ)
2000  FORMAT(/' RESIDUAL VECTOR'/(1X,10E12.5))
      WRITE(MP,2010) X1,I1
2010  FORMAT(/' MAX. RESIDUAL VALUE=',E12.5,' EQUATION',I5)
      RETURN 
      END


