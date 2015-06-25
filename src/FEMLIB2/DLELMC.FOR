      SUBROUTINE DLELMC(KLOCE,VDLG,VDIMP,VDLEC)
C****************************************************************
C     DLELMC.FOR
C****************************************************************
C     25/2/1992
C     Kompaniet Numerisk Lyd
C     Ulf R. Kristiansen
C     Programmert etter Dhatt og Touzot s.416
C****************************************************************
C     TO GENERATE ELEMENT D.O.F.
C================================================================
      IMPLICIT REAL*8(A-H,O-Z)
c      COMMON/RGDT/IEL,INUL(3),IDLE
c      COMMON/ES/M,MR,MP

      include 'femak.cmn'

      INTEGER KLOCE(1)
      REAL*8 VDIMP(1)
      COMPLEX*16 VDLG(1),VDLEC(1)
      DATA ZERO/0.D0/    
      
C----------------------------------------------------------------
      DO 40 ID=1,IDLE
      IL=KLOCE(ID)
      IF(IL) 10,20,30
10    VDLEC(ID)=VDIMP(-IL)
      GO TO 40
20    VDLEC(ID)=ZERO
      GO TO 40
30    VDLEC(ID)=VDLG(IL)
40    CONTINUE
      IF(M.GE.2) WRITE(MP,2000) IEL,(VDLEC(ID),ID=1,IDLE)
2000  FORMAT(' DEGREES OF FREEDOM OF ELEMENT ',I5/(1X,10E12.5))
      RETURN
      END
