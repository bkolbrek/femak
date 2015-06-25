      SUBROUTINE ASGRADC(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VSRE, 
     1  VSRC,VDLEC,VRES)
C********************************************************************** 
C     ASGRADC.FOR
C**********************************************************************
C     11/11 1994
C     NTH
C     Tonni F. Johansen
C     Programmert etter Dhatt og Touzot s.410
C**********************************************************************  
C     TO EVALUATE AND PRINT GRADIENTS ( STRESSES ) AT ELEMENT G.P.
C       ( ELEMENT FUNCTION 8)                  
CNOTE!!!!
C!!!!    
C  Test version only for ELEm03, will not work with any others.
C======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      include 'femak.cmn'

      COMPLEX*16 VIE(1),VSRE(1),VSRC(1),VDLEC(1)
      REAL*8 VDIMP(1),VCORE(1),VPRNE(1),VPREE(1),VRES(1)
      INTEGER KLD(1),KLOCE(1),KNE(1)
      
C----------------------------------------------------------------------
C---------------  REWIND ELEMENTS FILE (M1)
      REWIND M1
C--------------   LOOP OVER THE ELEMENTS
      DO 20 IE=1,NELT
C-------------    READ THE ELEMENT
      CALL RDELEM(M2,KLOCE,VCORE,VPRNE,VPREE,KNE)
C-------------    EVALUATE INTERPOLATION FUNCTION IF REQUIRED
      IF(ITPE.EQ.ITPE1) GO TO 10
      ICOD=2
      CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VKE,VKE1,VIE,VIE1,VFE,VSRE,
     1   VNDSRC,KTYPRSG,VPRSG)
C------------     FIND ELEMENT D.O.F.
10    CALL DLELMC(KLOCE,VSRC,VDIMP,VDLEC)
C------------     COMPUTE AND PRINT STRESSES OR GRADIENTS
      ICOD=8
      CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VKE,VKE1,VIE,VIE1,VFE,VSRE,
     1   VNDSRC,KTYPRSG,VPRSG)
20    ITPE1=ITPE
      RETURN
      END

