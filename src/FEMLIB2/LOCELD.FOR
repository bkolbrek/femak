      SUBROUTINE LOCELD(KDLNC,KNE,KNEQ,KLOCE,KLD)
C**********************************************************************C
C LOCELD.FOR							       C
C  14/2-1992	Kompaniet Numerisk lyd				       C
C  Tonni F. Johansen						       C
C**********************************************************************C
C  To form the element localization table (LOCE)		       C
C  and update coloumn heights fir a given element		       C
C  Programmert etter Dhatt & Touzot, s.397-398			       C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

c      COMMON/COOR/NDIM,NNT,NDLN
c      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
c     &		  IPG,ICODE,IDLE0,INEL0,IPG0

      include 'femak.cmn'    

      DIMENSION KDLNC(1),KNE(1),KNEQ(1),KLOCE(1),KLD(1)
      DATA NDLMAX/32000/

C----------------------------------------------------------------------
C---------------GENERATE KLOCE FROM KNEQ
      IDLE=0
      LOCMIN=NDLMAX
      DO 20 IN=1,INEL
	INN=KNE(IN)
	IF (INN.GT.NNT) CALL ERREUR(56,INN,NNT,1)
	IEQ=KDLNC(INN)
	IEQ1=KDLNC(INN+1)
10	IF (IEQ.GE.IEQ1) GO TO 20
	IEQ=IEQ+1
	IDLE=IDLE+1
	J=KNEQ(IEQ)
	KLOCE(IDLE)=J
	IF (J.LT.LOCMIN.AND.J.GT.0) LOCMIN=J
	GO TO 10
20    CONTINUE

C-------------UPDATE TABLE OF COLOUMN HEIGHTS (KLD)
      DO 30 ID=1,IDLE
	J=KLOCE(ID)
	IF (J.LE.0) GO TO 30
	IH=J-LOCMIN
	IF (IH.GT.KLD(J+1)) KLD(J+1)=IH
30    CONTINUE

      RETURN
      END
