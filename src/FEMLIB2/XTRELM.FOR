      SUBROUTINE XTRELM(IGPE,VCORG,VPRNG,VPREG,KNE,VCORE,VPRNE,VPREE)
C**********************************************************************C
C XTRELM.FOR							       C
C  14/2-1992	Kompaniet Numerisk lyd				       C
C  Tonni F. Johansen						       C
C**********************************************************************C
C  To generate element coordinates and properties from		       C
C  global arrays						       C
C    (IGPE:  group number for element properties)		       C
C  Programmert etter Dhatt & Touzot, s.398-399			       C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

c      COMMON/COOR/NDIM,NNT,NDLN
c      COMMON/PRND/NPRN
c      COMMON/PREL/NGPE,NPRE
c      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
c     &		  IPG,ICODE,IDLE0,INEL0,IPG0

      include 'femak.cmn'    

      DIMENSION VCORG(1),VPRNG(1),VPREG(1),KNE(1),VCORE(1),
     1	     VPRNE(1),VPREE(1)

C----------------------------------------------------------------------
C---------------GENERATE ELEMENT COORDINATES
      IPRNE=0
      ICE=0
      DO 30 IN=1,INEL
	IC=(KNE(IN)-1)*NDIM
	DO 10 I=1,NDIM
	  ICE=ICE+1
	  IC=IC+1
	  VCORE(ICE)=VCORG(IC)
10	CONTINUE

C--------------GENERATE ELEMENT NODAL PROPERTIES
	IF (NPRN.EQ.0) GO TO 30
	IC=(KNE(IN)-1)*NPRN
	DO 20 I=1,NPRN
	  IPRNE=IPRNE+1
	  IC=IC+1
	  VPRNE(IPRNE)=VPRNG(IC)
20	CONTINUE
30    CONTINUE

C---------GENERATE ELEMENT PROPERTIES
      IPREE=0
      IF (NPRE.EQ.0) GO TO 50
      IC=(IGPE-1)*NPRE
      DO 40 I=1,NPRE
	IPREE=IPREE+1
	IC=IC+1
	VPREE(IPREE)=VPREG(IC)
40    CONTINUE

50    RETURN
      END
