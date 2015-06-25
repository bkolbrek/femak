      SUBROUTINE BLCOOR
C**********************************************************************C
C BLCOOR.FOR							       C
C  18/1-1992   Kompaniet Numerisk lyd				       C
C  Tonni F. Johansen						       C
C**********************************************************************C
C  Call and execute  Block 'COOR'				       C
C  Programmert etter Dhatt & Touzot, s.385			       C
C**********************************************************************C                                		     		 		                         
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 TBL

c      COMMON/COOR/NDIM,NNT,NDLN,NDLT,FAC(3),IAX
c      COMMON/ES/M,MR,MP,M1
c      COMMON /ALLOC/NVA
c      COMMON/LOC/LCORG,LDLNC
      COMMON/TRVL/FAC1(3),IN(3)
      COMMON VA(1)
      include 'femak.cmn'

      DIMENSION TBL(2)
      DATA ZERO/0.D0/,TBL/4HCORG,4HDLNC/
C------------------------------------------------------------------------
C--------- BLOCK HEADING
c      IF (M1.EQ.0) M1=MR
c      READ(M1,1000)IN,FAC1
      READ(MR,1000)IN,FAC1
1000  FORMAT(3I6,3F10.0)

C----------------DEFAULT OPTIONS
      IF (IN(1).GT.0) NNT=IN(1)
      IF (IN(2).GT.0) NDLN=IN(2)
      IF (IN(3).GT.0) NDIM=IN(3)

C----------------Check for axisymmetry
      IF (NDIM.EQ.20) THEN
	NDIM=2
	IAX=1
      ENDIF

      DO 10 I=1,3
        IF(FAC1(I).NE.ZERO) FAC(I)=FAC1(I)
10    CONTINUE

C--------------PRINT BLOCK PARAMETERS
      WRITE(MP,2000)M,NNT,NDLN,NDIM,FAC,NVA
2000  FORMAT(//'INPUT OF NODES(M=',I2,')'/' ',18('=')/
     1  15X,'MAX. NUMBER OF NODES             (NNT)=',I5/     
     2  15X,'MAX. NUMBER OF D.O.F PER NODE   (NDLN)=',I5/
     3  15X,'DIMENSIONS OF THE PROBLEM       (NDIM)=',I5/
     4  15X,'COORDINATE SCALE FACTORS         (FAC)=',3E12.5/
     5  15X,'WORKSPACE IN REAL WORDS          (NVA)=',I10)

C--------------ALLOCATE SPACE
      IF (LCORG.EQ.1) CALL ESPACE(NNT*NDIM,1,TBL(1),LCORG)
      IF (LDLNC.EQ.1) CALL ESPACE(NNT+1,0,TBL(2),LDLNC)

C---------------EXECUTE THE BLOCK
      CALL EXCOOR(VA(LCORG),VA(LDLNC))     


      RETURN
      END
