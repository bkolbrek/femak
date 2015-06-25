      SUBROUTINE BLCOND
C**********************************************************************C
C BLCOND.FOR							       C
C  1/2-1992   Kompaniet Numerisk lyd				       C
C  Tonni F. Johansen						       C
C**********************************************************************C
C  Call and execute  Block 'COND'				       C
C  Programmed after Dhatt & Touzot, p.389			       C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 TBL

c      COMMON/COOR/NDIM,NNT,NDLN,NDLT,FAC(3)
c      COMMON/COND/NCLT,NCLZ,NCLNZ
c      COMMON /ALLOC/NVA,IVA
c      COMMON/ES/M,MR,MP,M1
c      COMMON/LOC/LCORG,LDLNC,LNEQ,LDIMP

      include 'femak.cmn'

      COMMON VA(1)

      DIMENSION TBL(2)
      DATA TBL/4HNEQ ,4HDIMP/

C------------------------------------------------------------------------
c      IF (M1.EQ.0) M1=MR
      WRITE(MP,2000)M
2000  FORMAT(//' INPUT OF BOUNDARY CONDITIONS (M=',I2,')'/' ',
     1	 33('=')/)
      IF (LNEQ.EQ.1) CALL ESPACE(NDLT,0,TBL(1),LNEQ)
      IF (LDIMP.EQ.1) CALL ESPACE(NDLT,1,TBL(2),LDIMP)
      CALL EXCOND(VA(LCORG),VA(LDLNC),VA(LNEQ),VA(LDIMP))
      CALL VIDE(LDIMP+NCLT,1,TBL(2))


      RETURN
      END
