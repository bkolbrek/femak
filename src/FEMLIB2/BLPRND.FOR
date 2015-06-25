      SUBROUTINE BLPRND
C**********************************************************************C
C BLPRND.FOR                                                           C
C  21/3-1995  NTH                                                      C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  Call and execute  Block 'PRND'                                      C
C  To read nodal properties                                            C
C  Programmed after Dhatt & Touzot, p.392                              C
C**********************************************************************C

      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'femak.cmn'
c $declare
      REAL*4 TBL
      
      COMMON VA(1)
      
      DATA TBL/4HPRNG/
      
      READ(MR,1000)NPRN
1000  FORMAT(I5)      
          
      
      WRITE(MP,2000) M,NPRN
2000  FORMAT(//' INPUT OF NODAL PROPERTIES (M=',I2,')'/' ',30('=')/
     1     15X,'NUMBER OF PROPERTIES PER NODE   (NPRN)=',I5)
     
      IF (LPRNG.EQ.1) CALL ESPACE(NNT*NPRN,1,TBL,LPRNG)
      CALL EXPRND(VA(LPRNG))    
          
      RETURN
      END
