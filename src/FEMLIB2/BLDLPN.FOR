      SUBROUTINE BLDLPN
C**********************************************************************C
C BLDLPN.FOR                                                           C
C  28/9-1994   Kompaniet Numerisk lyd                                  C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  Call   Block 'DLPN'                                                 C
C  Read number of DOF per node                                         C
C  Programed after Dhatt & Touzot, p.384                               C
C**********************************************************************C
      IMPLICIT REAL*8 (V-W)

      include 'femak.cmn'
      
      COMMON VA(1)
      
c      IF (M1.EQ.0) M1=MR
      WRITE(MP,2000)M
2000  FORMAT(//' INPUT OF D.O.F. (M=',I2,')'/' ',17('='))
      CALL EXDLPN(VA(LDLNC))      

      RETURN
      END
