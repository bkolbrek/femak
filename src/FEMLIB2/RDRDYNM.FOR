      SUBROUTINE RDRDYNM(KDLNC,VCORG,VDIMP,KNEQ,VSRC,VFREQ)
C**************************************************************
C     PRDYNM.FOR
C**************************************************************
C     4/9/1992
C     Kompaniet Numerisk Lyd
C     Tonni F. Johansen
C**************************************************************
C     TO PRINT SOLUTION for dynamic problem
C==============================================================


      include 'femak.cmn'

      REAL*8 VA(1)
      COMMON VA

      REAL*8 VDIMP(1),VCORG(1),VFREQ(1)
      COMPLEX*16 VSRC(1)
      INTEGER KDLNC(1),KNEQ(1)

      INTEGER I
      
      REWIND(M10)
      READ(M10)VFREQ(IFREQ)
      READ(m10)(VSRC(I),I=1,NEQ)


      RETURN
      END
