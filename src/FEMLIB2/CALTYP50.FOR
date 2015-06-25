      SUBROUTINE CALTYP50(VPREE,VNDSRC,VCORG,NODNR,INPRES)
C*******************************************************************C
C     CALTYP50.FOR                                                  C
C*******************************************************************C
C     21/4-1995                                                     C
C     NTH                                                           C
C     Tonni F. Johansen                                             C
C                                                                   C
C*******************************************************************C
C     To find incoming pressure in a specified coordinate           C
C     coordinate specified in VCORG (may be changed ?)              C
C*******************************************************************C
C
C
      COMPLEX*16 VNDSRC(1),INPRES
      REAL*8 VPREE(1),VCORG(1)
      INTEGER NODNR
      include 'femak.cmn'

C $declare

      REAL*8 KX,KY,A
      

      A=VNDSRC((IFREQ-1)*4+1)
      KX=REAL(VNDSRC((IFREQ-1)*4+2))
      KY=REAL(VNDSRC((IFREQ-1)*4+3))
      INPRES=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORG(NDIM*NODNR-1)*KX+
     1         VCORG(NDIM*NODNR)*KY))
          


      RETURN
      END     
      