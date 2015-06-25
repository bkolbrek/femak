      SUBROUTINE CALINWV(KLOCE,VCORE,VPRNE,VPREE,KNE,VNDSRC,
     1       KTYPRSG,VPRSG,VPRNG,VCORG,NODNR,PRES)
C*******************************************************************C
C     CALINWV.FOR                                                   C
C*******************************************************************C
C     21/4-1995                                                     C
C     NTH                                                           C
C     Tonni F. Johansen                                             C
C                                                                   C
C*******************************************************************C
C     TO calculate contribution to total field on ONE node from a   C
C     an incoming wave                                              C
C*******************************************************************C
C  Right now, only pressure                                         C
C*******************************************************************C
C
C

      COMPLEX*16 VNDSRC(1),VPRSG(1),PRES
      REAL*8 VCORE(1),VPRNE(1),VPREE(1)
      REAL*8 VPRNG(1),VCORG(1)
      INTEGER KLOCE(1),KNE(1),KTYPRSG(1),NODNR

C
      include 'femak.cmn'

c   $DECLARE                             
      INTEGER IE,GRPNO

      COMPLEX*16 INPRES
C----------------------------------------------------------------------
C------ REWIND SOURCE INFO FILE (M7), and reset element type
      REWIND M7
      ITPE1=0
C------LOOP OVER THE ELEMENTS      
      DO IE=1,NSRC
C------ READ AN SOURCE INFO ON FILE M7
        CALL RDSRCE(M7,KLOCE,VCORE,VPRNE,VPREE,KNE,VNDSRC,KTYPRSG)
        IF (KTYPRSG(IGSRC).EQ.50) THEN
          GRPNO=NINT(REAL(VPRSG((IGSRC-1)*NFREQ+IFREQ)))
          IF(NINT(VPRNG((NODNR-1)*NPRN+1)).EQ.GRPNO) THEN
            CALL CALTYP50(VPREE,VNDSRC,VCORG,NODNR,INPRES)
            PRES=PRES+INPRES
          ENDIF
C..............End IF (KTYPRSG.EQ.50)
        ELSEIF (KTYPRSG(IGSRC).EQ.55) THEN
          GRPNO=NINT(REAL(VPRSG((IGSRC-1)*NFREQ+IFREQ)))
          IF(NINT(VPRNG((NODNR-1)*NPRN+1)).EQ.GRPNO) THEN
C...........INTERFACE node
            CALL CALTYP50(VPREE,VNDSRC,VCORG,NODNR,INPRES)
            PRES=PRES+INPRES
          ELSEIF(NINT(REAL(VPRNG((NODNR-1)*NPRN+1))).EQ.GRPNO) THEN
C....NB!!!IMAG IS NOT STANDARDFORTRAN...
C.........  Node in WEEM
            CALL CALTYP50(VPREE,VNDSRC,VCORG,NODNR,INPRES)
            PRES=PRES+INPRES
          ENDIF
C..............End IF (KTYPRSG(IGSRC).EQ.55)        
        ENDIF       
      ENDDO       

      RETURN
      END

                               