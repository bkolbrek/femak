      SUBROUTINE BLSRCE
C****************************************************************
C     BLSRCE.FOR
C     Tonni F. Johansen, NTH
C     24/8-92
C***************************************************************
C     TO CALL BLOCK 'SRCE'
C     TO READ SRCE SPECIFICATIONS
C*************************************************************** 
c      IMPLICIT REAL*8(A-H,O-Z)

      include 'femak.cmn'

      REAL*8 VA(1)
      COMMON VA

      INTEGER NFREQ1

      REAL*4 TBL(3)
      DATA TBL/4HPRSG,4HTPRS,4HNDSC/

C--------------------------------------------------------------
C      IF(M1.EQ.0) M1=MR
C      READ(M1,1000)NGSRC,NFREQ1,NSRC
      READ(MR,1000)NGSRC,NFREQ1,NSRC
1000  FORMAT(3I5)
      IF (NFREQ.EQ.0) THEN
        NFREQ=NFREQ1
      ELSEIF (NFREQ.NE.NFREQ1) THEN
        WRITE(*,*)'FEIL ANTALL FREKVENSER'
      ENDIF

      WRITE(MP,2000)M,NGSRC,NFREQ,NSRC
2000  FORMAT(//' INPUT OF SOURCE SPECIFICATIONS (M=',I2,')'/' ',36('=')/
     115X,' NUMBER OF SOURCE GROUPS              (NGSRC) =',I2,/
     215X,' NUMBER OF FREQUENCIES                (NFREQ) =',I2,/
     315X,' TOTAL NUMBER OF VELOCITY SURFACES    (NSRC)  =',I2)

C------------Allocate space for global source arrays (KTYPRSG, VPRSG)
      IF (LPRSG.EQ.1) CALL ESPACE(NGSRC*NFREQ,2,TBL(1),LPRSG)
      IF (LTYPRSG.EQ.1) CALL ESPACE(NGSRC,0,TBL(2),LTYPRSG)

C----------Allocate space for local source array (VNDSRC)
      IF (LNDSRC.EQ.1) THEN
        IF (NDIM.EQ.2) THEN
C-------------for 2 dimentions: 3 nodes per side
c          CALL ESPACE(NGSRC*NFREQ*4,2,TBL(3),LNDSRC)
          CALL ESPACE(NFREQ*4,2,TBL(3),LNDSRC)
        ELSEIF (NDIM.EQ.3) THEN
C-------------for 3 dimentions: 8 nodes per side
c          CALL ESPACE(NGSRC*NFREQ*10,2,TBL(3),LNDSRC)
          CALL ESPACE(NFREQ*10,2,TBL(3),LNDSRC)
        ENDIF
      ENDIF

      CALL EXSRCE(VA(LPRSG),VA(LTYPRSG),VA(LNDSRC),VA(LLOCE),VA(LCORE),
     1    VA(LNE),VA(LPRNE),VA(LPREE),VA(LPRNG))
  
      RETURN
      END 
