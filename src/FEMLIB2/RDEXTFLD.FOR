       SUBROUTINE RDEXTFLD(VCOEXT)
C***********************************************************************C
C     RDEXTFLD.FOR							C
C***********************************************************************C
C     1/2-1993								C
C     Kompaniet Numerisk Lyd						C
C     Tonni F. Johansen 						C
C***********************************************************************C
C     TO READ COORDINATES FOR THE EXTERNAL FIELD POINTS 		C
C     The same syntax to read coordinates as used in BLOCK COOR		C
C***********************************************************************C

      REAL*8 VCOEXT(1)

c      INTEGER NDIM,NNT,NDLN,NDLT
c      COMMON/COOR/NDIM,NNT,NDLN,NDLT

c      INTEGER M,MR,MP,M1,M2,M3
c      COMMON/ES/M,MR,MP,M1,M2,M3

      REAL*8 X1(3),X2(3)
      COMMON/TRVL/X1,X2

c      INTEGER NGIMP,NGSRC,NFREQ,NSIMP,NSRC,PRPAR,VELPAR,INTPAR,NPEXT
c      COMMON/AKUS/NGIMP,NGSRC,NFREQ,NSIMP,NSRC,PRPAR,VELPAR,INTPAR,NPEXT

      include 'femak.cmn'    


      INTEGER I,IN1,IN2,INCR,IDLN,I1,I2,I3,IN
      REAL*8 SPECL
      DATA SPECL/1.23456789D31/

C..........Initialize coordinates
      DO I=1,NPEXT*NDIM
	VCOEXT(I)=SPECL
      END DO

C----------READ EXTERN FIELDPOINTS DATA CARDS
      IF (M.GT.0) WRITE(MP,2000)
2000  FORMAT(//' EXTERNAL FIELDPOINTS DATA CARDS'/)
C........Read a DATA CARD
C      READ(M1,1000) IN1,X1,IN2,X2,INCR,IDLN
      READ(MR,1000) IN1,X1,IN2,X2,INCR,IDLN
1000  FORMAT(2(I5,3F10.0),2I5)
      IF (M.GT.0) WRITE(MP,2010) IN1,X1,IN2,X2,INCR,IDLN
2010  FORMAT('>>>>>',2(I5,3E12.5),2I5)

C.........IN1.EQ.0 terminate input of points
      DO WHILE (IN1.NE.0)

C-----------DECODE THE CARD
      IF (IN1.GT.NPEXT) CALL ERREUR(11,IN1,NPEXT,0)
      IF (IN2.GT.NPEXT) CALL ERREUR(12,IN2,NPEXT,0)
      IF (IN2.LE.0) IN2=IN1
      IF (IDLN.GT.NDLN) CALL ERREUR(13,IDLN,NDLN,0)
      IF (IDLN.LE.0) IDLN=NDLN
      IF (INCR.EQ.0) INCR=1

      I1=(IN2-IN1)/INCR
      I2=IN1+I1*INCR
      IF (I1.EQ.0) I1=1
      IF (IN2.NE.I2) CALL ERREUR(14,IN2,IN2,0)

C--------------GENERATE NODES BY INTERPOLATION
      DO I=1,NDIM
        X1(I)=X1(I)
        X2(I)=X2(I)
        X2(I)=(X2(I)-X1(I))/I1
      END DO

      I1=0
      I2=(IN1-1)*NDIM+1
      I3=(INCR-1)*NDIM

      DO IN=IN1,IN2,INCR
        IF (VCOEXT(I2).NE.SPECL) CALL ERREUR(15,IN,IN,0)
        DO I=1,NDIM
          VCOEXT(I2)=X1(I)+X2(I)*I1
          I2=I2+1
        END DO
        I1=I1+1
        I2=I2+I3
      END DO

C........Read a DATA CARD
C      READ(M1,1000) IN1,X1,IN2,X2,INCR,IDLN
      READ(MR,1000) IN1,X1,IN2,X2,INCR,IDLN
      IF (M.GT.0) WRITE(MP,2010) IN1,X1,IN2,X2,INCR,IDLN
      END DO

C--------------CHECK FOR MISSING NODES
      I1=NPEXT*NDIM+1
      I2=0
      I3=NPEXT+1
      DO I=1,NPEXT
        I1=I1-NDIM
        I3=I3-1

      IF (VCOEXT(I1).NE.SPECL) THEN
        IF (I2.EQ.0) I2=I3
      ELSEIF (VCOEXT(I1).EQ.SPECL) THEN
        IF (I2.EQ.0) CALL ERREUR(16,I3,I3,0)
        IF (I2.NE.0) CALL ERREUR(17,I3,I3,1)
      ENDIF
      END DO

      IF (I2.NE.NPEXT) CALL ERREUR(18,NPEXT,I2,0)


C-----------OUTPUT
      IF (M.GE.2) THEN
	WRITE(MP,2020)
2020	FORMAT(/10X,'NODE ',5X,'X',11X,'Y',11X,'Z'/)
	I1=1
	I2=NDIM

	DO IN=1,NPEXT
	  WRITE(MP,2030)IN,(VCOEXT(I),I=I1,I2)
2030	  FORMAT(10X,I5,3E12.5)
	  I1=I1+NDIM
	  I2=I2+NDIM
	END DO

      ENDIF


      RETURN
      END
