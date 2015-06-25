      SUBROUTINE ESPACE(ILONG,IREEL,TBL,IDEB)
C****************************************************************************C
C   espace.for                                                               c
C   31/1-1992                                                                C
C   Kompaniet Numerisk lyd, Tonni F. Johansen                                C
C****************************************************************************C
C  To allocate a real or integer table in array VA                           C
C    Input                                                                   C
C       ILONG           Length of the table to be allocated                  C
C                       (in real or integer words)                           C
C       IREEL           Table type                                           C
C                          .EQ.0    Integer                                  C
C                          .EQ.1    Real                                     C
C                          .EQ.2    Complex                                  C
C       TBL             Name of the table (A4)                               C
C    Output                                                                  C
C       IDEB            Table to be allocated starts in VA(IDEB)             C
C                                                                            C
C        Programmed after Dhatt & Touzot p.368                               C
C****************************************************************************C

      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 TBL
      
c      COMMON/ES/M,MR,MP
c      COMMON/ALLOC/NVA,IVA,IVAMAX,NREEL

      include 'femak.cmn'

      COMMON VA(1)

      DIMENSION KA(1)
      EQUIVALENCE (VA(1),KA(1))
      DATA ZERO/0.D0/

C---------------------------------------------------------------------------
C----------------CALCULATE THE TABLE LENGTH IN REAL WORDS
	 
      IF (IREEL.EQ.0) THEN
C............Integer       
	ILGR=(ILONG+NREEL-1)/NREEL

      ELSEIF (IREEL.EQ.1) THEN
C............Real       
	ILGR=ILONG

      ELSEIF (IREEL.EQ.2) THEN
C............Complex       
	ILGR=ILONG*2
      ENDIF

      IVA1=IVA+ILGR

C------------------CHECK IF ENOUGH SPACE IS AVAILABLE`
      IF (IVA1.LE.NVA) GO TO 20

C............AUTOMATIC EXTENSION OF THE BLANK COMMON IF CORRESPONDING
C            SYSTEM COMMAND EXIST ON THE COMPUTER USED
C      CALL EXTEND(IVA1,IERR)
C      IF (IERR.EQ.1) GO TO 10
C      NVA=IVA1
C      GO TO 20

C------------ALLOCATION ERROR (NOT ENOUGH SPACE)
10    WRITE(MP,2000)TBL,IVA1,NVA
2000  FORMAT('**** ALLOCATION ERROR, TABLE ',A4/'REQUIRED SPACE:',I9,
     1  'REAL WORDS, AVAILABLE SPACE:',I9,'REAL WORDS')

      STOP
       
C-----------ALLOCATE TABLE
20    IDEB=IVA+1
      IVA=IVA1
      IF (IVA.GT.IVAMAX) IVAMAX=IVA
      IF (M.GT.0) WRITE(MP,2010)TBL,IDEB,IVA1
2010  FORMAT(60X,'TABLE ',A4,' GOES FROM VA(',I7,') TO VA(',I7,')')

C----------INITIALIZE THE ALLOCATED TABLE TO ZERO
      I1=IDEB
      IF (IREEL.EQ.0) THEN
	I1=(I1-1)*NREEL+1
	I2=I1+ILONG-1
	DO 30 I=I1,I2
	  KA(I)=0
30    CONTINUE

      ELSEIF (IREEL.EQ.1) THEN
       I2=I1+ILONG-1
	DO 40 I=I1,I2
	  VA(I)=ZERO
40    CONTINUE


       ELSEIF (IREEL.EQ.2) THEN
       I2=I1+2*ILONG-1
	DO 50 I=I1,I2
	  VA(I)=ZERO
50    CONTINUE

      ENDIF

      RETURN
      END
