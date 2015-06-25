      SUBROUTINE EXCIMP(VPRIG,KTYPRIG,KLOCE,VCORE,KNE,VPRNE,VPREE,
     1                  SPCELM)
C************************************************************
C     EXCIMP.FOR
C     Ulf R Kristiansen & Tonni F. Johansen
C     Kompaniet Numerisk Lyd
C     15/10-92
c
C************************************************************
C     TO EXECUTE BLOCK 'CIMP'
C     TO READ IMPEDANCE PROPERTIES
C************************************************************
      IMPLICIT REAL*8(A-H,O-Z)

      include 'femak.cmn'

      COMPLEX*16 VPRIG(1)
      REAL*8 VCORE(1),VPRNE(1),VPREE(1)
      INTEGER KLOCE(1),KNE(1),KTYPRIG(1)
      INTEGER SPCELM

      INTEGER IMPPAR,I,J,II,JJ,IK,I2,J2,IEL0

      II=0
      JJ=0

C-----------------------------------------------------------
C------- READ IMPEDANCE GROUP, IMPEDANCE PARAMETER AND
C                                IMPEDANCES
      DO IK=1,NGIMP
      READ(MR,1000)IGIMP,IMPPAR
1000    FORMAT(2I5)
      KTYPRIG(IGIMP)=IMPPAR
       
      IF (IMPPAR.EQ.0) THEN
        IF (MOD(NFREQ,4).EQ.0) THEN
          I2=NFREQ/4
        ELSE
          I2=NFREQ/4+1
        ENDIF
        DO 10 I=1,I2
          J2=MIN(4,NFREQ-(I-1)*4)
          READ(MR,1010)(VPRIG((IGIMP-1)*NFREQ+(I-1)*4+J),J=1,J2)
10      CONTINUE
        IF(M.GT.0) THEN
          CALL PRIMP(1,IMPPAR,VPRIG)
        ENDIF
1010      FORMAT(8E12.5)

      ELSEIF (IMPPAR.EQ.100) THEN
C............KAGAWA SPECIAL CIRCULAR HOLE IN INFINITE BAFFLE
        READ(MR,1020)NI100MOD,NI100NOD,I100RAD
1020      FORMAT(2I5,E12.5)
        NI100FL=(NI100NOD-1)/2
        IF(M.GT.0) THEN
          CALL PRIMP(1,IMPPAR,VPRIG)
          WRITE(MP,1022)NI100MOD,NI100NOD,I100RAD
1022        FORMAT('MAX.NUMBER OF MODES            (NI100MOD) :',I5/,
     1             'NUMBER OF NODES                (NI100NOD) :',I5/,
     2             'RADIUS OF HOLE                 (I100RAD)  :',E12.5)
        ENDIF


C............Set special impedance flag (SPCELM)
        SPCELM=1
      ENDIF

C!!!Husk at vi maa ha feilmeldinger !!
      ENDDO

C........Write text for next input
      IF(M.GT.0) THEN
      CALL PRIMP(0,0,VPRIG)
      ENDIF


C........Read info. about elem. with impedance. Put elem.prop. on file
C........NSIMP Number of Surfaces with IMPedance
      DO IK=1,NSIMP
        READ(MR,1150)IEL0,IFL,IGIMP
1150    FORMAT(3I5)
C.........No action if group has imppar = 100
C       IF (KTYPRIG(IGIMP).LT.100) THEN
C         CALL GETI100(VPRIG,KLOCE,VCORE,KNE,VPRNE,VPREE,IEL0)

C       ELSE
C...........Search for information about IEL in elementfile
          REWIND(M1)
          IEL=0
C...........Search "real-elements"
          NSEL=0
          DO WHILE ((IEL0.NE.IEL).AND.(NSEL.LT.NELTR))
            CALL RDELEM(M1,KLOCE,VCORE,VPRNE,VPREE,KNE)
            NSEL=NSEL+1
          ENDDO     

C............If not among "real-elements" searc "complex-elements"
          IF (IEL.NE.IEL0) THEN
            REWIND(M2)
            IEL=0
C...........Serch "complex-elements"
            NSEL=0
            DO WHILE ((IEL0.NE.IEL).AND.(NSEL.LT.NELTC))
              CALL RDELEM(M2,KLOCE,VCORE,VPRNE,VPREE,KNE)
              NSEL=NSEL+1
            ENDDO     
          
          ENDIF          
C...........If found, write element properties to imp.surf.file
          IF (IEL.EQ.IEL0) THEN
            CALL WRCIMP(M6,KLOCE,VCORE,VPRNE,VPREE,KNE,IGIMP)
            IF (M.GT.0) THEN
              CALL PRIMP(2,0,VPRIG)
            ENDIF
          ELSE
            WRITE(*,1160)IEL0
1160          FORMAT('ERROR; IEL0 ,'I4,' DOES NOT EXIST !')
          ENDIF


      ENDDO


      RETURN
      END
