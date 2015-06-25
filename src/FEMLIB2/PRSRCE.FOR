      SUBROUTINE PRSRCE(PRIPAR,SRCPAR,VPRSG,VNDSRC)
C*************************************************************
C    PRSRCE.FOR
C    Tonni F. Johansen,  NTH
C    9/9/92
C*************************************************************
C    TO PRINT DIFFERENT DATA OF IMPEDANCE GROUP
C     PRIPAR .EQ.0   Print a dashed line and headin for connection
C                    between element, side and group of impedances
C     PRIPAR .EQ.1   Print input cards defining the source group
C                    SRCPAR, and source values
C     PRIPAR .EQ.2 , Print connection between element, side of
C                    element and group of source
C*************************************************************

C----------------------------------------------------------------------

      INTEGER PRIPAR,SRCPAR
      COMPLEX*16 VPRSG(1),VNDSRC(1)
C----------------------------------------------------------------------

      include 'femak.cmn'

C----------------------------------------------------------------------

      INTEGER I2,I,J2,J
C----------------------------------------------------------------------
      IF (PRIPAR.EQ.0) THEN
C...........Print heading for connectivity
      WRITE(MP,1990)
1990  FORMAT(5x,72('-')/20X,' ELEMENT NO.',9X,' ELEMENT SIDE NO.',
     1    4X,' SOURCE GROUP NO.')


      ELSEIF (PRIPAR.EQ.1) THEN
C..........Print properties of a group of impedances
      WRITE(MP,2000)IGSRC,SRCPAR
2000  FORMAT(/5X,82('-')/5X,' SOURCE GROUP NO. :',I5,
     1  ' IS DEFINED USING READ PARAMETER  (SRCPAR)=',
     2  I5/5X,82('-'))

      IF (SRCPAR.EQ.0) THEN
        IF (MOD(NFREQ,4).EQ.0) THEN
          I2=NFREQ/4
        ELSE
          I2=NFREQ/4+1
        ENDIF
        DO I=1,I2
          J2=MIN(4,NFREQ-(I-1)*4)
          WRITE(MP,2005)(I-1)*4+1,
     1    (VPRSG((IGIMP-1)*NFREQ+(I-1)*4+J),J=1,J2)
        ENDDO
2005    FORMAT(10X,' SOURCE NO :',I5,2X,4('(',2E12.5,'),'))
      ENDIF

      ELSEIF (PRIPAR.EQ.2) THEN

C.........Print connection betweem element, side and imp.grp.
      WRITE(MP,2020)IEL,IFL,IGIMP
2020  FORMAT(20X,I5,20X,I5,20X,I5)
      IF (SRCPAR.EQ.20) THEN
        IF (NDIM.EQ.2) THEN
          DO I=1,NFREQ
            IF (I.EQ.1) THEN
              WRITE(MP,2025)I,(VNDSRC((I-1)*4+J),J=1,4)
            ELSE
              WRITE(MP,2026)I,(VNDSRC((I-1)*4+J),J=1,4)
            ENDIF
          ENDDO
2025      FORMAT(36X,' VELOCITY AT FREQUENCY :',I5,3('(',2E12.5,')'))
2026      FORMAT(60X,I5,3('(',2E12.5,')'))
        ENDIF
      ENDIF
      ENDIF



      RETURN
      END
