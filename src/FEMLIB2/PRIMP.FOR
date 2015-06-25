      SUBROUTINE PRIMP(PRIPAR,IMPPAR,VPRIG)
C*************************************************************
C    PRIMP.FOR
C    Ulf R. Kristiansen Tonni F.Johansen,  NTH
C    Kompaniet Numerisk Lyd 
C    12/8/92
C
C*************************************************************
C    TO PRINT DIFFERENT DATA OF IMPEDANCE GROUP
C     PRIPAR .EQ.0   Print a dashed line and headin for connection
C                    between element, side and group of impedances
C     PRIPAR .EQ.1   Print input cards defining the impedance group
C                    IMPPAR, and Impedance values
C     PRIPAR .EQ.2 , Print connection between element, side of
C                    element and group of imp.
C*************************************************************
C      IMPLICIT REAL*8(A-H,O-Z)

C----------------------------------------------------------------------

      INTEGER PRIPAR,IMPPAR
      COMPLEX*16 VPRIG(1)
C----------------------------------------------------------------------


      include 'femak.cmn'    

C----------------------------------------------------------------------

      INTEGER I2,I,J2,J
C----------------------------------------------------------------------
      IF (PRIPAR.EQ.0) THEN
C...........Print heading for connectivity
      WRITE(MP,1990)
1990  FORMAT(5x,72('-')/20X,' ELEMENT NO.',9X,' ELEMENT SIDE NO.',
     1  4X,' IMPEDANCE GROUP NO.')


      ELSEIF (PRIPAR.EQ.1) THEN
C..........Print properties of a group of impedances
      WRITE(MP,2000)IGIMP,IMPPAR
2000  FORMAT(/5X,82('-')/5X,' IMPEDANCE GROUP NO. :',I5,
     1  ' IS DEFINED USING READ PARAMETER  (IMPPAR)=',
     2  I5/5X,82('-'))

      IF (MOD(NFREQ,4).EQ.0) THEN
        I2=NFREQ/4
      ELSE
        I2=NFREQ/4+1
      ENDIF
      DO I=1,I2
        J2=MIN(4,NFREQ-(I-1)*4)
        WRITE(MP,2005)(I-1)*4+1,
     1    (VPRIG((IGIMP-1)*NFREQ+(I-1)*4+J),J=1,J2)
      ENDDO
2005  FORMAT(10X,' IMPEDANCE NO :',I5,2X,4('(',2E12.5,'),'))

      ELSEIF (PRIPAR.EQ.2) THEN

C.........Print connection betweem element, side and imp.grp.
      WRITE(MP,2020)IEL,IFL,IGIMP
2020  FORMAT(20X,I5,20X,I5,20X,I5)
      ENDIF



      RETURN
      END
