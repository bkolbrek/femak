      SUBROUTINE PRDYNM(KDLNC,VCORG,VDIMP,KNEQ,VSRC,VFREQ,KLOCE,VCORE,
     1     VPRNE,VPREE,KNE,VNDSRC,KTYPRSG,VPRSG,VPRNG)
C**************************************************************
C     PRDYNM.FOR
C**************************************************************
C     4/9/1992
C     Tonni F. Johansen,  NTH
C**************************************************************
C     TO PRINT SOLUTION for dynamic problem
C==============================================================


      INTEGER GRPMAX
      PARAMETER (GRPMAX=10)

      include 'femak.cmn'
c $DECLARE
      REAL*8 VA(1)
      COMMON VA

      REAL*8 VDIMP(1),VCORG(1),VFREQ(1),VCORE(1),VPRNE(1),VPREE(1)
      REAL*8 VPRNG(1)
      COMPLEX*16 VSRC(1),VNDSRC(1),VPRSG(1)
      INTEGER KDLNC(1),KNEQ(1),KLOCE(1),KNE(1),KTYPRSG(1)

      COMPLEX*16 V(6),CZERO,PRES
      REAL*8 ZERO,X1,X2,X3,AMPL,DBAMPL,PHASE

      INTEGER I,I1,I2,ID,ID1,ID2,J,JJ,IN,NINWVE,NGRPS(GRPMAX),IE,GRPNO
      INTEGER NODNR,II
      INTEGER IORES
      
      LOGICAL INGRP
      REAL*4 RF,STAR,EMPTY
      DATA CZERO/(0.D0,0.D0)/,ZERO/0.D0/
      DATA RF/4H *  /,EMPTY/4H    /
C--------------------------------------------------------------
      X2=ZERO
      X3=ZERO
      
      WRITE(M10)VFREQ(IFREQ)
      DO IE=1,NEQ
        WRITE(m10)VSRC(IE)
      ENDDO
C      WRITE(m10,IOSTAT=IORES)(VSRC(I),I=1,NEQ)
C      IF (IORES.NE.0) THEN
C        WRITE(*,*)'ERROR WHILE WRITING TO M10,  ERROR NO. ', IORES
C      ENDIF

      WRITE(MP,1990)VFREQ(IFREQ)
1990  FORMAT(' COMPUTATIONS AT FREQUENCY:',F12.3,'Hz'/72('='))

c      WRITE(MP,2000)
c2000  FORMAT(/3X' NODES',4X,'X',11X,'Y',11X,'Z',18X,'COMPLEX',
c     1	     12X,'AMPL',8X,'dB',8X,'PHASE(RAD)',2X,'*=PRESCRIBED DOF'/)


C........Check first whether it is any incoming wave
      IF (INWVE) THEN
C.......Find number of incoming waves, and which nodegroups  are affected.

C.........Reset variables, arrays and files.
        DO I=1,GRPMAX
          NGRPS(I)=0
        ENDDO
        REWIND M7
        ITPE1=0
        NINWVE=0

C------LOOP OVER THE SOURCES      
        DO IE=1,NSRC
C------ READ AN SOURCE INFO ON FILE M7
          CALL RDSRCE(M7,KLOCE,VCORE,VPRNE,VPREE,KNE,VNDSRC,KTYPRSG)
C..........PLANE WAVES
          IF (KTYPRSG(IGSRC).EQ.50) THEN
            GRPNO=NINT(DBLE(VPRSG((IGSRC-1)*NFREQ+IFREQ)))
            NINWVE=NINWVE+1                
            INGRP=.FALSE.
            I=1
            DO WHILE ((I.LE.NINWVE).AND.(.NOT.INGRP))
              IF(NGRPS(I).EQ.GRPNO) INGRP=.TRUE.
              I=I+1                
            ENDDO
            IF (.NOT.INGRP) THEN
              NGRPS(I-1)=GRPNO
            ENDIF
          ELSEIF (KTYPRSG(IGSRC).EQ.55) THEN
c............source type 55
c............first nodes of interface FEM-WEEM
            GRPNO=NINT(DBLE(VPRSG((IGSRC-1)*NFREQ+IFREQ)))
            NINWVE=NINWVE+1                
            INGRP=.FALSE.
            I=1
            DO WHILE ((I.LE.NINWVE).AND.(.NOT.INGRP))
              IF(NGRPS(I).EQ.GRPNO) INGRP=.TRUE.
              I=I+1                
            ENDDO
            IF (.NOT.INGRP) THEN
              NGRPS(I-1)=GRPNO
            ENDIF
c............second nodes in WEEM     , 
c...nb!!!IMAG IS NOT STANDARD FORTRAN, IS F90??
            GRPNO=NINT(AIMAG(VPRSG((IGSRC-1)*NFREQ+IFREQ)))
            NINWVE=NINWVE+1                
            INGRP=.FALSE.
            I=1
            DO WHILE ((I.LE.NINWVE).AND.(.NOT.INGRP))
              IF(NGRPS(I).EQ.GRPNO) INGRP=.TRUE.
              I=I+1                
            ENDDO
            IF (.NOT.INGRP) THEN
              NGRPS(I-1)=GRPNO
            ENDIF

          ENDIF          
        ENDDO       
C!!Note!!
C.....Must have error message for too many groups.
      ENDIF

      I2=0
      DO IN=1,NNT
        I1=I2+1
        I2=I2+NDIM
        ID1=KDLNC(IN)+1
        ID2=KDLNC(IN+1)
        ID=ID2-ID1+1
        IF (ID2.GE.ID1) THEN
          X1=VCORG(I1)
          IF (NDIM.GE.2) X2=VCORG(I1+1)
          IF (NDIM.GE.3) X3=VCORG(I1+2)
          J=ID1
          DO I=1,ID
            JJ=KNEQ(J)
            IF(JJ.LT.0) THEN
              V(I)=VDIMP(-JJ)
              STAR=RF
            ELSEIF(JJ.EQ.0) THEN
              V(I)=CZERO
            ELSEIF(JJ.GT.0) THEN
              V(I)=VSRC(JJ)
              STAR=EMPTY
            ENDIF

C.............Add contributions from incoming waves.
            IF (INWVE) THEN
C......check whether node is influenced
              INGRP=.FALSE.
              JJ=1
              DO WHILE ((JJ.LE.NINWVE).AND.(.NOT.INGRP))
                IF(NGRPS(JJ).EQ.NINT(VPRNG((IN-1)*NPRN+1))) INGRP=.TRUE.
C??                IF(NGRPS(JJ).EQ.NINT(IMAG(VPRNG((IN-1)*NPRN+1))))
C     1                                                  INGRP=.TRUE.
                JJ=JJ+1                
              ENDDO

C.......if influenced, add contributions from incoming wave
              IF (INGRP) THEN
                CALL CALINWV(KLOCE,VCORE,VPRNE,VPREE,KNE,VNDSRC,
     1            KTYPRSG,VPRSG,VPRNG,VCORG,IN,V(1))                               
              ENDIF
            ENDIF

C            AMPL=SQRT(REAL(V(I))*REAL(V(I))+IMAG(V(I))*IMAG(V(I)))
C            IF (AMPL.GT.EPS) THEN
C              DBAMPL=20*LOG10(AMPL/1.0)
C              PHASE=ATAN2(IMAG(V(I)),REAL(V(I)))             
C            ENDIF
            J=J+1
          ENDDO
c        WRITE(m10)(V(II),II=1,ID)
          WRITE(MP,2010)IN,X1,X2,X3,(V(II),II=1,ID)
c2010    FORMAT(1X,I5,3E12.5,5X,6('(',2E12.5,')'))
2010  FORMAT(1X,I5,3E12.5,5X,6(2E12.5))

c          IF (AMPL.GT.EPS) THEN
c            WRITE(MP,2015)IN,X1,X2,X3,V(1),AMPL,DBAMPL,PHASE,STAR
c          ELSE
c            WRITE(MP,2020)IN,X1,X2,X3,V(1),AMPL,STAR
c          ENDIF
c2015	  FORMAT(1X,I5,3E12.5,5X,5E12.5,5x,A4)
c2020	  FORMAT(1X,I5,3E12.5,5X,3E12.5,24('-'),5x,A4)
        ENDIF
      ENDDO

C..........CALCULATE AND PRINT VELOCITIES IF REQUIRED
      IF (VELPAR.EQ.100) THEN
      CALL KAGVEL(VA(LI100KM),VA(LI100T),VA(LI100Q),VA(LI100R),
     1    VA(LI100A),VA(LI100V),VA(LI100DM1),VA(LI100DM2),NI100MOD,
     2    NI100NOD,VA(LFREQ),VA(LPREE),VA(LDIMP),VA(LLOCE),
     3    VA(LCORE),VA(LPRNE),VA(LLD),VA(LNE),VA(LSRC))
      ENDIF

C..........calculate and print intensities if required
      IF (INTPAR.GT.0) THEN
c	CALL INT03()
      ENDIF


      RETURN
      END
