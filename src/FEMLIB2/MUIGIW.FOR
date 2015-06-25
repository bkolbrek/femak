      SUBROUTINE MUIGIW(VIGS,VIGD,VIGI,VSRC,KLD,EQNNO,PIN)
C*******************************************************************C
C     SRCTYP50.FOR                                                  C
C*******************************************************************C
C     22/3-1995                                                     C
C     NTH                                                           C
C     Tonni F. Johansen                                             C
C                                                                   C
C*******************************************************************C
C     TO ASSEMBLE SOURCE VECTOR FROM THE CONTRIBUTIONS FROM A PLANE C
C     WAVE.                                                         C
C*******************************************************************C
C
C
      COMPLEX*16 VIGS(1),VIGD(1),VIGI(1),VSRC(1)
      COMPLEX*16 PIN
      INTEGER KLD(1)
      INTEGER EQNNO

      include 'femak.cmn'

c  $declare
      INTEGER I,I3,COLH

C...........Find diagonal contribution
      VSRC(EQNNO)=VSRC(EQNNO)-VIGD(EQNNO)*PIN
C...........find contribution from upper triangle           
      COLH=KLD(EQNNO+1)-KLD(EQNNO)
c      I1=EQNNO-COLH
c      I2=EQNNO-1        
      I3=KLD(EQNNO)
c      DO I=I1,I2             
c        VSRC(I)=VIGS(I3+I2-I1)*PIN
      DO I=1,COLH
        VSRC(EQNNO-COLH+I-1)=VSRC(EQNNO-COLH+I-1)-VIGS(I3+I-1)*PIN
      ENDDO


C...........find contribution from lower triangle           
C      IF (NSYM.EQ.0) THEN
      
        I3=KLD(EQNNO)
c      DO I=I1,I2             
c        VSRC(I)=VIGS(I3+I2-I1)*PIN
        DO I=EQNNO+1,NEQ
          COLH=KLD(I+1)-KLD(I)
          IF (COLH.GE.I-EQNNO) THEN
            VSRC(I)=VSRC(I)-VIGS(KLD(I)+COLH-I+EQNNO)*PIN
          ENDIF
        ENDDO

C      ENDIF                                     
      


        
      RETURN
      END