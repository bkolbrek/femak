      SUBROUTINE MUIGIW55(VKGS,VKGD,VKGI,VMGS,VMGD,VMGI,VSRC,KLD,
     1            EQNNO,PIN)
C*******************************************************************C
C     MUIGIW55.FOR                                                  C
C*******************************************************************C
C     2/5-1995                                                      C
C     NTH                                                           C
C     Tonni F. Johansen                                             C
C                                                                   C
C*******************************************************************C
C     TO ASSEMBLE SOURCE VECTOR FROM THE CONTRIBUTIONS FROM A PLANE C
C     WAVE, SPECIAL VERSION FOR SRCTYP 55                           C
C*******************************************************************C
C
C
      COMPLEX*16 VSRC(1)
      COMPLEX*16 PIN
      REAL*8 VKGS(1),VKGD(1),VKGI(1),VMGS(1),VMGD(1),VMGI(1)
      INTEGER KLD(1)
      INTEGER EQNNO

      include 'femak.cmn'

c  $declare
      INTEGER I,I3,COLH

C...........Find diagonal contribution
      VSRC(EQNNO)=VSRC(EQNNO)-(VKGD(EQNNO)-W2*VMGD(EQNNO))*PIN
C...........find contribution from upper triangle           
      COLH=KLD(EQNNO+1)-KLD(EQNNO)
c      I1=EQNNO-COLH
c      I2=EQNNO-1        
      I3=KLD(EQNNO)
c      DO I=I1,I2             
c        VSRC(I)=VIGS(I3+I2-I1)*PIN
      DO I=1,COLH
        VSRC(EQNNO-COLH+I-1)=VSRC(EQNNO-COLH+I-1)-
     1           (VKGS(I3+I-1)-W2*VMGS(I3+I-1))*PIN
      ENDDO


C...........find contribution from lower triangle           
      IF (NSYM.EQ.0) THEN
      
        I3=KLD(EQNNO)
c      DO I=I1,I2             
c        VSRC(I)=VIGS(I3+I2-I1)*PIN
        DO I=EQNNO+1,NDLT
          COLH=KLD(I+1)-KLD(I)
          IF (COLH.GE.I-EQNNO) THEN
            VSRC(I)=VSRC(I)-(VKGS(KLD(I)+COLH-I+EQNNO)-
     1            W2*VMGS(KLD(I)+COLH-I+EQNNO))*PIN
          ENDIF
        ENDDO

      ELSE  

c!  nb......MUST BE CHANGED. dOES NOT TAKE ASYMMETRY INTO ACCOUNT.
        I3=KLD(EQNNO)
c      DO I=I1,I2             
c        VSRC(I)=VIGS(I3+I2-I1)*PIN
        DO I=EQNNO+1,NDLT
          COLH=KLD(I+1)-KLD(I)
          IF (COLH.GE.I-EQNNO) THEN
            VSRC(I)=VSRC(I)-(VKGS(KLD(I)+COLH-I+EQNNO)-
     1            W2*VMGS(KLD(I)+COLH-I+EQNNO))*PIN
          ENDIF
        ENDDO

      ENDIF                                     
      


        
      RETURN
      END
      