      SUBROUTINE SRCTYPE50(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,
     1  VKE1,VIE,VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VSRC,VDLE,VRES,VNDSRC,
     2  KTYPRSG,VPRSG,KTYPRIG,VPRIG,KNEQ,VPRNG,VCORG,KDLNC,SRCNR)
C********************************************************************C
C     SRCTYP50.FOR                                                   C
C********************************************************************C
C     22/3-1995                                                      C
C     NTH                                                            C
C     Tonni F. Johansen                                              C
C                                                                    C
C********************************************************************C
C     TO ASSEMBLE SOURCE VECTOR FROM THE CONTRIBUTIONS FROM A PLANE  C
C     WAVE.                                                          C
C********************************************************************C
C
C
      COMPLEX*16 VIE(1),VIE1(1),VSRE(1),VIGS(1),VIGD(1),VIGI(1),VSRC(1)
      COMPLEX*16 VNDSRC(1),VPRSG(1),VPRIG(1)
      REAL*8 VDIMP(1),VCORE(1),VPRNE(1),VPREE(1),VKE(1),VFE(1),VDLE(1)
      REAL*8 VRES(1),VKE1(1),VPRNG(1),VCORG(1)
      INTEGER KLD(1),KLOCE(1),KNE(1),KTYPRSG(1),KNEQ(1),KTYPRIG(1)
      INTEGER KDLNC(1)

      include 'femak.cmn'
      COMPLEX*16 VDLEC(1)
C  $declare
      COMPLEX*16 A,PIN(3),VEL(3),IMP
      REAL*8 KX,KY,RHOC
      
      INTEGER I,IE,GRPNO,FNDWV


      GRPNO=NINT(REAL(VPRSG((IGSRC-1)*NFREQ+IFREQ)))
      DO I=1,NNT
        IF (GRPNO.EQ.NINT(REAL(VPRNG(NPRN*(I-1)+1)))) THEN
C.........find pressures p_i
          A=VNDSRC((IFREQ-1)*4+1)
          IF (NDIM.EQ.2) THEN
          
            KX=REAL(VNDSRC((IFREQ-1)*4+2))
            KY=REAL(VNDSRC((IFREQ-1)*4+3))
            PIN(1)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORG(NDIM*I-1)*KX+
     1         VCORG(NDIM*I)*KY))
          ENDIF
C     CARRY OUT multiplications        
          CALL MUIGIW(VIGS,VIGD,VIGI,VSRC,KLD,KNEQ(KDLNC(I+1)),PIN(1))
          
        ENDIF
      ENDDO
          

      REWIND M6
      ITPE1=0
      IFL1=0
      DO IE=1,NSIMP

        FNDWV=0
C------ READ AN IMPEDANCE INFO ON FILE M6
        CALL RDCIMP(M6,KLOCE,VCORE,VPRNE,VPREE,KNE,IGIMP)

C.......Check whether its nodes are in the surface set.
        IF (NDIM.EQ.2) THEN
          IF(IFL.EQ.2) THEN
            IF ((GRPNO.EQ.NINT(VPRNG(KNE(7)))).OR.
     1           (GRPNO.EQ.NINT(VPRNG(KNE(8)))).OR.
     2           (GRPNO.EQ.NINT(VPRNG(KNE(1))))) THEN
              FNDWV=1
              A=VNDSRC((IFREQ-1)*4+1)
              KX=REAL(VNDSRC((IFREQ-1)*4+2))
              KY=REAL(VNDSRC((IFREQ-1)*4+3))
              PIN(1)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORE(13)*KX+
     1          VCORE(14)*KY))
              PIN(2)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORE(15)*KX+
     1          VCORE(16)*KY))
              PIN(3)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORE(1)*KX+
     1          VCORE(2)*KY))
            ENDIF
         ELSEIF(IFL.EQ.3) THEN
            IF ((GRPNO.EQ.NINT(VPRNG(KNE(3)))).OR.
     1           (GRPNO.EQ.NINT(VPRNG(KNE(4)))).OR.
     2           (GRPNO.EQ.NINT(VPRNG(KNE(5))))) THEN
              FNDWV=1
              A=VNDSRC((IFREQ-1)*4+1)
              KX=REAL(VNDSRC((IFREQ-1)*4+2))
              KY=REAL(VNDSRC((IFREQ-1)*4+3))
              PIN(1)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORE(5)*KX+
     1          VCORE(6)*KY))
              PIN(2)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORE(7)*KX+
     1          VCORE(8)*KY))
              PIN(3)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORE(9)*KX+
     1          VCORE(10)*KY))
            ENDIF
         ELSEIF(IFL.EQ.4) THEN
            IF ((GRPNO.EQ.NINT(VPRNG(KNE(1)))).OR.
     1           (GRPNO.EQ.NINT(VPRNG(KNE(2)))).OR.
     2           (GRPNO.EQ.NINT(VPRNG(KNE(3))))) THEN
              FNDWV=1
              A=VNDSRC((IFREQ-1)*4+1)
              KX=REAL(VNDSRC((IFREQ-1)*4+2))
              KY=REAL(VNDSRC((IFREQ-1)*4+3))
              PIN(1)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORE(1)*KX+
     1          VCORE(2)*KY))
              PIN(2)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORE(3)*KX+
     1          VCORE(4)*KY))
              PIN(3)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORE(5)*KX+
     1          VCORE(6)*KY))
            ENDIF
          ELSEIF(IFL.EQ.5) THEN
            IF ((GRPNO.EQ.NINT(VPRNG(KNE(5)))).OR.
     1           (GRPNO.EQ.NINT(VPRNG(KNE(6)))).OR.
     2           (GRPNO.EQ.NINT(VPRNG(KNE(7))))) THEN
              FNDWV=1
              A=VNDSRC((IFREQ-1)*4+1)
              KX=REAL(VNDSRC((IFREQ-1)*4+2))
              KY=REAL(VNDSRC((IFREQ-1)*4+3))
              PIN(1)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORE(9)*KX+
     1          VCORE(10)*KY))
              PIN(2)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORE(11)*KX+
     1          VCORE(12)*KY))
              PIN(3)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORE(13)*KX+
     1          VCORE(14)*KY))
            ENDIF
          ENDIF
        ENDIF  
      
        IF (FNDWV.EQ.1) THEN
C........Eliminate contributions first

C........Original members of VNDSRC are stored
          VEL(1)=VNDSRC((IFREQ-1)*4+1)
          VEL(2)=VNDSRC((IFREQ-1)*4+2)
          VEL(3)=VNDSRC((IFREQ-1)*4+3)
                        
C...........Find "equivalent normal velocity"
          RHOC=VPREE(1)*VPREE(2)
          IMP=VPRIG((IGIMP-1)*NFREQ+IFREQ)*RHOC
          VNDSRC((IFREQ-1)*4+1)=-PIN(1)/IMP        
          VNDSRC((IFREQ-1)*4+2)=-PIN(2)/IMP
          VNDSRC((IFREQ-1)*4+3)=-PIN(3)/IMP
          
          ICODE=2
          CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1          VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
          ICODE=10
          CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)

          CALL ASSELC(0,1,IDLE,NSYM,KLOCE,KLD,VIE,VSRE,VIGS,VIGD,VIGI,
     1            VSRC)        

C.......Find the new contributions, that is the direction of the wave
          VNDSRC((IFREQ-1)*4+1)=PIN(1)/RHOC
          VNDSRC((IFREQ-1)*4+2)=PIN(2)/RHOC
          VNDSRC((IFREQ-1)*4+3)=PIN(3)/RHOC
          VNDSRC((IFREQ-1)*4+4)=CMPLX(KX,KY)
          ICODE=15
          CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)

          CALL ASSELC(0,1,IDLE,NSYM,KLOCE,KLD,VIE,VSRE,VIGS,VIGD,VIGI,
     1            VSRC)        


C........find velocities
C...NB! might be dangerous trick: go through all earlier sources for 
C.......changes of source (type 20 and 21)        

C...........1. change other sources, sources
          REWIND(M7)
C...........Read all previously defined sources
          DO I=1,SRCNR-1
            CALL RDSRCE(M7,KLOCE,VCORE,VPRNE,VPREE,KNE,VNDSRC,KTYPRSG)

            IF ((KTYPRSG(IGSRC).EQ.20).OR.(KTYPRSG(IGSRC).EQ.21)) THEN
C............only sorce types 20 or 21 are changed
              A=VNDSRC((IFREQ-1)*4+1)
              KX=REAL(VNDSRC((IFREQ-1)*4+2))
              KY=REAL(VNDSRC((IFREQ-1)*4+3))
              IF(IFL.EQ.2) THEN
                IF ((GRPNO.EQ.NINT(VPRNG(KNE(7)))).OR.
     1           (GRPNO.EQ.NINT(VPRNG(KNE(8)))).OR.
     2           (GRPNO.EQ.NINT(VPRNG(KNE(1))))) THEN
                  VNDSRC((IFREQ-1)*4+1)=A*EXP(CMPLX(0.0,-W/VPREE(2))
     1                    *(VCORE(13)*KX+VCORE(14)*KY))/RHOC
                  VNDSRC((IFREQ-1)*4+2)=A*EXP(CMPLX(0.0,-W/VPREE(2))
     1                    *(VCORE(13)*KX+VCORE(14)*KY))/RHOC
                  VNDSRC((IFREQ-1)*4+1)=A*EXP(CMPLX(0.0,-W/VPREE(2))
     1                    *(VCORE(13)*KX+VCORE(14)*KY))/RHOC
                  VNDSRC((IFREQ-1)*4+4)=CMPLX(KX,KY)
                ENDIF
              ELSEIF(IFL.EQ.3) THEN
                IF ((GRPNO.EQ.NINT(VPRNG(KNE(3)))).OR.
     1             (GRPNO.EQ.NINT(VPRNG(KNE(4)))).OR.
     2             (GRPNO.EQ.NINT(VPRNG(KNE(5))))) THEN
                  VNDSRC((IFREQ-1)*4+1)=A*EXP(CMPLX(0.0,-W/VPREE(2))*
     1                (VCORE(5)*KX+VCORE(6)*KY))/RHOC
                  VNDSRC((IFREQ-1)*4+2)=A*EXP(CMPLX(0.0,-W/VPREE(2))*
     1                (VCORE(7)*KX+VCORE(8)*KY))/RHOC
                  VNDSRC((IFREQ-1)*4+3)=A*EXP(CMPLX(0.0,-W/VPREE(2))*
     1                (VCORE(9)*KX+VCORE(10)*KY))/RHOC
                  VNDSRC((IFREQ-1)*4+4)=CMPLX(KX,KY)
                ENDIF
              ELSEIF(IFL.EQ.4) THEN
                IF ((GRPNO.EQ.NINT(VPRNG(KNE(1)))).OR.
     1            (GRPNO.EQ.NINT(VPRNG(KNE(2)))).OR.
     2            (GRPNO.EQ.NINT(VPRNG(KNE(3))))) THEN
                  VNDSRC((IFREQ-1)*4+1)=A*EXP(CMPLX(0.0,-W/VPREE(2))*
     1                (VCORE(1)*KX+VCORE(2)*KY))/RHOC
                  VNDSRC((IFREQ-1)*4+2)=A*EXP(CMPLX(0.0,-W/VPREE(2))*
     1                (VCORE(3)*KX+VCORE(4)*KY))/RHOC
                  VNDSRC((IFREQ-1)*4+3)=A*EXP(CMPLX(0.0,-W/VPREE(2))*
     1                (VCORE(5)*KX+VCORE(6)*KY))/RHOC
                  VNDSRC((IFREQ-1)*4+4)=CMPLX(KX,KY)
                ENDIF
              ELSEIF(IFL.EQ.5) THEN
                IF ((GRPNO.EQ.NINT(VPRNG(KNE(5)))).OR.
     1              (GRPNO.EQ.NINT(VPRNG(KNE(6)))).OR.
     2              (GRPNO.EQ.NINT(VPRNG(KNE(7))))) THEN
                  VNDSRC((IFREQ-1)*4+1)=A*EXP(CMPLX(0.0,-W/VPREE(2))*
     1              (VCORE(9)*KX+VCORE(10)*KY))/RHOC
                  VNDSRC((IFREQ-1)*4+2)=A*EXP(CMPLX(0.0,-W/VPREE(2))*
     1              (VCORE(11)*KX+VCORE(12)*KY))/RHOC
                  VNDSRC((IFREQ-1)*4+3)=A*EXP(CMPLX(0.0,-W/VPREE(2))*
     1              (VCORE(13)*KX+VCORE(14)*KY))/RHOC
                  VNDSRC((IFREQ-1)*4+4)=CMPLX(KX,KY)
                ENDIF
              ENDIF
              ICODE=15
              CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,
     1            VIE1,VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)

              CALL ASSELC(0,1,IDLE,NSYM,KLOCE,KLD,VIE,VSRE,VIGS,VIGD,
     1            VIGI,VSRC)        
            ENDIF
          ENDDO



C........call new type of source (takes any velocities and return contributions from
C........normal component v.hj.a. elem03
        
C........Restore Original members of VNDSRC 
          VNDSRC((IFREQ-1)*4+1)=VEL(1)
          VNDSRC((IFREQ-1)*4+2)=VEL(2)
          VNDSRC((IFREQ-1)*4+3)=VEL(3)
        ENDIF
      ENDDO

C      DO FOR ALL IMPEDANCES
C        IF IMPEDANCE HAS NODE IN WANTES SET:
C          CARRY OUT CHANGES OF CONTRIBUTION IN IMPEDANCE                
C        ENDIF                
C      ENDDO

      RETURN
      END
      