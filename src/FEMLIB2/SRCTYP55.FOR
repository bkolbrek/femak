      SUBROUTINE SRCTYPE55(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,
     1  VKE1,VIE,VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     2  VMGI,VSRC,VDLE,VRES,VNDSRC,
     3  KTYPRSG,VPRSG,KTYPRIG,VPRIG,KNEQ,VPRNG,VCORG,KDLNC,SRCNR)
C*******************************************************************C
C     SRCTYP55.FOR                                                  C
C*******************************************************************C
C     2/5-1995                                                      C
C     NTH                                                           C
C     Tonni F. Johansen                                             C
C                                                                   C
C*******************************************************************C
C     TO ASSEMBLE SOURCE VECTOR FROM THE CONTRIBUTIONS FROM A PLANE C
C     WAVE AT JUNCTION BETWEEN WAVE ENVELOPE ELEMENT AND COMMON     C
C     ACOUSTIC ELEMENT (TYPE 3)                                     C
C*******************************************************************C
C
C
      COMPLEX*16 VIE(1),VIE1(1),VSRE(1),VIGS(1),VIGD(1),VIGI(1),VSRC(1)
      COMPLEX*16 VNDSRC(1),VPRSG(1),VPRIG(1)
      REAL*8 VDIMP(1),VCORE(1),VPRNE(1),VPREE(1),VKE(1),VFE(1),VDLE(1)
      REAL*8 VRES(1),VKE1(1),VPRNG(1),VCORG(1)
      REAL*8 VKGS(1),VKGD(1),VKGI(1),VMGS(1),VMGD(1),VMGI(1)
      INTEGER KLD(1),KLOCE(1),KNE(1),KTYPRSG(1),KNEQ(1),KTYPRIG(1)
      INTEGER KDLNC(1)
      REAL*8 ARG

      include 'femak.cmn'
      COMPLEX*16 VDLEC(1)
c $declare
      COMPLEX*16 A,PIN(3),VEL(3),IMP
      REAL*8 KX,KY,RHOC

      INTEGER I,IE,GRPNO,FNDWV,NSIDES,SIDE(6)

      ARG = VPRSG((IGSRC-1)*NFREQ+IFREQ)
      GRPNO=NINT(ARG)
      DO I=1,NNT
        ARG = VPRNG(NPRN*(I-1)+1)
        IF (GRPNO.EQ.NINT(ARG)) THEN
C.........find pressures p_i
          A=VNDSRC((IFREQ-1)*4+1)
          IF (NDIM.EQ.2) THEN

            KX=REAL(VNDSRC((IFREQ-1)*4+2))
            KY=REAL(VNDSRC((IFREQ-1)*4+3))
            PIN(1)=A*EXP(CMPLX(0.0,-W/VPREE(2))*(VCORG(NDIM*I-1)*KX+
     1         VCORG(NDIM*I)*KY))
          ENDIF
C     CARRY OUT multiplications, note only for KG and w**2MG
          CALL MUIGIW55(VKGS,VKGD,VKGI,VMGS,VMGD,VMGI,VSRC,KLD,
     1             KDLNC(I+1),PIN(1))

        ENDIF
      ENDDO
C--------Finished loop through all elements





C...........Loop through all elements to find contributions to source-terms
      REWIND M1
      ITPE1=0
      IFL1=0
      DO IE=1,NELTR

        FNDWV=0
C------ READ AN IMPEDANCE INFO ON FILE M6
        CALL RDELEM(M1,KLOCE,VCORE,VPRNE,VPREE,KNE)

C..........2D VERSION
        NSIDES=0
        IF ((GRPNO.EQ.NINT(VPRNG(KNE(7)))).AND.
     1           (GRPNO.EQ.NINT(VPRNG(KNE(8)))).AND.
     2           (GRPNO.EQ.NINT(VPRNG(KNE(1))))) THEN
          NSIDES=NSIDES+1
          SIDE(NSIDES)=2
        ENDIF
        IF ((GRPNO.EQ.NINT(VPRNG(KNE(1)))).AND.
     1           (GRPNO.EQ.NINT(VPRNG(KNE(2)))).AND.
     2           (GRPNO.EQ.NINT(VPRNG(KNE(3))))) THEN
          NSIDES=NSIDES+1
          SIDE(NSIDES)=4
        ENDIF
        IF ((GRPNO.EQ.NINT(VPRNG(KNE(3)))).AND.
     1           (GRPNO.EQ.NINT(VPRNG(KNE(4)))).AND.
     2           (GRPNO.EQ.NINT(VPRNG(KNE(5))))) THEN
          NSIDES=NSIDES+1
          SIDE(NSIDES)=3
        ENDIF
        IF ((GRPNO.EQ.NINT(VPRNG(KNE(5)))).AND.
     1           (GRPNO.EQ.NINT(VPRNG(KNE(6)))).AND.
     2           (GRPNO.EQ.NINT(VPRNG(KNE(7))))) THEN
          NSIDES=NSIDES+1
          SIDE(NSIDES)=5
        ENDIF

C.......Loop through all sides of the elements with "incoming wave sources"
        DO I=1,NSIDES
          IF (NDIM.EQ.2) THEN
            IF(SIDE(I).EQ.2) THEN
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
            ELSEIF(SIDE(I).EQ.3) THEN
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
            ELSEIF(SIDE(I).EQ.4) THEN
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
            ELSEIF(SIDE(I).EQ.5) THEN
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

C.......Find the new contributions, that is the direction of the wave
          IFL=SIDE(I)
          RHOC=VPREE(1)*VPREE(2)
          VEL(1)=VNDSRC((IFREQ-1)*4+1)
          VEL(2)=VNDSRC((IFREQ-1)*4+2)
          VEL(3)=VNDSRC((IFREQ-1)*4+3)

          VNDSRC((IFREQ-1)*4+1)=PIN(1)/RHOC
          VNDSRC((IFREQ-1)*4+2)=PIN(2)/RHOC
          VNDSRC((IFREQ-1)*4+3)=PIN(3)/RHOC
          VNDSRC((IFREQ-1)*4+4)=CMPLX(KX,KY)
          ICODE=15
          CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)

          VNDSRC((IFREQ-1)*4+1)=VEL(1)
          VNDSRC((IFREQ-1)*4+2)=VEL(2)
          VNDSRC((IFREQ-1)*4+3)=VEL(3)
          CALL ASSELC(0,1,IDLE,NSYM,KLOCE,KLD,VIE,VSRE,VIGS,VIGD,VIGI,
     1            VSRC)
        ENDDO
      ENDDO

C      DO FOR ALL IMPEDANCES
C        IF IMPEDANCE HAS NODE IN WANTES SET:
C          CARRY OUT CHANGES OF CONTRIBUTION IN IMPEDANCE
C        ENDIF
C      ENDDO

      RETURN
      END
