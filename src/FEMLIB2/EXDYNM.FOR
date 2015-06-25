      SUBROUTINE EXDYNM(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1   VFE,VKGS,VKGD,VKGI,VMGS,VMGD,VMGI,VIGS,VIGD,VIGI,VFG,VSRC,
     2   VCORG,KDLNC,KNEQ,VRES,VDLE,VIE,VIE1,VSRE,VFREQ,KTYPRIG,VPRIG,
     3   VNDSRC,KTYPRSG,VPRSG,VPRNG)
C***********************************************************************
C     EXDYNM.FOR
C***********************************************************************
C     3/1-1995
C     NTH
C     Tonni F. Johansen  & Ulf R. Kristiansen
C***********************************************************************
C     TO EXECUTE BLOCK 'DYNM'
C     ASSEMBLE AND SOLVE A DYNAMIC PROBLEM IN CORE
c     Designed primarely for acoustic problems.
C***********************************************************************
      REAL*8 PI
      PARAMETER (PI=3.1415926)

      include 'femak.cmn'
C $DECLARE      
      REAL*8 VDIMP(1),VCORE(1),VPRNE(1),VPREE(1),VKE(1),VFE(1),VFG(1)
      REAL*8 VKGS(1),VKGD(1),VKGI(1),VMGS(1),VMGD(1),VMGI(1)
      REAL*8 VCORG(1),VRES(1),VDLE(1),VFREQ(1),VKE1(1),VPRNG(1)
      COMPLEX*16 VIGS(1),VIGD(1),VIGI(1),VSRC(1),VIE(1),VSRE(1)
      COMPLEX*16 VPRIG(1),VNDSRC(1),VPRSG(1),VIE1(1)
      INTEGER KLD(1),KLOCE(1),KNE(1),KDLNC(1),KNEQ(1),KTYPRSG(1)
      INTEGER KTYPRIG(1)

      INTEGER I,J,FREQPAR,I2,J2
      REAL*8 ENERG
      COMPLEX*16 CZERO
      DATA CZERO/(0.0,0.0)/
      
C------------------------------------------------------------

C.........Read frequencies
      READ(MR,1000)FREQPAR
1000  FORMAT(I5)

      IF (FREQPAR.EQ.0) THEN
        IF (MOD(NFREQ,10).EQ.0) THEN
          I2=NFREQ/10
            ELSE
          I2=NFREQ/10+1
        ENDIF
        DO  I=1,I2
          J2=MIN(10,NFREQ-(I-1)*10)
          READ(MR,1010)(VFREQ((I-1)*10+J),J=1,J2)
1010      FORMAT(10E12.5)
          IF(M.GT.0) THEN
            WRITE(MP,1020)(VFREQ((I-1)*10+J),J=1,J2)
1020        FORMAT(' FREQUENCIES :',10E12.5)
          ENDIF
        ENDDO
      ENDIF

C initialisation pour DSP
      kxpas = 30.11111D0
      kxinc = 1.333333D0
      radis = 0.D0      

C-------- LOOP THROUGH ALL FREQUENCIES
      DO IFREQ=1,NFREQ
C-----------Set system matrix to zero
        DO I=1,NKG
          VIGS(I)=CZERO
        ENDDO

        DO I=1,NEQ
          VIGD(I)=CZERO
        ENDDO

        IF (NSYM.EQ.1) THEN
          DO I=1,NKG
            VIGI(I)=CZERO
          ENDDO
        ENDIF

        DO I=1,NEQ
          VSRC(I)=CZERO
        ENDDO



C---------Calculate the frequency (or table)
        W=2*PI*VFREQ(IFREQ)
        W2=W*W

        IF (UPDTFRQ.OR.(IFREQ.EQ.1)) THEN

C
C-------- ASSEMBLE KG AND MG
          CALL ASKG(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1      VIE,VIE1,VFE,VKGS,VKGD,VKGI,VFG,VDLE,VRES,
     2      VNDSRC,KTYPRSG,VPRSG)
          CALL ASMG(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1      VIE,VIE1,VFE,VMGS,VMGD,VMGI,VFG,VDLE,VRES,
     2      VNDSRC,KTYPRSG,VPRSG)
C
        ENDIF
        IF (CMAT) THEN
C
C-------- ASSEMBLE COMPLEX KG AND MG
          CALL ASKGC(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1      VIE,VIE1,VFE,VSRE,VKGS,VKGD,VKGI,VIGS,VIGD,VIGI,VFG,VSRC,
     2      VDLE,VRES,VNDSRC,KTYPRSG,VPRSG)
          CALL ASMGC(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1      VIE,VIE1,VFE,VSRE,VKGS,VKGD,VKGI,VIGS,VIGD,VIGI,VFG,VSRC,
     2      VDLE,VRES,VNDSRC,KTYPRSG,VPRSG)
C
        ENDIF
C---------ASSEMBLE IG(IMPEDANCEMATRIX)
        CALL ASIG(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1  VIE,VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VSRC,VDLE,VRES,VNDSRC,
     2  KTYPRIG,VPRIG,VFREQ)
      
C---------ASSEMBLE LOSSMATRIX
        CALL ASLGC(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1    VIE,VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VSRC,VDLE,VRES,VNDSRC,
     2    KTYPRSG,VPRSG,VFREQ)
     
C---------ASSEMBLE SYSTEMMATRIX  IG(I,J)=KG-w**2*MG+IG(I,J)
c        write(mp,*)'vigs kun impedans '
c        do i=1,nkg
c        write(mp,*)i,vigs(i)
c        enddo

        DO I=1,NKG
          VIGS(I)=VKGS(I)-W2*VMGS(I)+VIGS(I)
        ENDDO
c        write(mp,*)'vigs fra total'
c        do i=1,nkg
c        write(mp,*)i,vigs(i)
c        enddo

c        write(mp,*)'vigd kun fra impedans'
c        do i=1,neq
c        write(mp,*)i,vigd(i)
c        enddo
        DO I=1,NEQ
          VIGD(I)=VKGD(I)-W2*VMGD(I)+VIGD(I)
        ENDDO
c        write(mp,*)'vkgd'
c        do i=1,neq
c        write(mp,*)i,vigd(i)
c        enddo

        IF (NSYM.EQ.1) THEN
          DO I=1,NKG
            VIGI(I)=VKGI(I)-W2*VMGI(I)+VIGI(I)
          ENDDO
c        print*,'vkgi fra exdynm'
c        do i=1,nkg
c        print*,i,vkgi(i)
c        enddo
        ENDIF

C-----------Assemble source vector
cMARC ======================================
C calcul de densite spectrale de puissance ?      
       IF (excit.NE.0) THEN

        CALL DSP(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,VIE,
     1           VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     2           VMGI,VSRC,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG,KTYPRIG,
     3           VPRIG,KNEQ,VPRNG,VCORG,KDLNC,ENERG,VFREQ)
     
cMARC ======================================     
       ELSE

        CALL ASRC(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1    VIE,VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,VMGI
     2    ,VSRC,VDLE,VRES,VNDSRC,
     3    KTYPRSG,VPRSG,KTYPRIG,VPRIG,KNEQ,VPRNG,VCORG,KDLNC)

C......Kutter ut mulitplikqsjon med jw
c        DO I=1,NEQ
c           if (m.ge.2) then
c             write(*,*)'vsrc=',vsrc(i)
c           endif
c          VSRC(I)=DCMPLX(0.,-W)*VSRC(I)
c        ENDDO


C
C--------   SOLVE  COMPLEKS
c        write(MP,*)'utskrift fra EXDYNM'
c        write(MP,*)'vigs:'
c        write(MP,*)(vigs(i),i=1,nkg)
C
c        write(MP,*)'vigd:'
c        write(MP,*)(vigd(i),i=1,neq)
C
c        write(MP,*)'vigi:'
c        write(MP,*)(vigi(i),i=1,nkg)
C
c       write(MP,*)'vsrc:'
c       write(MP,*)(vsrc(i),i=1,neq)

        CALL SOLC(VIGS,VIGD,VIGI,VSRC,KLD,NEQ,MP,1,1,NSYM,
     1        ENERG)

C        write(MP,*)'vsrc:'
C        write(MP,*)(vsrc(i),i=1,neq)


C...........PRINT SOLUTIONS
        CALL PRDYNM(KDLNC,VCORG,VDIMP,KNEQ,VSRC,VFREQ,KLOCE,VCORE,
     1     VPRNE,VPREE,KNE,VNDSRC,KTYPRSG,VPRSG,VPRNG)

       ENDIF

      ENDDO



      RETURN
      END

