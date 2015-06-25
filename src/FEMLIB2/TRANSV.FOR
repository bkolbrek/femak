      SUBROUTINE TRANSV(kxinf,kxsup,Ndsp,Nbi,
     1                  Vrep,Vdsp,Vdsi,Vcte,Ures,Vspa,Vsim,Verr,
     2                  VGrep,VGdsp,VGdsi,VGcte,UGres,VGspa,VGsim,
     3            KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,VIE,
     4            VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     5            VMGI,VSRC,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG,KTYPRIG,
     6            VPRIG,KNEQ,VPRNG,VCORG,KDLNC,ENERG,VFREQ)

      include 'femak.cmn'

c      COMMON/CMDSP/dir_xy,excit,idsp1,idsp2,idsp3,idsp4,idsp5,
c     +             vit_son1,rhof1,vit_infi,del_cl,kspace


      COMPLEX*16 VIE(1),VIE1(1),VSRE(1),VIGS(1),VIGD(1),VIGI(1),VSRC(1)
      COMPLEX*16 VNDSRC(1),VPRSG(1),VPRIG(1)
      REAL*8 VDIMP(1),VCORE(1),VPRNE(1),VPREE(1),VKE(1),VFE(1),VDLE(1)
      REAL*8 VRES(1),VKE1(1),VPRNG(1),VCORG(1),VFREQ(1)
      REAL*8 VKGS(1),VKGD(1),VKGI(1),VMGS(1),VMGD(1),VMGI(1)
      INTEGER*4 KLD(1),KLOCE(1),KNE(1),KTYPRSG(1),KNEQ(1),KTYPRIG(1)
      INTEGER*4 KDLNC(1)
      COMPLEX*16 ENERG, ARG


      INTEGER*4   Nbp,Ndsp,Nbi,ii
      REAL*8      ori,ext,pai,imp,del
      REAL*8      kxinf,kxsup,pomme,precision
      REAL*8      Vcte(1),Ures(1),Vspa(1),Vsim(1)
      REAL*8      Vdsi(1),Vdsp(1),Vrep(1),Verr(1)
      REAL*8      VGcte(1),UGres(1),VGspa(1),VGsim(1)
      REAL*8      VGdsi(1),VGdsp(1),VGrep(1)


      PARAMETER (precision = 1.0D-2 )


C INITIALISATION :

      Nbp = 2
      Nbi = Nbp
      ori = kxinf
      ext = kxsup
      del = ext - ori

      CALL INSTAN(ori,del,Nbp,Ndsp,Vdsi,VGdsi,
     1            KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,VIE,
     2            VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     3            VMGI,VSRC,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG,KTYPRIG,
     4            VPRIG,KNEQ,VPRNG,VCORG,KDLNC,ENERG,VFREQ)

      DO ii=1,Ndsp
       Vcte(ii) = Vdsi(ii)
      END DO
      DO ii=1,NNT
       VGcte(ii) = VGdsi(ii)
      END DO

C INTEGRATION A 8 POINTS :

      Nbp = 3
      Nbi = Nbi + Nbp
      del = ( ext - ori ) / DBLE(Nbp+1)
      pai = ori + del

      CALL INSTAN(pai,del,Nbp,Ndsp,Vdsi,VGdsi,
     1            KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,VIE,
     2            VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     3            VMGI,VSRC,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG,KTYPRIG,
     4            VPRIG,KNEQ,VPRNG,VCORG,KDLNC,ENERG,VFREQ)

      DO ii=1,Ndsp
       Vspa(ii) = Vdsi(ii)
      END DO
      DO ii=1,NNT
       VGspa(ii) = VGdsi(ii)
      END DO

      Nbp = 4
      Nbi = Nbi + Nbp
      imp = ori + del / 2.D0

      CALL INSTAN(imp,del,Nbp,Ndsp,Vdsi,VGdsi,
     1            KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,VIE,
     2            VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     3            VMGI,VSRC,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG,KTYPRIG,
     4            VPRIG,KNEQ,VPRNG,VCORG,KDLNC,ENERG,VFREQ)

      DO ii=1,Ndsp
       Vsim(ii) = Vdsi(ii)
       Ures(ii) = Vcte(ii) + 2.D0 * Vspa(ii) + 4.D0 * Vsim(ii)
       Ures(ii) = Ures(ii) * del / 6.D0
      END DO
      DO ii=1,NNT
       VGsim(ii) = VGdsi(ii)
       UGres(ii) = VGcte(ii) + 2.D0 * VGspa(ii) + 4.D0 * VGsim(ii)
       UGres(ii) = UGres(ii) * del / 6.D0
      END DO

C INTEGRATION A 2*i*Nb POINTS :

      pomme = 2.D0 * precision
      DO WHILE ((pomme.GT.precision).OR.(Nbp.LT.32))
       Nbp = 2 * Nbp
       Nbi = Nbi + Nbp
       del = ( ext - ori ) / DBLE(Nbp)
       imp = ori + del / 2.D0

       DO ii=1,Ndsp
        Vspa(ii) = Vspa(ii) + Vsim(ii)
       END DO
       DO ii=1,NNT
        VGspa(ii) = VGspa(ii) + VGsim(ii)
       END DO

       CALL INSTAN(imp,del,Nbp,Ndsp,Vdsi,VGdsi,
     1            KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,VIE,
     2            VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     3            VMGI,VSRC,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG,KTYPRIG,
     4            VPRIG,KNEQ,VPRNG,VCORG,KDLNC,ENERG,VFREQ)

       pomme = 0.D0
       DO ii=1,Ndsp
        Vsim(ii) = Vdsi(ii)
        Vrep(ii) = Ures(ii)
        Ures(ii) = Vcte(ii) + 2.D0 * Vspa(ii) + 4.D0 * Vsim(ii)
        Ures(ii) = Ures(ii) * del / 6.D0
        ARG = cmplx( (Ures(ii)-Vrep(ii)) / (Vdsp(ii)+Vrep(ii)) )
        Verr(ii) = CDABS( ARG )
        pomme   = DMAX1(pomme,Verr(ii))
       END DO
       DO ii=1,NNT
        VGsim(ii) = VGdsi(ii)
        VGrep(ii) = UGres(ii)
        UGres(ii) = VGcte(ii) + 2.D0 * VGspa(ii) + 4.D0 * VGsim(ii)
        UGres(ii) = UGres(ii) * del / 6.D0
       END DO

      END DO

      DO ii=1,Ndsp
       Vdsp(ii) = Vdsp(ii) + Ures(ii)
      END DO
      DO ii=1,NNT
       VGdsp(ii) = VGdsp(ii) + UGres(ii)
      END DO

      RETURN
      END
   
