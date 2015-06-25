      SUBROUTINE DSP(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,VIE,
     1            VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     2            VMGI,VSRC,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG,KTYPRIG,
     3            VPRIG,KNEQ,VPRNG,VCORG,KDLNC,ENERG,VFREQ) 

      include 'femak.cmn'
      
c      COMMON/CMDSP/dir_xy,excit,idsp1,idsp2,idsp3,idsp4,idsp5,
c     +             vit_son1,rhof1,vit_infi,del_cl,kspace,radis,kxpas,kxinc 


 
      COMPLEX*16 VIE(1),VIE1(1),VSRE(1),VIGS(1),VIGD(1),VIGI(1),VSRC(1)
      COMPLEX*16 VNDSRC(1),VPRSG(1),VPRIG(1)
      REAL*8 VDIMP(1),VCORE(1),VPRNE(1),VPREE(1),VKE(1),VFE(1),VDLE(1)
      REAL*8 VRES(1),VKE1(1),VPRNG(1),VCORG(1),VFREQ(1)
      REAL*8 VKGS(1),VKGD(1),VKGI(1),VMGS(1),VMGD(1),VMGI(1)
      INTEGER*4 KLD(1),KLOCE(1),KNE(1),KTYPRSG(1),KNEQ(1),KTYPRIG(1)
      INTEGER*4 KDLNC(1)
      COMPLEX*16 ENERG
      
      REAL*8 Vdsp(5),Vant(5)
      REAL*8 Vrep(5),Vdsi(5),Vcte(5),Ures(5),Vspa(5),Vsim(5),Verr(5)
      REAL*8 VGdsp(498)
      REAL*8 VGrep(498),VGdsi(498),VGcte(498)
      REAL*8 UGres(498),VGspa(498),VGsim(498)
     
      REAL*8 kxinf,kxsup,delkx,kxmax,kxmin,kacou1,freq
      REAL*8 epsilon,banane,precision,pi 
      REAL*8 x1,x2,x3
      INTEGER*4 Nbi,Ndsp,ii,navet
      DATA pi/3.14159D0/
      OPEN(99,FILE='marcout.dat')

C initialisation
      Ndsp = 5
      freq = VFREQ(IFREQ)
      kxmax = 250.D0
      kxmin = -50.D0
      epsilon   = 0.001D0
      precision = 0.05D0
      Kacou1    =  2*pi*freq/vit_son1
      DO i=1,Ndsp
         Vdsp(i) = 0.D0
         Vant(i) = 0.D0
         Vrep(i) = 0.D0
         Vdsi(i) = 0.D0
         Vcte(i) = 0.D0
         Ures(i) = 0.D0
         Vspa(i) = 0.D0
         Vsim(i) = 0.D0
         Verr(i) = 0.D0
      ENDDO   
      DO i=1,NNT
         VGdsp(i) = 0.D0
         VGrep(i) = 0.D0
         VGdsi(i) = 0.D0
         VGcte(i) = 0.D0
         UGres(i) = 0.D0
         VGspa(i) = 0.D0
         VGsim(i) = 0.D0
      ENDDO   
      
C      WRITE(99,*) 'FREQUENCE :',freq,'Hz'
      
C------------------------------------------------
C INTEGRATION SUR LE DOMAINE ACOUSTIQUE : 0<kx<ka
C------------------------------------------------

      kxinf = epsilon
      kxsup = Kacou1 - epsilon

c      write(99,*) "Integration Kx/Ka entre",kxinf/kacou1," et "
c     &,kxsup/kacou1

      CALL TRANSV(kxinf,kxsup,Ndsp,Nbi,
     1            Vrep,Vdsp,Vdsi,Vcte,Ures,Vspa,Vsim,Verr,
     2            VGrep,VGdsp,VGdsi,VGcte,UGres,VGspa,VGsim,
     3            KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,VIE,
     4            VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     5            VMGI,VSRC,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG,KTYPRIG,
     6            VPRIG,KNEQ,VPRNG,VCORG,KDLNC,ENERG,VFREQ) 

c      write(99,"(A,I4)") "Nombre d'integrations = ",Nbi
c      write(99,*)
      
      navet = Nbi

      DO ii=1,Ndsp
       Vant(ii) = Vdsp(ii)
      END DO

C-------------------------------------------------
C INTEGRATION SUR LE DOMAINE ACOUSTIQUE : -ka<kx<0
C-------------------------------------------------

      kxinf = - Kacou1 + epsilon
      kxsup = - epsilon

C      write(99,*) "Integration Kx/Ka entre",kxinf/kacou1," et "
C     &,kxsup/kacou1

      CALL TRANSV(kxinf,kxsup,Ndsp,Nbi,
     1            Vrep,Vdsp,Vdsi,Vcte,Ures,Vspa,Vsim,Verr,
     2            VGrep,VGdsp,VGdsi,VGcte,UGres,VGspa,VGsim,
     3            KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,VIE,
     4            VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     5            VMGI,VSRC,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG,KTYPRIG,
     6            VPRIG,KNEQ,VPRNG,VCORG,KDLNC,ENERG,VFREQ) 
                  
c      write(99,*) "Nombre d'integrations = ",Nbi," Pas = "
c      write(99,*)

      navet = navet + Nbi

      banane = 0.D0
      DO ii=1,Ndsp
       Vant(ii) = Vant(ii) + Vdsp(ii)
      END DO

C-----------------------------------------------
C INTEGRATION SUR LE DOMAINE CONVECTIF (kx>ka>0)
C-----------------------------------------------

      delkx = kxpas
      kxsup = kacou1 + epsilon
      banane = precision * 2.0D0
      DO WHILE ((banane.GT.precision).OR.(kxsup.LT.(kxmax*kacou1)))
       kxinf = kxsup
       kxsup = kxinf + delkx * kacou1

c       write(99,*) "Integration Kx/Ka entre",kxinf/kacou1," et "
c     &  ,kxsup/kacou1

      CALL TRANSV(kxinf,kxsup,Ndsp,Nbi,
     1            Vrep,Vdsp,Vdsi,Vcte,Ures,Vspa,Vsim,Verr,
     2            VGrep,VGdsp,VGdsi,VGcte,UGres,VGspa,VGsim,
     3            KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,VIE,
     4            VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     5            VMGI,VSRC,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG,KTYPRIG,
     6            VPRIG,KNEQ,VPRNG,VCORG,KDLNC,ENERG,VFREQ) 
     
c       write(99,*) "Nombre d'integrations = ",Nbi," Pas = ",delkx 
c       write(99,*)

       IF (Nbi.LT.257) delkx = delkx * 2.D0
       IF (Nbi.GT.257) delkx = delkx / 2.D0

       navet = navet + Nbi

       banane = 0.D0
       DO ii=1,Ndsp
        banane  = DMAX1(banane,DABS((Vdsp(ii)-Vant(ii))/Vant(ii)))
        Vant(ii) = Vdsp(ii)
       END DO
      END DO

C----------------------------------------------------
C INTEGRATION SUR LES FREQUENCES NEGATIVES (kx<-ka<0)
C----------------------------------------------------

      delkx = kxpas
      kxinf = -kacou1- epsilon
      banane = precision * 2.0D0
      DO WHILE ((banane.GT.precision).OR.(kxinf.GT.(kxmin*kacou1)))
       kxsup = kxinf
       kxinf = kxsup - delkx * kacou1

c       write(99,*) "Integration Kx/Ka entre",kxinf/kacou1," et "
c     &  ,kxsup/kacou1

      CALL TRANSV(kxinf,kxsup,Ndsp,Nbi,
     1            Vrep,Vdsp,Vdsi,Vcte,Ures,Vspa,Vsim,Verr,
     2            VGrep,VGdsp,VGdsi,VGcte,UGres,VGspa,VGsim,
     3            KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,VIE,
     4            VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     5            VMGI,VSRC,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG,KTYPRIG,
     6            VPRIG,KNEQ,VPRNG,VCORG,KDLNC,ENERG,VFREQ) 

c       write(99,*) "Nombre d'integrations = ",Nbi," Pas = ",delkx
c       write(99,*)

       IF (Nbi.LT.257) delkx = delkx * 2.D0
       IF (Nbi.GT.257) delkx = delkx / 2.D0

       navet = navet + Nbi

       banane = 0.D0
       DO ii=1,Ndsp
        banane = DMAX1(banane,DABS((Vdsp(ii)-Vant(ii))/Vant(ii)))
        Vant(ii) = Vdsp(ii)
       END DO
      END DO

      WRITE(99,*)
     &"Frequence = ",freq,",  nombre de points = ",navet,",  pas = "
     &,kxpas
  
      IF ((navet.GT.radis).AND.(radis.NE.0)) kxinc = - kxinc
      IF ((kxpas+kxinc).GE.0.1D0) kxpas = kxpas + kxinc
      radis = navet


C print results only at the five nodes:
c       WRITE(MP,'(F10.3)') VFREQ(IFREQ)
c       WRITE(MP,'(I5,5X,F15.5)') idsp1,10*DLOG10(Vdsp(1))      
c       WRITE(MP,'(I5,5X,F15.5)') idsp2,10*DLOG10(Vdsp(2))      
c       WRITE(MP,'(I5,5X,F15.5)') idsp3,10*DLOG10(Vdsp(3))      
c       WRITE(MP,'(I5,5X,F15.5)') idsp4,10*DLOG10(Vdsp(4))      
c       WRITE(MP,'(I5,5X,F15.5)') idsp5,10*DLOG10(Vdsp(5))
c       WRITE(MP,*) '=================================================='                      '
c
C print results at all the nodes like in PRDYNM:
      WRITE(MP,1000)VFREQ(IFREQ)
1000  FORMAT(' COMPUTATIONS AT FREQUENCY:',F12.3,'Hz'/72('='))       
      DO II=1,NNT
          x1=vcorg(ndim*(ii-1)+1)
          IF (ndim.GE.2) x2=vcorg(ndim*(ii-1)+2)
          IF (ndim.GE.3) x3=vcorg(ndim*(ii-1)+3)
          WRITE(MP,'(1X,I5,3E12.5,5X,E12.5)') II,X1,X2,X3,VGdsp(II)
      END DO      
C      CLOSE(UNIT=99)

      RETURN
      END