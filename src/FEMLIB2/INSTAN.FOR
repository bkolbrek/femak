      SUBROUTINE INSTAN(ori,del,Nbp,Ndsp,Vdsi,VGdsi,
     1            KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,VIE,
     2            VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     3            VMGI,VSRC,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG,KTYPRIG,
     4            VPRIG,KNEQ,VPRNG,VCORG,KDLNC,ENERG,VFREQ)

      include 'femak.cmn'
      
c      COMMON/CMDSP/dir_xy,excit,idsp1,idsp2,idsp3,idsp4,idsp5,
c     +             vit_son1,rhof1,vit_infi,del_cl,kspace 
     

      REAL*8  Vdsi(1),VGdsi(1)
      COMPLEX*16 VIE(1),VIE1(1),VSRE(1),VIGS(1),VIGD(1),VIGI(1),VSRC(1)
      COMPLEX*16 VNDSRC(1),VPRSG(1),VPRIG(1)
      REAL*8 VDIMP(1),VCORE(1),VPRNE(1),VPREE(1),VKE(1),VFE(1),VDLE(1)
      REAL*8 VRES(1),VKE1(1),VPRNG(1),VCORG(1),VFREQ(1)
      REAL*8 VKGS(1),VKGD(1),VKGI(1),VMGS(1),VMGD(1),VMGI(1)
      INTEGER*4 KLD(1),KLOCE(1),KNE(1),KTYPRSG(1),KNEQ(1),KTYPRIG(1)
      INTEGER*4 KDLNC(1)
      COMPLEX*16 ENERG,V(6)

      REAL*8      CORCOS,CHASE,CHASEM
      REAL*8      ori,del
      REAL*8      densite,freq
      INTEGER*4   ii,kk,Ndsp,Nbp
      INTEGER*4 I2,IN,I1,ID1,ID2,ID,J,JJ
      REAL*8 X1,X2,X3
      COMPLEX*16 CZERO
      REAL*4 STAR,RF,EMPTY

      DATA CZERO/(0.D0,0.D0)/
      DATA RF/4H *  /,EMPTY/4H    / 
cmarc      OPEN(98,FILE='contrib.dat') 
      freq  = vfreq(ifreq)       
C INITIALISATION :
      DO ii=1,Ndsp
       Vdsi(ii) = 0.D0
      END DO
      DO ii=1,NNT
       VGdsi(ii) = 0.D0
      END DO

C RESOLUTION :

c      write(99,"(I4,A)") Nbp," resolutions."

      DO kk=1,Nbp

C initialisation
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
        
        CALL ASIG(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1  VIE,VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VSRC,VDLE,VRES,VNDSRC,
     2  KTYPRIG,VPRIG,VFREQ)
     
        DO I=1,NKG
          VIGS(I)=VKGS(I)-W2*VMGS(I)+VIGS(I)
        ENDDO
        DO I=1,NEQ
          VIGD(I)=VKGD(I)-W2*VMGD(I)+VIGD(I)
        ENDDO
        IF (NSYM.EQ.1) THEN
          DO I=1,NKG
            VIGI(I)=VKGI(I)-W2*VMGI(I)+VIGI(I)
          ENDDO
        ENDIF                               
              
       kspace = ori + DBLE(kk-1) * del

       IF (excit.EQ.1) THEN
         densite = CORCOS(kspace,freq)
       ELSE IF (excit.EQ.2) THEN
         densite = CHASEM(kspace,freq)
       ELSE
         STOP "excitation : CORCOS ou CHASE ?"
       END IF
cmarc       WRITE(98,'(F6.2,1X,E12.5,1X,E12.5)') freq,kspace,densite       
       
       CALL ASRC (KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,VIE,
     1            VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     2            VMGI,VSRC,VDLE,VRES,VNDSRC,KTYPRSG,VPRSG,KTYPRIG,
     3            VPRIG,KNEQ,VPRNG,VCORG,KDLNC)
     
       CALL SOLC (VIGS,VIGD,VIGI,VSRC,KLD,NEQ,MP,1,1,NSYM,ENERG)
       
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
              J=J+1
            ENDDO
                      
            IF (IN.EQ.idsp1) THEN
               Vdsi(1) = Vdsi(1) + DREAL( V(1)*densite*DCONJG(V(1)) )
            ELSEIF (IN.EQ.idsp2) THEN   
               Vdsi(2) = Vdsi(2) + DREAL( V(1)*densite*DCONJG(V(1)) )
            ELSEIF (IN.EQ.idsp3) THEN   
               Vdsi(3) = Vdsi(3) + DREAL( V(1)*densite*DCONJG(V(1)) )
            ELSEIF (IN.EQ.idsp4) THEN   
               Vdsi(4) = Vdsi(4) + DREAL( V(1)*densite*DCONJG(V(1)) )      
            ELSEIF (IN.EQ.idsp5) THEN   
               Vdsi(5) = Vdsi(5) + DREAL( V(1)*densite*DCONJG(V(1)) )
            ENDIF
            
            VGdsi(in) = VGdsi(in) + DREAL( V(1)*densite*DCONJG(V(1)) )
          
          ENDIF
        ENDDO

      ENDDO

      RETURN
      END
      

      
      REAL*8 FUNCTION CORCOS(kx,freq)
C===============================================================
C Calcul de la densite spectrale des fluctuations de pression
C sur une plaque plane : model de Corcos.
C La fonction dimensionne le spectre
C===============================================================

      include 'femak.cmn'
c      COMMON/CMDSP/dir_xy,excit,idsp1,idsp2,idsp3,idsp4,idsp5,
c     +             vit_son1,rhof1,vit_infi,del_cl,kspace 

      
      REAL*8 vit_frot,vit_conv
      REAL*8 kx,kc,freq
      REAL*8 aa1,aa2,aap,gamma
      REAL*8 csup,cinf
      REAL*8 pi

      DATA pi/3.14159D0/

c constante de Corcos:
      aa1   = 0.09D0
      aa2   = 7.D0 * aa1
      aap   = 0.766D0
      gamma = 0.389D0

      vit_frot = 0.03D0 * vit_infi
      vit_conv = 0.70D0 * vit_infi
      kc = 2*pi*freq / vit_conv
      
c calcul de la densite
      csup = aap*(1.D0+gamma)*rhof1*rhof1*vit_frot**4/(2*pi*freq)
      cinf = aa2*pi*pi*((kx-kc)**2+(aa1*kc)**2)
      CORCOS = aa1*csup/cinf
      
      RETURN
      END



      REAL*8 FUNCTION CHASE(kx,freq)
C===============================================================
C Calcul de la densite spectrale des fluctuations de pression
C sur une plaque plane : model de Chase.
C La fonction dimensionne le spectre
C===============================================================

      include 'femak.cmn'
c      COMMON/CMDSP/dir_xy,excit,idsp1,idsp2,idsp3,idsp4,idsp5,
c     +             vit_son1,rhof1,vit_infi,del_cl,kspace 

      
      REAL*8 vit_frot,vit_conv
      REAL*8 kx,freq
      REAL*8 kc,kp,c2,c3,b,h,Ct,Cm
      REAL*8 term1,term2,term3,term4,term5
      REAL*8 pi

      DATA pi/3.14159D0/

c constante de Chase:
      c2 = 1.D0/6.D0
      c3 = c2
      b  = 0.75D0
      h  = 3.D0
      Ct = 0.014D0/h
      Cm = 0.466D0/h      
      vit_frot = 0.03D0 * vit_infi
      vit_conv = 0.70D0 * vit_infi
      kp = (2*pi*freq-vit_conv*kx) / (h*vit_frot)      
      kc =  DSQRT(DABS( kx**2 - (2*pi*freq/vit_son1)**2 ))
      
c calcul de la densite
      term1 = rhof1*rhof1*vit_frot*vit_frot*vit_frot
      term2 = (kp*kp + 1/(b*b*del_cl*del_cl))**2.5D0
      term3 = c2*(kc/kx)**2 + c3*(kx/kc)**2 + 1.D0 - c2 - c3
      term4 = Ct*kx*kx*
     &     ((kp*kp + 1.D0/(b*del_cl)**2) / (kx*kx + 1.D0/(b*del_cl)**2))
      term5 = Cm * kx*kx*kx*kx / kc**2 
      
      CHASE = term1/term2 * ( term3*term4 + term5 )
      
      RETURN
      END
      
      REAL*8 FUNCTION CHASEM(kx,freq)
C===============================================================
C Calcul de la densite spectrale des fluctuations de pression
C sur une plaque plane : model de Chase du pgm Mazzoni.
C La fonction dimensionne le spectre
C===============================================================

      include 'femak.cmn'
c      COMMON/CMDSP/dir_xy,excit,idsp1,idsp2,idsp3,idsp4,idsp5,
c     +             vit_son1,rhof1,vit_infi,del_cl,kspace 

      
      REAL*8 vit_frot,vit_conv
      REAL*8 kx,freq
      REAL*8 bb,hh,aa1,aa2,aa3,epsilon,ct,cm
      REAL*8 kx2,bb2,kacou,kacou2,delta2,kplus2
      REAL*8 pm,pt
      REAL*8 pi

      DATA pi/3.14159D0/

c constante de Chasem

      bb=0.75D0
      hh=3.D0
      vit_conv = 0.70D0 * vit_infi
      vit_frot = 0.03D0 * vit_infi
      aa1=2.D0/3.D0
      aa2=1.D0/6.D0
      aa3=1.D0/6.D0
      epsilon = 0.04D0
      ct=0.0047D0
      cm=0.1553D0            
! CALCUL DU SPECTRE DE HOWE :

      kx2 = kx * kx
      bb2  = 1.D0 / bb / bb
      kacou = 2*pi*freq / vit_son1
      kacou2 = kacou * kacou
      delta2  = del_cl * del_cl

      kplus2 = ( 2*pi*freq - vit_conv * kx ) / hh / vit_frot
      kplus2 = ( kplus2 * kplus2 + kx2 ) * delta2 + bb2

      pm = DABS( kx2 - kacou2 )

      pt = ( aa1
     &     + aa2 * pm / kx2 
     &     + aa3 * kx2 / ( pm + epsilon * kacou2 ) )
     &     * ct  * kx2 * delta2 * kplus2
     &     / ( kx2 * delta2 + bb2 )

      pm = cm * kx2 * delta2 * kx2 / (pm + epsilon * kacou2)

      CHASEM  = ( pt + pm ) / ( kplus2**2.5D0 ) * del_cl * delta2
     &        * rhof1 * rhof1 * vit_frot * vit_frot * vit_frot
     
c      write(6,*) 'pt=',pt
c      write(6,*) 'pm=',pm
c      write(6,*) 'chasem=',chasem
      
      RETURN
      END      