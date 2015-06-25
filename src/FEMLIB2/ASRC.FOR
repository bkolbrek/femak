      SUBROUTINE ASRC(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1  VIE,VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,VMGI,
     2  VSRC,VDLE,VRES,VNDSRC,
     3  KTYPRSG,VPRSG,KTYPRIG,VPRIG,KNEQ,VPRNG,VCORG,KDLNC)
C*******************************************************************C
C     ASRC.FOR                                                      C
C*******************************************************************C
C     25/1-1995                                                     C
C     Tonni F. Johansen    NTH                                      C
C                                                                   C
C*******************************************************************C
C     TO ASSEMBLE SOURCE VECTOR                                     C
C*******************************************************************C
C NOTE! Diffent types of sources may be modelled. For some of them  C
C       the vectors/matrices are multiplied by "source strengths"   C
C       in THIS routine (not standard procedure).                   C
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
      REAL*8 DNSRC
      INTEGER i1,id,indinf,indsup,in
      REAL*8 xinf, yinf, zinf, xsup, ysup, zsup

C
      COMPLEX*16 VDLEC(1)
      include 'femak.cmn'

C  $DECLARE                             
      INTEGER I,IE,IKE,DOFNR

      COMPLEX*16 CZERO      
      DATA CZERO/(0.0,0.0)/
C----------------------------------------------------------------------
C------ REWIND SOURCE INFO FILE (M7), and reset element type
      REWIND M7
      ITPE1=0
C normalisation of the coordinate of the source distribution
      IF (excit.NE.0) THEN
      
         i1=0
         in=1
         indinf = 0
         indsup = 0
         DO WHILE ( (in.LE.nnt) .AND. ((indinf.NE.1).OR.(indsup.NE.1)) )
             i1 = ndim*(in-1) +1
             id = kdlnc(in)+1
             ifl1 = kneq(id)
             IF (ifl1.EQ.idlinf) THEN
                   xinf = VCORG(i1)
                   IF (ndim.GE.2) yinf = VCORG(i1+1)
                   IF (ndim.GE.3) zinf = VCORG(i1+2)
                   indinf = 1
             ELSEIF (ifl1.EQ.idlsup) THEN      
                   xsup = VCORG(i1)
                   IF (ndim.GE.2) ysup = VCORG(i1+1)
                   IF (ndim.GE.3) zsup = VCORG(i1+2)
                   indsup = 1
             ENDIF
             in = in + 1
          END DO 
          IF (indinf.NE.1) STOP 'error message from ASRC.FOR'
          IF (indsup.NE.1) STOP 'error message from ASRC.FOR'

       ENDIF                  
C------LOOP OVER THE ELEMENTS      
      DO IE=1,NSRC
C------ READ AN SOURCE INFO ON FILE M7
        CALL RDSRCE(M7,KLOCE,VCORE,VPRNE,VPREE,KNE,VNDSRC,KTYPRSG)
        IF (KTYPRSG(IGSRC).EQ.0) THEN
C---------------Velocity source, constant source strength
          DOFNR=KNEQ(KDLNC(IFL)+1)
          VSRC(DOFNR)=CMPLX(0.0,-VPREE(1)*W)
     1                     *VPRSG((IGSRC-1)*NFREQ+IFREQ)

CMARC
C==============test TYPE 1 in 2D !!
        ELSEIF (ktyprsg(igsrc).EQ.1) THEN
        
          pi   = 3.14159D0
          dnsrc = DBLE(nsrc)
          
          DO i=1,inel
             in = kne(i)
             id = kdlnc(in)+1
             ifl1 = kneq(id)
             IF (ifl.EQ.ifl1) icore = ndim*(i-1)+1
          END DO
       
          IF (dir_xy.EQ.1) THEN
             vsrc(ifl)=(DABS(xsup-xinf)/dnsrc)*CDEXP(2*pi*DCMPLX(0,1.D0)
     .                      *kspace*( vcore(icore) - (xinf+xsup)/2.D0 ))
          ELSEIF (dir_xy.EQ.2) THEN
             vsrc(ifl)=(DABS(ysup-yinf)/dnsrc)*CDEXP(2*pi*DCMPLX(0,1.D0)
     .                    *kspace*( vcore(icore+1) - (yinf+ysup)/2.D0 ))
          ENDIF 
C==============END test TYPE 1 in 2D !!          
CMARC                      

        ELSEIF (KTYPRSG(IGSRC).EQ.3) THEN
C---------------Nodal load
          DOFNR=KNEQ(KDLNC(IFL)+IEL)
          VSRC(DOFNR)=VPRSG((IGSRC-1)*NFREQ+IFREQ)


C--------START TEST FOR NEW TYPE SOURCE FOR ELEM05 NB! 23/9-1994
         ELSEIF ((KTYPRSG(IGSRC).EQ.25).OR.(KTYPRSG(IGSRC).EQ.26)) THEN 
C........Set element array to zero
          IF (NDIM.EQ.2) THEN
            DO I=1,IDLE
              VSRE(I)=CZERO
            ENDDO
          ELSEIF (NDIM.EQ.3) THEN
            DO I=1,IDLE
              VSRE(I)=CZERO
            ENDDO
          ENDIF

          IF(ITPE.NE.ITPE1) THEN    
            ICODE=2         
            CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1          VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
          ENDIF
C------ EVALUATE ELEMENT VECTOR         
        ICODE=10
        IF (NSYM.EQ.0) THEN
            CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
          ELSE
            CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
        ENDIF
C      ENDIF



C--------END TEST FOR NEW TYPE SOURCE FOR ELEM05 NB! 23/9-1994


      ELSEIF ((KTYPRSG(IGSRC).EQ.20).OR.(KTYPRSG(IGSRC).EQ.21)) THEN 
C........Set element array to zero
        IF (NDIM.EQ.2) THEN
          DO I=1,8
            VSRE(I)=CZERO
          ENDDO
        ELSEIF (NDIM.EQ.3) THEN
          DO I=1,20
            VSRE(I)=CZERO
          ENDDO
        ENDIF

C------ EVALUATE INTERPOLATION FUNCTIONS IF REQUIRED    
        IF(ITPE.NE.ITPE1) THEN    
          ICODE=2         
            CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
C.............Beregner antall plasser i element matrise (IKE)         
          IF (NSYM.EQ.0) THEN       
            IKE=IDLE*(IDLE+1)/2           
          ELSE      
            IKE=IDLE*IDLE         
          ENDIF        
        ENDIF
C------ EVALUATE ELEMENT VECTOR         
        ICODE=9
        IF (NSYM.EQ.0) THEN
          CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
        ELSE
          CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE1,VKE,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRSG,VPRSG)
        ENDIF
C.........Discreet velocity value for each node. Each side (IFL) has to
C.........be treated separately         
      IF (NDIM.EQ.2) THEN       
      IF (NSYM.EQ.0) THEN
        IF (IFL.EQ.2) THEN        
          VSRE(1)=DCMPLX(0.,-W)*VPREE(1)*(VKE(1)*VNDSRC((IFREQ-1)*4+3)+     
     1                   VKE(22)*VNDSRC((IFREQ-1)*4+1)+     
     2                   VKE(29)*VNDSRC((IFREQ-1)*4+2))          
          VSRE(7)=DCMPLX(0.,-W)*VPREE(1)*(VKE(22)*VNDSRC((IFREQ-1)*4+3)+     
     1                  VKE(28)*VNDSRC((IFREQ-1)*4+1)+     
     2                   VKE(35)*VNDSRC((IFREQ-1)*4+2))          
          VSRE(8)=DCMPLX(0.,-W)*VPREE(1)*(VKE(29)*VNDSRC((IFREQ-1)*4+3)+     
     1                  VKE(35)*VNDSRC((IFREQ-1)*4+1)+     
     2                  VKE(36)*VNDSRC((IFREQ-1)*4+2))        
        ELSEIF (IFL.EQ.3) THEN            
          VSRE(3)=DCMPLX(0.,-W)*VPREE(1)*(VKE(6)*VNDSRC((IFREQ-1)*4+1)+     
     1                  VKE(9)*VNDSRC((IFREQ-1)*4+2)+     
     2                  VKE(13)*VNDSRC((IFREQ-1)*4+3))          
          VSRE(4)=DCMPLX(0.,-W)*VPREE(1)*(VKE(9)*VNDSRC((IFREQ-1)*4+1)+     
     1                   VKE(10)*VNDSRC((IFREQ-1)*4+2)+     
     2                    VKE(14)*VNDSRC((IFREQ-1)*4+3))          
          VSRE(5)=DCMPLX(0.,-W)*VPREE(1)*(VKE(13)*VNDSRC((IFREQ-1)*4+1)+     
     1                  VKE(14)*VNDSRC((IFREQ-1)*4+2)+     
     2                    VKE(15)*VNDSRC((IFREQ-1)*4+3))
        ELSEIF (IFL.EQ.4) THEN            
          VSRE(1)=DCMPLX(0.,-W)*VPREE(1)*(VKE(1)*VNDSRC((IFREQ-1)*4+1)+     
     1                   VKE(2)*VNDSRC((IFREQ-1)*4+2)+     
     2                     VKE(4)*VNDSRC((IFREQ-1)*4+3))           
          VSRE(2)=DCMPLX(0.,-W)*VPREE(1)*(VKE(2)*VNDSRC((IFREQ-1)*4+1)+     
     1                   VKE(3)*VNDSRC((IFREQ-1)*4+2)+     
     2                     VKE(5)*VNDSRC((IFREQ-1)*4+3))           
	    VSRE(3)=DCMPLX(0.,-W)*VPREE(1)*(VKE(4)*VNDSRC((IFREQ-1)*4+1)+     
     1                   VKE(5)*VNDSRC((IFREQ-1)*4+2)+     
     2                     VKE(6)*VNDSRC((IFREQ-1)*4+3))
        ELSEIF (IFL.EQ.5) THEN            
	    VSRE(5)=DCMPLX(0.,-W)*VPREE(1)*(VKE(15)*VNDSRC((IFREQ-1)*4+1)+     
     1                  VKE(20)*VNDSRC((IFREQ-1)*4+2)+     
     2                    VKE(26)*VNDSRC((IFREQ-1)*4+3))          
	    VSRE(6)=DCMPLX(0.,-W)*VPREE(1)*(VKE(20)*VNDSRC((IFREQ-1)*4+1)+     
     1                 VKE(21)*VNDSRC((IFREQ-1)*4+2)+     
     2                 VKE(27)*VNDSRC((IFREQ-1)*4+3))          
          VSRE(7)=DCMPLX(0.,-W)*VPREE(1)*(VKE(26)*VNDSRC((IFREQ-1)*4+1)+     
     1                  VKE(27)*VNDSRC((IFREQ-1)*4+2)+     
     2                    VKE(28)*VNDSRC((IFREQ-1)*4+3))
C..............Finished NDIM = 2
        ENDIF
        ELSE
C......note dangerous trick using VKE1 !!!          
        IF (IFL.EQ.2) THEN        
          VSRE(1)=DCMPLX(0.,-W)*VPREE(1)*(VKE(1)*VNDSRC((IFREQ-1)*4+3)+     
     1                   VKE(22)*VNDSRC((IFREQ-1)*4+1)+     
     2                   VKE(29)*VNDSRC((IFREQ-1)*4+2))          
          VSRE(7)=DCMPLX(0.,-W)*VPREE(1)*(VKE(22)*VNDSRC((IFREQ-1)*4+3)+     
     1                  VKE(28)*VNDSRC((IFREQ-1)*4+1)+     
     2                   VKE(35)*VNDSRC((IFREQ-1)*4+2))          
          VSRE(8)=DCMPLX(0.,-W)*VPREE(1)*(VKE(29)*VNDSRC((IFREQ-1)*4+3)+     
     1                  VKE(35)*VNDSRC((IFREQ-1)*4+1)+     
     2                  VKE(36)*VNDSRC((IFREQ-1)*4+2))        
        ELSEIF (IFL.EQ.3) THEN            
          VSRE(3)=DCMPLX(0.,-W)*VPREE(1)*(VKE(6)*VNDSRC((IFREQ-1)*4+1)+     
     1                  VKE(9)*VNDSRC((IFREQ-1)*4+2)+     
     2                  VKE(13)*VNDSRC((IFREQ-1)*4+3))          
          VSRE(4)=DCMPLX(0.,-W)*VPREE(1)*(VKE(9)*VNDSRC((IFREQ-1)*4+1)+     
     1                   VKE(10)*VNDSRC((IFREQ-1)*4+2)+     
     2                    VKE(14)*VNDSRC((IFREQ-1)*4+3))          
          VSRE(5)=DCMPLX(0.,-W)*VPREE(1)*(VKE(13)*VNDSRC((IFREQ-1)*4+1)+     
     1                  VKE(14)*VNDSRC((IFREQ-1)*4+2)+     
     2                    VKE(15)*VNDSRC((IFREQ-1)*4+3))
        ELSEIF (IFL.EQ.4) THEN            
          VSRE(1)=DCMPLX(0.,-W)*VPREE(1)*(VKE(1)*VNDSRC((IFREQ-1)*4+1)+     
     1                   VKE(2)*VNDSRC((IFREQ-1)*4+2)+     
     2                     VKE(4)*VNDSRC((IFREQ-1)*4+3))           
          VSRE(2)=DCMPLX(0.,-W)*VPREE(1)*(VKE(2)*VNDSRC((IFREQ-1)*4+1)+     
     1                   VKE(3)*VNDSRC((IFREQ-1)*4+2)+     
     2                     VKE(5)*VNDSRC((IFREQ-1)*4+3))           
	    VSRE(3)=DCMPLX(0.,-W)*VPREE(1)*(VKE(4)*VNDSRC((IFREQ-1)*4+1)+     
     1                   VKE(5)*VNDSRC((IFREQ-1)*4+2)+     
     2                     VKE(6)*VNDSRC((IFREQ-1)*4+3))
        ELSEIF (IFL.EQ.5) THEN            
	    VSRE(5)=DCMPLX(0.,-W)*VPREE(1)*(VKE(15)*VNDSRC((IFREQ-1)*4+1)+     
     1                  VKE(20)*VNDSRC((IFREQ-1)*4+2)+     
     2                    VKE(26)*VNDSRC((IFREQ-1)*4+3))          
	    VSRE(6)=DCMPLX(0.,-W)*VPREE(1)*(VKE(20)*VNDSRC((IFREQ-1)*4+1)+     
     1                 VKE(21)*VNDSRC((IFREQ-1)*4+2)+     
     2                 VKE(27)*VNDSRC((IFREQ-1)*4+3))          
          VSRE(7)=DCMPLX(0.,-W)*VPREE(1)*(VKE(26)*VNDSRC((IFREQ-1)*4+1)+     
     1                  VKE(27)*VNDSRC((IFREQ-1)*4+2)+     
     2                    VKE(28)*VNDSRC((IFREQ-1)*4+3))
C..............Finished NDIM = 2
        ENDIF
        ENDIF
      ELSEIF (NDIM.EQ.3) THEN
        IF (NSYM.EQ.0) THEN
        IF (IFL.EQ.2) THEN
C...........IFL = 2  (NDIM=3 / NSYM=0)
            VSRE(1)=VPREE(1)*(VKE(1)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(37)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(79)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(191)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(172)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(67)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(22)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(29)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(9)=VPREE(1)*(VKE(37)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(45)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(87)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(199)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(180)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(75)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(43)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(44)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(13)=VPREE(1)*(VKE(79)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(87)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(91)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(203)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(184)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(90)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(85)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(86)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(20)=VPREE(1)*(VKE(191)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(199)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(203)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(210)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(209)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(202)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(197)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(198)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(19)=VPREE(1)*(VKE(172)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(180)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(184)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(209)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(190)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(183)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(178)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(179)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(12)=VPREE(1)*(VKE(67)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(75)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(90)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(202)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(183)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(78)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(73)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(74)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(7)=VPREE(1)*(VKE(22)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(43)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(85)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(197)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(178)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(73)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(28)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(35)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(8)=VPREE(1)*(VKE(29)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(44)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(86)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(198)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(179)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(74)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(35)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(36)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
        ELSEIF (IFL.EQ.3) THEN
C...........IFL = 3  (NDIM=3 / NSYM=0)
          VSRE(3)=DCMPLX(0.,-W)*VPREE(1)*(VKE(6)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(9)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(13)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(58)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(139)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(123)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(108)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(48)*VNDSRC((IFREQ-1)*10+8))
          VSRE(4)=DCMPLX(0.,-W)*VPREE(1)*(VKE(9)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(10)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(14)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(59)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(140)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(124)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(109)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(49)*VNDSRC((IFREQ-1)*10+8))
        VSRE(5)=DCMPLX(0.,-W)*VPREE(1)*(VKE(13)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(14)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(15)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(60)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(141)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(125)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(110)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(50)*VNDSRC((IFREQ-1)*10+8))
        VSRE(11)=DCMPLX(0.,-W)*VPREE(1)*(VKE(58)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(59)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(60)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(66)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(147)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(131)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(116)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(65)*VNDSRC((IFREQ-1)*10+8))
       VSRE(17)=DCMPLX(0.,-W)*VPREE(1)*(VKE(139)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(140)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(141)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(147)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(153)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(152)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(151)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(146)*VNDSRC((IFREQ-1)*10+8))
       VSRE(16)=DCMPLX(0.,-W)*VPREE(1)*(VKE(123)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(124)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(125)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(131)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(152)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(136)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(135)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(130)*VNDSRC((IFREQ-1)*10+8))
       VSRE(15)=DCMPLX(0.,-W)*VPREE(1)*(VKE(108)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(109)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(110)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(116)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(151)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(135)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(120)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(115)*VNDSRC((IFREQ-1)*10+8))
        VSRE(10)=DCMPLX(0.,-W)*VPREE(1)*(VKE(48)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(49)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(50)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(65)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(146)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(130)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(115)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(55)*VNDSRC((IFREQ-1)*10+8))
        ELSEIF (IFL.EQ.4) THEN
C...........IFL = 4  (NDIM=3 / NSYM=0)
          VSRE(3)=DCMPLX(0.,-W)*VPREE(1)*(VKE(6)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(48)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(108)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(94)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(81)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(39)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(4)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(5)*VNDSRC((IFREQ-1)*10+8))
        VSRE(10)=DCMPLX(0.,-W)*VPREE(1)*(VKE(48)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(55)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(115)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(101)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(88)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(54)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(46)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(47)*VNDSRC((IFREQ-1)*10+8))
       VSRE(15)=DCMPLX(0.,-W)*VPREE(1)*(VKE(108)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(115)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(120)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(119)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(118)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(114)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(106)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(107)*VNDSRC((IFREQ-1)*10+8))
        VSRE(14)=DCMPLX(0.,-W)*VPREE(1)*(VKE(94)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(101)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(119)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(105)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(104)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(100)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(93)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(92)*VNDSRC((IFREQ-1)*10+8))
        VSRE(13)=DCMPLX(0.,-W)*VPREE(1)*(VKE(81)*VNDSRC((IFREQ-1)*10+8)+
     +                  VKE(88)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(118)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(104)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(91)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(87)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(79)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(80)*VNDSRC((IFREQ-1)*10+8))
         VSRE(9)=DCMPLX(0.,-W)*VPREE(1)*(VKE(39)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(54)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(114)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(100)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(87)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(45)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(37)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(38)*VNDSRC((IFREQ-1)*10+8))
          VSRE(1)=DCMPLX(0.,-W)*VPREE(1)*(VKE(4)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(46)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(106)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(92)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(79)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(37)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(1)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(2)*VNDSRC((IFREQ-1)*10+8))
          VSRE(2)=DCMPLX(0.,-W)*VPREE(1)*(VKE(5)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(47)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(107)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(93)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(80)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(38)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(2)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(3)*VNDSRC((IFREQ-1)*10+8))
C
        ELSEIF (IFL.EQ.5) THEN
C...........IFL = 5  (NDIM=3 / NSYM=0)
	   VSRE(7)=DCMPLX(0.,-W)*VPREE(1)*(VKE(28)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(73)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(178)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(160)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(143)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(62)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(26)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(27)*VNDSRC((IFREQ-1)*10+8))
        VSRE(12)=DCMPLX(0.,-W)*VPREE(1)*(VKE(73)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(78)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(183)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(165)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(148)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(77)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(71)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(72)*VNDSRC((IFREQ-1)*10+8))
       VSRE(19)=DCMPLX(0.,-W)*VPREE(1)*(VKE(178)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(183)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(190)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(189)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(188)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(182)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(176)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(177)*VNDSRC((IFREQ-1)*10+8))
       VSRE(18)=DCMPLX(0.,-W)*VPREE(1)*(VKE(160)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(165)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(189)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(171)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(170)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(164)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(158)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(159)*VNDSRC((IFREQ-1)*10+8))
       VSRE(17)=DCMPLX(0.,-W)*VPREE(1)*(VKE(143)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(148)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(188)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(170)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(153)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(147)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(141)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(142)*VNDSRC((IFREQ-1)*10+8))
        VSRE(11)=DCMPLX(0.,-W)*VPREE(1)*(VKE(62)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(77)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(182)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(164)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(147)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(66)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(60)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(61)*VNDSRC((IFREQ-1)*10+8))
         VSRE(5)=DCMPLX(0.,-W)*VPREE(1)*(VKE(26)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(71)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(176)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(158)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(141)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(60)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(15)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(20)*VNDSRC((IFREQ-1)*10+8))
         VSRE(6)=DCMPLX(0.,-W)*VPREE(1)*(VKE(27)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(72)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(177)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(159)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(142)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(61)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(20)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(21)*VNDSRC((IFREQ-1)*10+8))
C
	    ELSEIF (IFL.EQ.6) THEN
C...........IFL = 6  (NDIM=3 / NSYM=0)
            VSRE(1)=VPREE(1)*(VKE(1)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(29)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(22)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(16)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(11)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(7)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(4)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(2)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(8)=VPREE(1)*(VKE(29)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(36)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(35)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(34)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(33)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(32)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(31)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(30)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(7)=VPREE(1)*(VKE(22)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(35)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(28)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(27)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(26)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(25)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(24)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(23)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(6)=VPREE(1)*(VKE(16)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(34)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(27)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(21)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(20)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(19)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(18)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(17)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(5)=VPREE(1)*(VKE(11)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(33)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(26)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(20)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(15)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(14)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(13)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(12)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(4)=VPREE(1)*(VKE(7)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(32)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(25)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(19)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(14)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(10)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(9)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(8)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(3)=VPREE(1)*(VKE(4)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(31)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(24)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(18)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(13)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(9)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(6)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(5)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(2)=VPREE(1)*(VKE(2)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(30)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(23)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(17)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(12)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(8)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(5)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(3)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
C
	    ELSEIF (IFL.EQ.7) THEN
C...........IFL = 7  (NDIM=3 / NSYM=0)
	      VSRE(13)=VPREE(1)*(VKE(91)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(104)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(118)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(133)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(149)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(166)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(184)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(203)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(14)=VPREE(1)*(VKE(104)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(105)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(119)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(134)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(150)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(167)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(185)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(204)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(15)=VPREE(1)*(VKE(118)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(119)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(120)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(135)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(151)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(168)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(186)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(205)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(16)=VPREE(1)*(VKE(133)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(134)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(135)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(136)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(152)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(169)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(187)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(206)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(17)=VPREE(1)*(VKE(149)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(150)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(151)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(152)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(153)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(170)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(188)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(207)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(18)=VPREE(1)*(VKE(166)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(167)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(168)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(169)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(170)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(171)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(189)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(208)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(19)=VPREE(1)*(VKE(184)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(185)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(186)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(187)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(188)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(189)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(190)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(209)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(20)=VPREE(1)*(VKE(203)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(204)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(205)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(206)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(207)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(208)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(209)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(210)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
          ENDIF
          
      ELSE
C......note dangerous trick using VKE1 !!!          
        IF (IFL.EQ.2) THEN
C...........IFL = 2  (NDIM=3 / NSYM=0)
            VSRE(1)=VPREE(1)*(VKE(1)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(37)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(79)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(191)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(172)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(67)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(22)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(29)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(9)=VPREE(1)*(VKE(37)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(45)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(87)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(199)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(180)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(75)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(43)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(44)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(13)=VPREE(1)*(VKE(79)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(87)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(91)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(203)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(184)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(90)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(85)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(86)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(20)=VPREE(1)*(VKE(191)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(199)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(203)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(210)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(209)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(202)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(197)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(198)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(19)=VPREE(1)*(VKE(172)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(180)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(184)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(209)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(190)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(183)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(178)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(179)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(12)=VPREE(1)*(VKE(67)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(75)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(90)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(202)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(183)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(78)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(73)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(74)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(7)=VPREE(1)*(VKE(22)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(43)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(85)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(197)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(178)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(73)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(28)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(35)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(8)=VPREE(1)*(VKE(29)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(44)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(86)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(198)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(179)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(74)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(35)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(36)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
        ELSEIF (IFL.EQ.3) THEN
C...........IFL = 3  (NDIM=3 / NSYM=0)
          VSRE(3)=DCMPLX(0.,-W)*VPREE(1)*(VKE(6)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(9)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(13)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(58)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(139)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(123)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(108)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(48)*VNDSRC((IFREQ-1)*10+8))
          VSRE(4)=DCMPLX(0.,-W)*VPREE(1)*(VKE(9)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(10)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(14)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(59)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(140)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(124)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(109)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(49)*VNDSRC((IFREQ-1)*10+8))
        VSRE(5)=DCMPLX(0.,-W)*VPREE(1)*(VKE(13)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(14)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(15)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(60)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(141)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(125)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(110)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(50)*VNDSRC((IFREQ-1)*10+8))
        VSRE(11)=DCMPLX(0.,-W)*VPREE(1)*(VKE(58)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(59)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(60)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(66)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(147)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(131)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(116)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(65)*VNDSRC((IFREQ-1)*10+8))
       VSRE(17)=DCMPLX(0.,-W)*VPREE(1)*(VKE(139)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(140)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(141)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(147)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(153)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(152)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(151)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(146)*VNDSRC((IFREQ-1)*10+8))
       VSRE(16)=DCMPLX(0.,-W)*VPREE(1)*(VKE(123)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(124)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(125)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(131)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(152)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(136)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(135)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(130)*VNDSRC((IFREQ-1)*10+8))
       VSRE(15)=DCMPLX(0.,-W)*VPREE(1)*(VKE(108)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(109)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(110)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(116)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(151)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(135)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(120)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(115)*VNDSRC((IFREQ-1)*10+8))
        VSRE(10)=DCMPLX(0.,-W)*VPREE(1)*(VKE(48)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(49)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(50)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(65)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(146)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(130)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(115)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(55)*VNDSRC((IFREQ-1)*10+8))
        ELSEIF (IFL.EQ.4) THEN
C...........IFL = 4  (NDIM=3 / NSYM=0)
          VSRE(3)=DCMPLX(0.,-W)*VPREE(1)*(VKE(6)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(48)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(108)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(94)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(81)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(39)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(4)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(5)*VNDSRC((IFREQ-1)*10+8))
        VSRE(10)=DCMPLX(0.,-W)*VPREE(1)*(VKE(48)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(55)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(115)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(101)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(88)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(54)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(46)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(47)*VNDSRC((IFREQ-1)*10+8))
       VSRE(15)=DCMPLX(0.,-W)*VPREE(1)*(VKE(108)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(115)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(120)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(119)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(118)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(114)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(106)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(107)*VNDSRC((IFREQ-1)*10+8))
        VSRE(14)=DCMPLX(0.,-W)*VPREE(1)*(VKE(94)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(101)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(119)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(105)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(104)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(100)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(93)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(92)*VNDSRC((IFREQ-1)*10+8))
        VSRE(13)=DCMPLX(0.,-W)*VPREE(1)*(VKE(81)*VNDSRC((IFREQ-1)*10+8)+
     +                  VKE(88)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(118)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(104)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(91)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(87)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(79)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(80)*VNDSRC((IFREQ-1)*10+8))
         VSRE(9)=DCMPLX(0.,-W)*VPREE(1)*(VKE(39)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(54)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(114)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(100)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(87)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(45)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(37)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(38)*VNDSRC((IFREQ-1)*10+8))
          VSRE(1)=DCMPLX(0.,-W)*VPREE(1)*(VKE(4)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(46)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(106)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(92)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(79)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(37)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(1)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(2)*VNDSRC((IFREQ-1)*10+8))
          VSRE(2)=DCMPLX(0.,-W)*VPREE(1)*(VKE(5)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(47)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(107)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(93)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(80)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(38)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(2)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(3)*VNDSRC((IFREQ-1)*10+8))
C
        ELSEIF (IFL.EQ.5) THEN
C...........IFL = 5  (NDIM=3 / NSYM=0)
	   VSRE(7)=DCMPLX(0.,-W)*VPREE(1)*(VKE(28)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(73)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(178)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(160)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(143)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(62)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(26)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(27)*VNDSRC((IFREQ-1)*10+8))
        VSRE(12)=DCMPLX(0.,-W)*VPREE(1)*(VKE(73)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(78)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(183)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(165)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(148)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(77)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(71)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(72)*VNDSRC((IFREQ-1)*10+8))
       VSRE(19)=DCMPLX(0.,-W)*VPREE(1)*(VKE(178)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(183)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(190)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(189)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(188)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(182)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(176)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(177)*VNDSRC((IFREQ-1)*10+8))
       VSRE(18)=DCMPLX(0.,-W)*VPREE(1)*(VKE(160)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(165)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(189)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(171)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(170)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(164)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(158)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(159)*VNDSRC((IFREQ-1)*10+8))
       VSRE(17)=DCMPLX(0.,-W)*VPREE(1)*(VKE(143)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(148)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(188)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(170)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(153)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(147)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(141)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(142)*VNDSRC((IFREQ-1)*10+8))
        VSRE(11)=DCMPLX(0.,-W)*VPREE(1)*(VKE(62)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(77)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(182)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(164)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(147)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(66)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(60)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(61)*VNDSRC((IFREQ-1)*10+8))
         VSRE(5)=DCMPLX(0.,-W)*VPREE(1)*(VKE(26)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(71)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(176)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(158)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(141)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(60)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(15)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(20)*VNDSRC((IFREQ-1)*10+8))
         VSRE(6)=DCMPLX(0.,-W)*VPREE(1)*(VKE(27)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(72)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(177)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(159)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(142)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(61)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(20)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(21)*VNDSRC((IFREQ-1)*10+8))
C
	    ELSEIF (IFL.EQ.6) THEN
C...........IFL = 6  (NDIM=3 / NSYM=0)
            VSRE(1)=VPREE(1)*(VKE(1)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(29)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(22)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(16)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(11)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(7)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(4)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(2)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(8)=VPREE(1)*(VKE(29)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(36)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(35)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(34)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(33)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(32)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(31)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(30)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(7)=VPREE(1)*(VKE(22)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(35)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(28)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(27)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(26)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(25)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(24)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(23)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(6)=VPREE(1)*(VKE(16)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(34)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(27)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(21)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(20)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(19)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(18)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(17)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(5)=VPREE(1)*(VKE(11)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(33)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(26)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(20)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(15)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(14)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(13)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(12)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(4)=VPREE(1)*(VKE(7)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(32)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(25)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(19)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(14)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(10)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(9)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(8)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(3)=VPREE(1)*(VKE(4)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(31)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(24)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(18)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(13)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(9)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(6)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(5)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(2)=VPREE(1)*(VKE(2)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(30)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(23)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(17)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(12)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(8)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(5)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(3)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
C
	    ELSEIF (IFL.EQ.7) THEN
C...........IFL = 7  (NDIM=3 / NSYM=0)
	      VSRE(13)=VPREE(1)*(VKE(91)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(104)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(118)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(133)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(149)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(166)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(184)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(203)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(14)=VPREE(1)*(VKE(104)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(105)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(119)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(134)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(150)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(167)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(185)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(204)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(15)=VPREE(1)*(VKE(118)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(119)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(120)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(135)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(151)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(168)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(186)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(205)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(16)=VPREE(1)*(VKE(133)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(134)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(135)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(136)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(152)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(169)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(187)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(206)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(17)=VPREE(1)*(VKE(149)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(150)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(151)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(152)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(153)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(170)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(188)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(207)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(18)=VPREE(1)*(VKE(166)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(167)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(168)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(169)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(170)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(171)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(189)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(208)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
	      VSRE(19)=VPREE(1)*(VKE(184)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(185)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(186)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(187)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(188)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(189)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(190)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(209)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
            VSRE(20)=VPREE(1)*(VKE(203)*VNDSRC((IFREQ-1)*10+1)+
     +                  VKE(204)*VNDSRC((IFREQ-1)*10+2)+
     +                  VKE(205)*VNDSRC((IFREQ-1)*10+3)+
     +                  VKE(206)*VNDSRC((IFREQ-1)*10+4)+
     +                  VKE(207)*VNDSRC((IFREQ-1)*10+5)+
     +                  VKE(208)*VNDSRC((IFREQ-1)*10+6)+
     +                  VKE(209)*VNDSRC((IFREQ-1)*10+7)+
     +                  VKE(210)*VNDSRC((IFREQ-1)*10+8))*DCMPLX(0.,-W)
          ENDIF
          
        ENDIF      
C
C..............Finished NDIM = 3
C
        ENDIF



C..............End IF (KTYPRSG.EQ.20)


       ELSEIF (KTYPRSG(IGSRC).EQ.50) THEN 
         CALL SRCTYPE50(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,
     1  VKE1,VIE,VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VSRC,VDLE,VRES,VNDSRC,
     2  KTYPRSG,VPRSG,KTYPRIG,VPRIG,KNEQ,VPRNG,VCORG,KDLNC,IE)
C..............End IF (KTYPRSG.EQ.50)
       
       ELSEIF (KTYPRSG(IGSRC).EQ.55) THEN 
         CALL SRCTYPE55(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,
     1  VKE1,VIE,VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VKGS,VKGD,VKGI,VMGS,VMGD,
     2  VMGI,VSRC,VDLE,VRES,VNDSRC,
     3  KTYPRSG,VPRSG,KTYPRIG,VPRIG,KNEQ,VPRNG,VCORG,KDLNC,IE)
C..............End IF (KTYPRSG.EQ.55)
       
      ENDIF
C..............End IF (KTYPRSG.EQ...)
C------ PRINT ELEMENT VECTOR VFEC       
      IF(M.GE.1) WRITE(MP,2000) IEL,(VSRE(I),I=1,IDLE)
2000     FORMAT(/' VECTOR (VSRE) , ELEMENT:',I5/(10X,10E12.5))

C------ ASSEMBLE IF NOT VELOCITY SOURCE
cmarc type 1
      IF ((KTYPRSG(IGSRC).NE.0).AND.(ktyprsg(igsrc).NE.1).AND.
     .    (KTYPRSG(IGSRC).NE.50).AND.(KTYPRSG(IGSRC).NE.55)) THEN
        CALL ASSELC(0,1,IDLE,NSYM,KLOCE,KLD,VIE,VSRE,VIGS,VIGD,VIGI,
     1            VSRC)        
      ENDIF
cmarc type 1 fin
C      IF ((KTYPRSG(IGSRC).NE.0).AND.(KTYPRSG(IGSRC).NE.50).AND.
C     1    (KTYPRSG(IGSRC).NE.55)) THEN
C        CALL ASSELC(0,1,IDLE,NSYM,KLOCE,KLD,VIE,VSRE,VIGS,VIGD,VIGI,
C     1            VSRC)        
C      ENDIF

C-------------Finished an element with element      
      ENDDO       

      RETURN
      END

                               