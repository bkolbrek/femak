      SUBROUTINE ASIG(KLD,VDIMP,KLOCE,VCORE,VPRNE,VPREE,KNE,VKE,VKE1,
     1  VIE,VIE1,VFE,VSRE,VIGS,VIGD,VIGI,VSRC,VDLE,VRES,VNDSRC,
     2  KTYPRIG,VPRIG,VFREQ)
C****************************************************************** 
C     ASIG.FOR
C****************************************************************** 
C     2/1-1995
C     NTH
C     Tonni F. Johansen
C*******************************************************************
C     TO ASSEMBLE GLOBAL MATRIX AG , ADMITTANCE FUNCTION
C     COMPLEX VARIABLES BASED ON THE SAME IDEAS AS ASKG() AND
C     ASMG()  
C******************************************************************  
C NOTE! In this routine the elements are called with ICODE=9,
C       impedance evaluation. However, matrices are multiplied with
C       the impedance in THIS routine (not standard procedure).
C******************************************************************  
      IMPLICIT REAL*8(A-H,O-Z )

      include 'femak.cmn'

C $DECLARE
      REAL*8 VA(1)
      COMMON VA

      COMPLEX*16 VIE(1),VIE1(1),VSRE(1),VIGS(1),VIGD(1),VIGI(1),VSRC(1)
      REAL*8 VDIMP(1),VCORE(1),VPRNE(1),VPREE(1),VKE(1),VFE(1),VDLE(1)
      REAL*8 VRES(1),VFREQ(1),VKE1(1)
      INTEGER KLD(1),KLOCE(1),KNE(1)

      COMPLEX*16 VNDSRC(1),VPRIG(1)
      INTEGER KTYPRIG(1)


      INTEGER IKE,IE,I,modnr,nodnr
      COMPLEX*16 IMP
C---------------------------------------------------------------------- 
C------ REWIND IMPEDANCE INFO FILE (M6), and reset element type and side
      REWIND M6
      ITPE1=0
      IFL1=0

C------LOOP OVER THE impedance surfaces

      DO IE=1,NSIMP

C------ READ AN IMPEDANCE INFO ON FILE M6
        CALL RDCIMP(M6,KLOCE,VCORE,VPRNE,VPREE,KNE,IGIMP)
c         for debugger
C        Nifl=ifl

C--------Check whether normal or special impedance element
        IF (ITPE.LT.100) THEN

C-------Not special element:
C------ EVALUATE INTERPOLATION FUNCTIONS IF REQUIRED 
          IF(ITPE.NE.ITPE1) THEN
            ICODE=2
            CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRIG,VPRIG)
c           if (itpe.NE.11) then
c            CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VKE,VFE)
c           ELSE
c              CALL ELEMLBC(VCORE,VPRNE,VPREE,VDLE,VKE,VIE,VFE,VSRE,
c     1      VNDSRC,KTYPRSG,VPRSG)              
c           ENDIF
C.............Beregner antall plasser i element matrise (IKE)
           IF (NSYM.EQ.0) THEN
             IKE=IDLE*(IDLE+1)/2
           ELSE
             IKE=IDLE*IDLE
           ENDIF
        ENDIF
C------ EVALUATE ELEMENT VECTOR 
        ICODE=9
C        IF (NSYM.EQ.0) THEN

C          CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VKE,VKE1,VIE,VIE1,VFE,
C     1            VSRE,VNDSRC,KTYPRIG,VPRIG)
c            if (itpe.NE.11) then
c        CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VKE,VFE)
c            ELSE
c              CALL ELEMLBC(VCORE,VPRNE,VPREE,VDLE,VKE,VIE,VFE,VSRE,
c     1      VNDSRC,KTYPRSG,VPRSG)              
c            ENDIF
c          ELSE
C...NB!
C........MIDLERTIDIG !!  3/11-1993
C..NB! 
          CALL ELEMLB(VCORE,VPRNE,VPREE,VDLE,VDLEC,VKE,VKE1,VIE,VIE1,
     1            VFE,VSRE,VNDSRC,KTYPRIG,VPRIG)
C-------MULTIPLISER MED ADIMTTANSEN 

        IMP=VPRIG((IGIMP-1)*NFREQ+IFREQ)
C...NB!
C........MIDLERTIDIG !!  3/11-1993
C..NB! 
C          IF (ITPE.NE.11) THEN
C            DO I=1,IKE
C              VIE(I)=DCMPLX(0.,W*VKE(I))/IMP/VPREE(2)
C            ENDDO
C          ELSE
            DO I=1,IKE

              if (cflag.eq.1) then
                VIE(I)=DCMPLX(0.0,W)*VIE(I)/IMP/VPREE(2)
              else
                VIE(I)=DCMPLX(0.0,W)*VKE(I)/IMP/VPREE(2)
              endif




C              VIE(I)=DCMPLX(0.0,W)*VIE1(I)/IMP/VPREE(2)
            ENDDO
C          ENDIF
c
c        
c                   do i=1,ike
c                   print*,'fra asig i,vie',i,vie(i)
c                  end do


 
C------ PRINT ELEMENT VECTOR VFE
	  IF(M.GE.2) WRITE(MP,2000) IEL,(VIE(I),I=1,IKE)
2000	  FORMAT(/' VECTOR (IE) , ELEMENT:',I5/(10X,10E12.5))

C------ ASSEMBLE
        CALL ASSELC(1,0,IDLE,NSYM,KLOCE,KLD,VIE,VSRE,VIGS,VIGD,VIGI,
     1        VSRC)
C------------Finished a normal impedance element

      ELSEIF (ITPE.EQ.100) THEN
C-----------Special imp.element no.100
        modnr=ni100mod
        nodnr=ni100nod
        CALL IMP100(VA(LI100KM),VA(LI100E),VA(LI100T),VA(LI100Y),
     1   VA(LI100Q),VA(LI100DM1),VA(LI100DM2),MODnr,
     2   nodnr,VFREQ,VPREE,VDIMP,KLOCE,VCORE,VPRNE,VIGS,
     3   VIGD,VIGI,KLD,KNE)

C-----------Finished special imp.element no.100

        ENDIF

      ENDDO
c----------End loop through all impedance surfaces


      RETURN
      END
