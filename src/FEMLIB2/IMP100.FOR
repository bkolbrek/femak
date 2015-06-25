      SUBROUTINE IMP100(VI100KM,VI100E,VI100T,VI100Y,VI100Q,VI100DM1,
     1   VI100DM2,MODNR,NODNR,VFREQ,VPREE,VDIMP,KLOCE,VCORE,VPRNE,
     2   VIGS,VIGD,VIGI,KLD,KNE)
C********************************************************************C
C     IMP100.FOR                                                     C
C     Tonni F. Johansen , NTH                                        C
C     5/6/92                                                         C
C                                                                    C
C********************************************************************C
C     To assemble "admittance matrix" for a circular hole in an      C
C     infinite baffle. As described by Kagawa et.al. J.S.&V. vol.69  C
C     1980 p.207-228.  Assemble into global matrix in this routine.  C
C********************************************************************C

      INTEGER MODNR,NODNR
      COMPLEX*16 VI100KM(MODNR),VI100E(MODNR,MODNR),VI100T(MODNR,MODNR)
      COMPLEX*16 VI100Y(NODNR,NODNR),VI100Q(NODNR,MODNR)
      COMPLEX*16 VI100DM1(NODNR,NODNR),VI100DM2(NODNR,NODNR)
      COMPLEX*16 VIGS(1),VIGD(1),VIGI(1)
      INTEGER KLD(1),KNE(1),KLOCE(1)
      REAL*8 VFREQ(1),VPREE(1),VCORE(1),VDIMP(1),VPRNE(1)

      include 'femak.cmn'

C      INTEGER NDIM,NNT
C      COMMON/COOR/NDIM,NNT

C      INTEGER NGIMP,NGSRC,NFREQ,NSIMP,NSRC
C      COMMON/AKUS/NGIMP,NGSRC,NFREQ,NSIMP,NSRC

C      INTEGER IFL,IGIMP,IFREQ,IGSRC,IFL1,NI100NOD,NI100MOD,NI100FL
C      REAL*8 I100RAD
C      COMMON/IMPEL/IFL,IGIMP,IFREQ,IGSRC,IFL1,NI100NOD,NI100MOD,
C     1           NI100FL,I100RAD

C      INTEGER IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
C     &       IPG,ICODE,IDLE0,INEL0,IPG0
C      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
C     &       IPG,ICODE,IDLE0,INEL0,IPG0


C      INTEGER M,MR,MP,M1
C      COMMON/ES/M,MR,MP,M1

      INTEGER KD(30)
      COMPLEX*16 DET
      INTEGER I,J
      REAL*8 K

      REAL*8 GAMMA(10),PI
C...........Zero in first order bessel function Gamma
      DATA GAMMA/0.0,3.8317059702,7.0155866698,10.1734681351,
     1       13.3236919363,16.4706300509,19.6158585105,
     2       22.7600843806,25.9036720876,29.0468285349/
      DATA PI/3.141592654/


      W=2*PI*VFREQ(IFREQ)
      K=W/VPREE(2)



C..........Find modal wavenumbers
      DO I=1,NI100MOD
        IF (K.GE.GAMMA(I)/I100RAD) THEN
          VI100KM(I)=DCMPLX(SQRT(K*K-GAMMA(I)*GAMMA(I)/I100RAD/I100RAD)
     1        ,0.D0)
        ELSE
          VI100KM(I)=DCMPLX(0.D0,-SQRT(-K*K+GAMMA(I)*GAMMA(I)
     1            /I100RAD/I100RAD))
        ENDIF
      ENDDO

c      write(*,*)'imp100 utskrift1'

C...........Assemble matrix E
      CALL KAGINT(MODnr,VI100E,K,VI100KM)
      IF (M.GE.2) THEN
       write(mp,2002)
        DO I=1,NI100MOD
            write(mp,2000)i,(vi100e(i,j),j=1,ni100mod)
        ENDDO      
      ENDIF
2000  format(i5,20e12.5)
2002  format('e-matrix:')

      DO I=1,NI100MOD
        DO J=1,NI100MOD
          IF (J.NE.I) THEN
            VI100DM1(J,I)=VI100E(J,I)
          ELSE
            VI100DM1(J,I)=VI100E(J,I)-DCMPLX(1.,0.)
            VI100E(J,I)=VI100E(J,I)+DCMPLX(1.,0.)
          ENDIF
        ENDDO
      ENDDO
c      write(*,*)'imp100 utskrift 2'


      CALL INVERSC(VI100DM1,NI100MOD,NI100NOD,KD,DET)
c      write(*,*)'imp100 utskrift 3'
      IF (M.GE.2) THEN
       write(mp,2004)
        DO I=1,NI100MOD
            write(mp,2006)i,(vi100dm1(i,j),j=1,ni100mod)
        ENDDO      
      ENDIF
2006  format(i5,20e12.5)
2004  format('invers of e-matrix + 1:')

C.......Assemble T-matrix
      CALL MULTABC(VI100DM1,NI100MOD,NI100MOD,NI100NOD,NI100NOD,VI100E,
     1     NI100MOD,NI100MOD,NI100MOD,NI100MOD,
     1     VI100T,NI100MOD,NI100MOD,0)

      IF (M.GE.2) THEN
       write(mp,2008)
        DO I=1,NI100MOD
            write(mp,2010)i,(vi100t(i,j),j=1,ni100mod)
        ENDDO      
      ENDIF
2010  format(i5,20e12.5)
2008  format('t-matrix:')

c      write(*,*)'imp100 utskrift 4'

C..........Assemble matrix Q
      CALL KAGSRC(KLOCE,VCORE,VPRNE,VPREE,KNE,VI100Q,NI100MOD,NI100NOD)
c      write(*,*)'imp100 utskrift 5'
      IF (M.GE.2) THEN
       write(mp,2012)
        DO I=1,NI100NOD
            write(mp,2014)i,(vi100q(i,j),j=1,ni100mod)
        ENDDO      
      ENDIF
2014  format(i5,20e12.5)
2012  format('q-matrix:')

C.........Assemble admittance matrix

C.........VI100DM1()=INV(T+I)
      DO I=1,NI100MOD
        DO J=1,NI100MOD
          IF (J.NE.I) THEN
            VI100DM1(I,J)=VI100T(I,J)
          ELSE
            VI100DM1(I,J)=VI100T(I,J)+DCMPLX(1.,0.)
          ENDIF
        ENDDO
      ENDDO

      CALL INVERSC(VI100DM1,NI100MOD,NI100NOD,KD,DET)
c      write(*,*)'imp100 utskrift 6'

C...........VI100DM2()=INV(T+I)*Q'
      CALL MULTABC(VI100DM1,NI100MOD,NI100MOD,NI100NOD,NI100NOD,VI100Q,
     1     NI100NOD,NI100MOD,NI100NOD,NI100MOD,VI100DM2,
     2     NI100NOD,NI100NOD,2)
c      write(*,*)'imp100 utskrift 7'
      IF (M.GE.2) THEN
       write(mp,2016)
        DO I=1,NI100nOD
            write(mp,2018)i,(vi100dm2(i,j),j=1,ni100nod)
        ENDDO      
      ENDIF
2018  format(i5,70e12.5)
2016  format('inv(t+i)*q-matrix:')

C.........VI100DM1()=(T-I)
      DO I=1,NI100MOD
        DO J=1,NI100MOD
          IF (J.NE.I) THEN
            VI100DM1(I,J)=VI100T(I,J)
          ELSE
            VI100DM1(I,J)=VI100T(I,J)-DCMPLX(1.,0.)
          ENDIF
        ENDDO
      ENDDO

C.........VI100Y()=(T-I)*INV(T+I)*Q'  (Intermediate use of VI100Y)
      CALL MULTABC(VI100DM1,NI100MOD,NI100MOD,NI100NOD,NI100NOD,
     1    VI100DM2,NI100MOD,NI100NOD,NI100NOD,NI100NOD,
     2    VI100Y,NI100NOD,NI100NOD,0)
c      write(*,*)'imp100 utskrift 8'


C.........VI100DM1()=K*(T-I)*INV(T+I)*Q'
      DO I=1,NI100MOD
        DO J=1,NI100NOD
          VI100DM1(I,J)=VI100Y(I,J)*VI100KM(I)
        ENDDO
      ENDDO

C.........VI100Y()=1/(RO*W)*Q*K*(T-I)*INV(T+I)*Q'
      CALL MULTABC(VI100Q,NI100NOD,NI100MOD,NI100NOD,NI100MOD,VI100DM1,
     1     NI100MOD,NI100NOD,NI100NOD,NI100NOD,VI100Y,NI100NOD,
     2     NI100NOD,0)
c      write(*,*)'imp100 utskrift 9'
      DO I=1,NI100NOD
        DO J=1,NI100NOD
c	  VI100Y(I,J)=VI100Y(I,J)/W
c	  VI100Y(I,J)=dcmplx(0,-1)*VI100Y(I,J)
          VI100Y(I,J)=dcmplx(0.0,1.0)*VI100Y(I,J)
        ENDDO
      ENDDO

      IF (M.GE.2) THEN
       write(mp,2020)
        DO I=1,NI100nOD
            write(mp,2022)i,(vi100y(i,j),j=1,ni100nod)
        ENDDO      
      ENDIF
2022  format(i5,70e12.5)
2020  format('y-matrix:')



        


C........Assemble into global impedance matrices
C........Diagonal terms
      DO I=1,NI100NOD
        VIGD(KNE(I))=VIGD(KNE(I))+VI100Y(I,I)
      ENDDO

      DO I=2,NI100NOD
        DO J=1,I-1
          IF (KNE(I).GT.KNE(J)) THEN
C.............upper half of matrix
            VIGS(KLD(KNE(I)+1)-(KNE(I)-KNE(J)))=
     1           VIGS(KLD(KNE(I)+1)-(KNE(I)-KNE(J)))+VI100Y(I,J)
          ENDIF
        ENDDO
      ENDDO

      IF (NSYM.NE.0) THEN
        DO I=2,NI100NOD
          DO J=1,I-1
            IF (KNE(I).GT.KNE(J)) THEN
C.............lower half of matrix
              VIGI(KLD(KNE(I)+1)-(KNE(I)-KNE(J)))=
     1           VIGI(KLD(KNE(I)+1)-(KNE(I)-KNE(J)))+VI100Y(I,J)
            ENDIF
          ENDDO
        ENDDO
      ENDIF


      RETURN
      END 
