      SUBROUTINE NI01(VKPG,VNI)
C**********************************************************************C
C NI01.FOR                                                             C
C  21/2-1992     Kompaniet Numerisk lyd                                C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  To evaluate the interpolation functions and their derivatives       C
C  D(N)/D(KSI) D(N)/D(ETA) by the gerneral PN-inverse method           C
C  For dimensions 1, 2 or 3 dimensional quadratic elements             C
C    Input                                                             C
C        VKPG    coordinates at which N is to be evaluated             C
C        IPG     number of points                                      C
C        INEL    number of functions N (or of nodes)       INEL.LE.20  C
C        NDIM    number of dimensions                      NDIM.LE.3   C
C    Output                                                            C
C        VNI     funcions N and derivatives                            C
C**********************************************************************C
C  Programmert etter Dhatt & Touzot, s.197-198                         C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

c      COMMON/COOR/NDIM
c      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
c     &            IPG,ICODE,IDLE0,INEL0,IPG0
c      COMMON/TRVL/VKSI,VPN,VP,KEXP,KDER,K1

      include 'femak.cmn'    

      DIMENSION VKPG(1),VNI(1)
      DIMENSION VKSI1(3),KEXP1(3),VKSI2(16),KEXP2(16),VKSI3(60)
      DIMENSION KEXP3(60)

C-----------INFORMATION TO DEFINE THE 3 REFERENCE ELEMENTS
C     DIMENSION VKSI(NDIM*INEL),KEXP(NDIM*INEL),KDER(NDIM)
      DIMENSION VKSI(60),KEXP(60),KDER(3)
C     DIMENSION VPN(INEL*INEL),VP(INEL)
      DIMENSION VPN(400),VP(20)
C     DIMENSION K1(INEL)
      DIMENSION K1(20)

C------------Characteristics for 1, 2 and 3 dimensional reference elements
      DATA VKSI1/-1.D0,0.D0,1.D0/,KEXP1/0,1,2/
      DATA VKSI2/-1.D0,-1.D0,+0.D0,-1.D0,+1.D0,-1.D0,+1.D0,+0.D0,
     1           +1.D0,+1.D0,+0.D0,+1.D0,-1.D0,+1.D0,-1.D0,+0.D0/
      DATA KEXP2/0,0,1,0,0,1,2,0,1,1,0,2,2,1,1,2/,IDEGR/2/
      DATA VKSI3/-1.D0,-1.D0,-1.D0,+0.D0,-1.D0,-1.D0,
     1           +1.D0,-1.D0,-1.D0,+1.D0,+0.D0,-1.D0,
     2           +1.D0,+1.D0,-1.D0,+0.D0,+1.D0,-1.D0,
     3           -1.D0,+1.D0,-1.D0,-1.D0,+0.D0,-1.D0,
     4           -1.D0,-1.D0,+0.D0,+1.D0,-1.D0,+0.D0,
     5           +1.D0,+1.D0,+0.D0,-1.D0,+1.D0,+0.D0,
     6           -1.D0,-1.D0,+1.D0,+0.D0,-1.D0,+1.D0,
     7           +1.D0,-1.D0,+1.D0,+1.D0,+0.D0,+1.D0,
     8           +1.D0,+1.D0,+1.D0,+0.D0,+1.D0,+1.D0,
     9           -1.D0,+1.D0,+1.D0,-1.D0,+0.D0,+1.D0/
      DATA KEXP3/0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,
     1   1,1,0,0,1,1,1,0,1,2,0,0,0,2,0,0,0,2,
     2   2,1,0,2,0,1,2,1,1,1,2,0,0,2,1,1,2,1,
     3   1,0,2,0,1,2,1,1,2/

C---------------------------------------------------
      IDEG=IDEGR
C---------SELECT TABLES VKSI AND KEXP ACCORDING TO NDIM
      I1=NDIM*INEL
      DO 5 I=1,I1
        IF (NDIM.EQ.1) THEN
          VKSI(I)=VKSI1(I)
          KEXP(I)=KEXP1(I)
        ELSEIF (NDIM.EQ.2) THEN
          VKSI(I)=VKSI2(I)
          KEXP(I)=KEXP2(I)
        ELSEIF (NDIM.EQ.3) THEN
          VKSI(I)=VKSI3(I)
          KEXP(I)=KEXP3(I)
        ENDIF
5     CONTINUE

C----------EVALUATE THE PN-INVERSE MATRIX
      CALL PNINV(VKSI,KEXP,VP,K1,VPN)

C----------EVALUATE N, D(N)/D(KSI), D(N)/D(ETA) AT G.P.
      I1=1
      I2=1
      DO 10 IG=1,IPG
        KDER(1)=0
        KDER(2)=0
        KDER(3)=0
        CALL NI(VKPG(I1),KEXP,KDER,VP,VPN,VNI(I2))
        I2=I2+INEL
        KDER(1)=1
        CALL NI(VKPG(I1),KEXP,KDER,VP,VPN,VNI(I2))
        I2=I2+INEL
        IF (NDIM.EQ.1) GO TO 10
        KDER(1)=0
        KDER(2)=1
        CALL NI(VKPG(I1),KEXP,KDER,VP,VPN,VNI(I2))
        I2=I2+INEL
        IF (NDIM.EQ.2) GO TO 10
        KDER(2)=0
        KDER(3)=1
        CALL NI(VKPG(I1),KEXP,KDER,VP,VPN,VNI(I2))
        I2=I2+INEL
10	I1=I1+NDIM

      RETURN
      END
