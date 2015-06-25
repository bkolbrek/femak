      SUBROUTINE NI02(VKPG,VNI)
C**********************************************************************C
C NI02.FOR								C
C  19/3-1992	  Kompaniet Numerisk lyd				C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  To evaluate the interpolation functions and their derivatives       C
C  D(N)/D(KSI) D(N)/D(ETA) by the gerneral PN-inverse method           C
C  For dimensions 1, 2 or 3 dimensional quadratic elements             C
C    Input                                                             C
C        VKPG    coordinates at which N is to be evaluated             C
C        IPG     number of points                                      C
C	 INEL	 number of functions N (or of nodes)	   INEL.EQ.8   C
C	 NDIM	 number of dimensions			   NDIM.EQ.2   C
C    Output                                                            C
C        VNI     funcions N and derivatives                            C
C**********************************************************************C
C  Programmert etter Dhatt & Touzot, s.197-198                         C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

c      COMMON/COOR/NDIM
c      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
c     &            IPG,ICODE,IDLE0,INEL0,IPG0

      include 'femak.cmn'    

      DIMENSION VKPG(1),VNI(1)

C........INFORMATIONS REALATED TO THE 8 NODES REFERENCE SQUARE ELEMENT

C     DIMENSION VKSI(NDIM*INEL),KEXP(NDIM*INEL),KDER(NDIM)
      DIMENSION VKSI(16),KEXP(16),KDER(2)
C     DIMENSION VPN(INEL*INEL),VP(INEL),K1(INEL)
      DIMENSION VPN(64),VP(8),K1(8)

C------------NODAL Coordinates of the reference lement
      DATA VKSI/-1.D0,-1.D0,+0.D0,-1.D0,+1.D0,-1.D0,+1.D0,+0.D0,
     1           +1.D0,+1.D0,+0.D0,+1.D0,-1.D0,+1.D0,-1.D0,+0.D0/
      DATA KEXP/0,0,1,0,0,1,2,0,1,1,0,2,2,1,1,2/,IDEGR/2/

C---------------------------------------------------
      IDEG=IDEGR

C----------EVALUATE THE PN-INVERSE MATRIX
      CALL PNINV(VKSI,KEXP,VP,K1,VPN)

C----------EVALUATE N, D(N)/D(KSI), D(N)/D(ETA) AT G.P.
      I1=1
      I2=1
      DO 10 IG=1,IPG
        KDER(1)=0
        KDER(2)=0
        CALL NI(VKPG(I1),KEXP,KDER,VP,VPN,VNI(I2))
        I2=I2+INEL
        KDER(1)=1
        CALL NI(VKPG(I1),KEXP,KDER,VP,VPN,VNI(I2))
        I2=I2+INEL
        KDER(1)=0
        KDER(2)=1
        CALL NI(VKPG(I1),KEXP,KDER,VP,VPN,VNI(I2))
        I2=I2+INEL
10	I1=I1+NDIM

      RETURN
      END
