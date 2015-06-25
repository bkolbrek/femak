      SUBROUTINE NI03(VKPG,VNI)
C**********************************************************************C
C NI03.FOR                                                             C
C  23/5-1995   Kompaniet Numerisk lyd                                  C
C  Bruno Brouard                                                       C
C**********************************************************************C
C  To evaluate the interpolation functions and their derivatives       C
C  D(N)/D(KSI) D(N)/D(ETA) by the gerneral PN-inverse method           C
C  For dimensions 1, 2 or 3 dimensional quadratic elements             C
C    Input                                                             C
C        VKPG    coordinates at which N is to be evaluated             C
C        IPG     number of points                                      C
C        INEL    number of functions N (or of nodes)       INEL.EQ.16  C
C        NDIM    number of dimensions                      NDIM.EQ.3   C
C    Output                                                            C
C        VNI     functions N and derivatives                           C
C**********************************************************************C
C  The routine Ni02.for has been the model for this routine            C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

c      COMMON/COOR/NDIM
c      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
c     &            IPG,ICODE,IDLE0,INEL0,IPG0

      include 'femak.cmn'    

      DIMENSION VKPG(1),VNI(1)

C........INFORMATIONS RELATED TO THE 6 NODES REFERENCE SQUARE ELEMENT
      DIMENSION VKSI2(12),KEXP2(12),VKSI3(48),KEXP3(48)
C     DIMENSION VKSI(NDIM*INEL),KEXP(NDIM*INEL),KDER(NDIM)
      DIMENSION VKSI(48),KEXP(48),KDER(3)
C     DIMENSION VPN(INEL*INEL),VP(INEL),K1(INEL)
      DIMENSION VPN(256),VP(16),K1(16)

C------------NODAL Coordinates of the reference lement
      DATA VKSI2/-1.D0,-1.D0,-1.D0,0.D0,-1.D0,+1.D0,+1.D0,-1.D0,
     1           +1.D0,0.D0,+1.D0,+1.D0/
      DATA KEXP2/0,0,1,0,0,1,1,1,0,2,1,2/,IDEGR/2/
      DATA VKSI3/-1.D0,-1.D0,-1.D0,+0.D0,-1.D0,-1.D0,
     1           +1.D0,-1.D0,-1.D0,+1.D0,+0.D0,-1.D0,
     2           +1.D0,+1.D0,-1.D0,+0.D0,+1.D0,-1.D0,
     3           -1.D0,+1.D0,-1.D0,-1.D0,+0.D0,-1.D0,
     4           -1.D0,-1.D0,+1.D0,+0.D0,-1.D0,+1.D0,
     5           +1.D0,-1.D0,+1.D0,+1.D0,+0.D0,+1.D0,
     6           +1.D0,+1.D0,+1.D0,+0.D0,+1.D0,+1.D0,
     7           -1.D0,+1.D0,+1.D0,-1.D0,+0.D0,+1.D0/
      DATA KEXP3/0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,
     1           1,1,0,0,1,1,1,0,1,2,0,0,0,2,0,
     2           2,1,0,1,2,0,2,0,1,0,2,1,2,1,1,
     3           1,2,1/

C---------------------------------------------------
      IDEG=IDEGR
C---------SELECT TABLES VKSI AND KEXP ACCORDING TO NDIM
      I1=NDIM*INEL
      DO 5 I=1,I1
        IF (NDIM.EQ.2) THEN
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
        KDER(1)=0
        KDER(2)=1
        CALL NI(VKPG(I1),KEXP,KDER,VP,VPN,VNI(I2))
        I2=I2+INEL
        IF (NDIM.EQ.2) GO TO 10
        KDER(2)=0
        KDER(3)=1
        CALL NI(VKPG(I1),KEXP,KDER,VP,VPN,VNI(I2))
        I2=I2+INEL
10    I1=I1+NDIM

      RETURN
      END
