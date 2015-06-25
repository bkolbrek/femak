      SUBROUTINE D20(VPREE,VDE)
C********************************************************************C
C D20.FOR                                                            C
C  29/5-1996    NTH                                                 C
C  Tonni F. Johansen                                                 C
C********************************************************************C
C  To form matrix D (2 and 3 dimensional elasticity)                 C
C    Input                                                           C
C     VPREE element properties                                       C
C     VPREE(1)  YOUNGS' MODULUS                                      C
C     VPREE(2)  Poisson's coefficient                                C
C    Output                                                          C
C     VDE   Matrix D (full)                                          C
C********************************************************************C
C  Programmert etter Dhatt & Touzot, s.206                           C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
      REAL*8 VPREE(1),VDE(36)


      DATA ZERO/0.D0/,UN/1.0D0/,DEUX/2.0D0/

C---------------------------------------------------
      E=VPREE(1)
      X=VPREE(2)
      C1=E*(UN-X)/((UN+X)*(UN-2.*X))
      C2=C1*X/(UN-X)
      C3=E/(DEUX*(UN+X))
      IF (NDIM.EQ.2) THEN
        IF (IAX.NE.0) THEN
C..........Axisymmetrical cas      
          VDE(1)=C1
          VDE(2)=C2
          VDE(3)=C2
          VDE(4)=ZERO
          VDE(5)=C2
          VDE(6)=C1
          VDE(7)=C2
          VDE(8)=ZERO
          VDE(9)=C2
          VDE(10)=C2
          VDE(11)=C1
          VDE(12)=ZERO
          VDE(13)=ZERO
          VDE(14)=ZERO
          VDE(15)=ZERO
          VDE(16)=C3
        ELSE
C..........2-D case      
          VDE(1)=C1
          VDE(2)=C2
          VDE(3)=ZERO
          VDE(4)=C2
          VDE(5)=C1
          VDE(6)=ZERO
          VDE(7)=ZERO
          VDE(8)=ZERO
          VDE(9)=C3
        ENDIF
      ELSEIF (NDIM.EQ.3) THEN
C..........3-D case      
          VDE(1)=C1
          VDE(2)=C2
          VDE(3)=C2
          VDE(4)=ZERO
          VDE(5)=ZERO
          VDE(6)=ZERO
          VDE(7)=C2
          VDE(8)=C1
          VDE(9)=C2
          VDE(10)=ZERO
          VDE(11)=ZERO
          VDE(12)=ZERO
          VDE(13)=C2
          VDE(14)=C2
          VDE(15)=C1
          VDE(16)=ZERO
          VDE(17)=ZERO
          VDE(18)=ZERO
          VDE(19)=ZERO
          VDE(20)=ZERO
          VDE(21)=ZERO
          VDE(22)=C3
          VDE(23)=ZERO
          VDE(24)=ZERO
          VDE(25)=ZERO
          VDE(26)=ZERO
          VDE(27)=ZERO
          VDE(28)=ZERO
          VDE(29)=C3
          VDE(30)=ZERO
          VDE(31)=ZERO
          VDE(32)=ZERO
          VDE(33)=ZERO
          VDE(34)=ZERO
          VDE(35)=ZERO
          VDE(36)=C3
      ENDIF
      
      RETURN
      END
