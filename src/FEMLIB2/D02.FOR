      SUBROUTINE D02(VPREE,VDE)
C********************************************************************C
C D02.FOR                                                            C
C  19/3-1992                                                         C
C  Tonni F. Johansen NTH                                             C
C********************************************************************C
C  To form matrix D (2 dimensional elasticity)                       C
C    Input                                                           C
C     VPREE element properties                                       C
C     VPREE(1)  YOUNGS' MODULUS                                      C
C     VPREE(2)  Poisson's coefficient                                C
C     VPREE(3)  EQ.0 Plane stresses                                  C
C               EQ.1 Plane strains                                   C
C    Output                                                          C
C     VDE   Matrix D (full)                                          C
C********************************************************************C
C  Programmed after Dhatt & Touzot, p.206                            C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
      REAL*8 VPREE(1),VDE(16)


      DATA ZERO/0.D0/,UN/1.0D0/,DEUX/2.0D0/

C---------------------------------------------------
      E=VPREE(1)
      X=VPREE(2)
      A=VPREE(3)
      C1=E*(UN-A*X)/((UN+X)*(UN-X-A*X))
      C2=C1*X/(UN-A*X)
      C3=E/(DEUX*(UN+X))
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
      
      RETURN
      END
