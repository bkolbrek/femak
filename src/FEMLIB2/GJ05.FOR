      SUBROUTINE GJ05(VPREE,GTYPE,GJ)
C**********************************************************************C
C GJ05.FOR                                                             C
C  21/9-1994    NTH                                                    C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  To find frequency dependt functio for viscous and termal losses     C
C  for  Biot type porous material                                      C
C    Input                                                             C
C      VPREE  element properties                                       C
C         VPREE(1)  Density of fluid                                   C
C         VPREE(2)  Fluid's Bulk Modulus.                              C
C         VPREE(3)  Prandtl's number                                   C
C         VPREE(4)  Viscosity of fluid                                 C
C         VPREE(5)  Gamma for the fluid, type of gas (for liquid ??)   C
C         VPREE(6)  Density of frame.                                  C
C         VPREE(7)  Shear modulus of frame, real part.                 C
C         VPREE(8)  Shear modulus of frame, imaginary part.            C
C         VPREE(9)  Poisson ratio of frame, real part.                 C
C         VPREE(10) Poisson ratio of frame, imaginary part.            C
C         VPREE(11) Porosity                                           C
C         VPREE(12) Tortuosity                                         C
C         VPREE(13) Flow resistivity                                   C
C         VPREE(14) Alpha, parameter for viscous losses                C
C         VPREE(15) Alpha-dot, parameter for thermal losses.           C
C      GTYPE      Flag identifying the type of Gj to be assembled      C
C         GTYPE.EQ.1  Gj for viscous losses                            C
C         GTYPE.EQ.2  Gj for thermal losses  (Gj')                     C
C    Output                                                            C
C      Gj   Complex number                                             C
C**********************************************************************C
c      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
c  $DECLARE

      REAL*8 VPREE(1)
      COMPLEX*16 GJ
      INTEGER GTYPE        
      
      COMPLEX*16 IMONE
      DATA IMONE/(0.0,1.0)/      
      
C---------------------------------------------------
      IF (GTYPE.EQ.1) THEN
C--------------------------------Viscous losses
        GJ=SQRT(1+(4*IMONE*VPREE(12)*VPREE(12)*VPREE(4)*VPREE(1)*
     1      W/(VPREE(13)*VPREE(13)*VPREE(14)*VPREE(14)*
     2      VPREE(11)*VPREE(11))))
      ELSEIF (GTYPE.EQ.2) THEN
C--------------------------------Thermal losses
        GJ=SQRT(1+(IMONE*VPREE(1)*W*VPREE(3)*VPREE(15)*VPREE(15)/
     1      (16*VPREE(4))))

      ENDIF
      
      
      RETURN
      END
