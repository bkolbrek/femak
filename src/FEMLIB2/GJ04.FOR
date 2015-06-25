      SUBROUTINE GJ04(VPREE,GTYPE,GJ)
C**********************************************************************C
C GJ04.FOR                                                             C
C  30/5-1995    NTH                                                    C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  To find frequency dependt functio for viscous and termal losses     C
C  for  Johnson-Allard model of porous material                        C
C  used with ELEM04                                                    C
C    Input                                                             C
C      VPREE  element properties                                       C
C         VPREE(1)  2 (this type of model)
C         VPREE(2)  Density of fluid                                   C
C         VPREE(3)  Fluid's Bulk Modulus.                              C
C         VPREE(4)  Prandtl's number                                   C
C         VPREE(5)  Viscosity of fluid                                 C
C         VPREE(6)  Gamma for the fluid, type of gas (for liquid ??)   C
C         VPREE(7) Porosity                                           C
C         VPREE(8) Tortuosity                                         C
C         VPREE(9) Flow resistivity                                   C
C         VPREE(10) Alpha, parameter for viscous losses                C
C         VPREE(11) Alpha-dot, parameter for thermal losses.           C
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
        GJ=SQRT(1+(4*IMONE*VPREE(8)*VPREE(8)*VPREE(5)*VPREE(2)*W/
     1      (VPREE(9)*VPREE(9)*VPREE(10)*VPREE(10)*VPREE(7)*VPREE(7))))



      ELSEIF (GTYPE.EQ.2) THEN
C--------------------------------Thermal losses
        GJ=SQRT(1+(IMONE*VPREE(2)*W*VPREE(4)*VPREE(11)*VPREE(11)/
     1      (16*VPREE(5))))

      ENDIF
      
      
      RETURN
      END
