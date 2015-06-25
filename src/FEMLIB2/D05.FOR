      SUBROUTINE D05(VPREE,BIOTM,VDEC)
C********************************************************************C
C D05.FOR                                                            C
C  19/6-1995                                                         C
C  Tonni F. Johansen  NTH                                            C
C********************************************************************C
C  To form matrix D-for Biot type porous material                    C
C    Input                                                           C
C      VPREE  element properties                                     C
C         VPREE(1)  Density of fluid                                 C
C         VPREE(2)  Fluid's Bulk Modulus.                            C
C         VPREE(3)  Prandtl's number                                 C
C         VPREE(4)  Viscosity of fluid                               C
C         VPREE(5)  Gamma for the fluid, type of gas (for liquid ??) C
C         VPREE(6)  Density of frame.                                C
C         VPREE(7)  Shear modulus of frame, real part.               C
C         VPREE(8)  Shear modulus of frame, imaginary part.          C
C         VPREE(9)  Poisson ratio of frame, real part.               C
C         VPREE(10) Poisson ratio of frame, imaginary part.          C
C         VPREE(11) Porosity                                         C
C         VPREE(12) Tortuosity                                       C
C         VPREE(13) Flow resistivity                                 C
C         VPREE(14) Alpha, parameter for viscous losses              C
C         VPREE(15) Alpha-dot, parameter for thermal losses.         C
C      BIOTM      Biot elasticity coefficient                        C
C                                  (also named Q in som refs)        C
C    Output                                                          C
C      VDEC   Matrix D (full)                                        C
C********************************************************************C
C  Based on D02 in Dhatt & Touzot, p.206                             C
C********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8 VPREE(1)
      COMPLEX*16 BIOTM,VDEC(36)
      
      INCLUDE 'femak.cmn'
      
      COMPLEX*16 E,NY,C1,C2,C3
      COMPLEX*16 CZERO

      DATA CZERO/(0.D0,0.D0)/

C---------------------------------------------------
      NY=CMPLX(VPREE(9),VPREE(10))
      E=CMPLX(VPREE(7),VPREE(8))*2*(1+NY)


      IF (NDIM.EQ.2) THEN
        IF (IAX.NE.0) THEN
C-----------Axisymmetric - case
        
          C1=E*(1.0-NY)/((1.0+NY)*(1.0-2.0*NY))
          C2=C1*NY/(1.0-NY)
          C3=CMPLX(VPREE(7),VPREE(8))
          VDEC(1)=C1+BIOTM
          VDEC(2)=C2+BIOTM
          VDEC(3)=C2+BIOTM
c          VDEC(3)=C2
          VDEC(4)=CZERO
          VDEC(5)=C2+BIOTM
          VDEC(6)=C1+BIOTM
          VDEC(7)=C2+BIOTM
c          VDEC(7)=C2
          VDEC(8)=CZERO
          VDEC(9)=C2+BIOTM
          VDEC(10)=C2+BIOTM
          VDEC(11)=C1+BIOTM
c          VDEC(9)=C2
c          VDEC(10)=C2
c          VDEC(11)=C1
          VDEC(12)=CZERO
          VDEC(13)=CZERO
          VDEC(14)=CZERO
          VDEC(15)=CZERO
          VDEC(16)=C3
        ELSE
C-----------2D - case
          C1=E*(1.0-NY)/((1.0+NY)*(1.0-2.0*NY))
          C2=C1*NY/(1.0-NY)
          C3=CMPLX(VPREE(7),VPREE(8))
          VDEC(1)=C1+BIOTM
          VDEC(2)=C2+BIOTM
          VDEC(3)=CZERO
          VDEC(4)=C2+BIOTM
          VDEC(5)=C1+BIOTM
          VDEC(6)=CZERO
          VDEC(7)=CZERO
          VDEC(8)=CZERO
          VDEC(9)=C3
        ENDIF
      ELSEIF (NDIM.EQ.3) THEN
C-----------3D - case
        C1=E*(1.0-NY)/((1.0+NY)*(1.0-2.0*NY))
        C2=C1*NY/(1.0-NY)
        C3=CMPLX(VPREE(7),VPREE(8))
        VDEC(1)=C1+BIOTM
        VDEC(2)=C2+BIOTM
        VDEC(3)=C2+BIOTM
        VDEC(4)=CZERO
        VDEC(5)=CZERO
        VDEC(6)=CZERO
        VDEC(7)=C2+BIOTM
        VDEC(8)=C1+BIOTM
        VDEC(9)=C2+BIOTM
        VDEC(10)=CZERO
        VDEC(11)=CZERO
        VDEC(12)=CZERO
        VDEC(13)=C2+BIOTM
        VDEC(14)=C2+BIOTM
        VDEC(15)=C1+BIOTM
        VDEC(16)=CZERO
        VDEC(17)=CZERO
        VDEC(18)=CZERO
        VDEC(19)=CZERO
        VDEC(20)=CZERO
        VDEC(21)=CZERO
        VDEC(22)=C3
        VDEC(23)=CZERO
        VDEC(24)=CZERO
        VDEC(25)=CZERO
        VDEC(26)=CZERO
        VDEC(27)=CZERO
        VDEC(28)=CZERO
        VDEC(29)=C3
        VDEC(30)=CZERO
        VDEC(31)=CZERO
        VDEC(32)=CZERO
        VDEC(33)=CZERO
        VDEC(34)=CZERO
        VDEC(35)=CZERO
        VDEC(36)=C3
      ENDIF
      RETURN
      END
