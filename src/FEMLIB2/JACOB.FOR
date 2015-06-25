      SUBROUTINE JACOB(VNI,VCORE,NDIM,INEL,VJ,VJ1,DETJ)
C******************************************************************************C
C    JACOB.FOR                                                                 C
C     28/2-1992    Tonni F. Johansen                                           C
c    Kompaniet numerisk lyd						       C
C   To evaluate the jacobian matrix, its determinant and its inverse           C
C   (1, 2 and 3 dimensions)                                                    C
C    input								       C
C	VNI    derivatives of interpolation funciont W.R.T.		       C
C       VCORE   element nodal coordinates                                      C 
C       NDIM    number of dimensions (1, 2 or 3)                               C 
C       INEL    number of nodes per element                                    C
C     output                                                                   C 
C	 VJ	 Jacobian matrix					       C
C	 VJ1	 Inverse of jacobian matrix				       C
C	 DETJ	 Determinant of jacobian matrix				       C
C                     programmert etter Dhatt & Touzot s.63                    C
C******************************************************************************C

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION VNI(INEL,1),VCORE(NDIM,1),VJ(1),VJ1(1)
      DATA ZERO/0.D0/,UN/1.D0/

C-------------------------------------------------------
C--------------FORM THE JACOBIAN MATRIX
      J=1
      DO 20 JJ=1,NDIM
        DO 20 II=1,NDIM
          C=ZERO
          DO 10 IJ=1,INEL
            C=C+VNI(IJ,II)*VCORE(JJ,IJ)
10        CONTINUE
          VJ(J)=C
          J=J+1
20    CONTINUE

      GO TO (40,50,60),NDIM

C----------------1 DIMENSIONAL INVERSION
40    DETJ=VJ(1)
      IF (DETJ.LT.ZERO) THEN
	DETJ=-DETJ
      ENDIF
      IF (DETJ.EQ.ZERO) RETURN
      VJ1(1)=UN/DETJ
      RETURN

C----------------2 DIMENSIONAL INVERSION
50    DETJ=VJ(1)*VJ(4)-VJ(2)*VJ(3)
      IF (DETJ.LT.ZERO) THEN
	DETJ=-DETJ
      ENDIF
      IF (DETJ.EQ.ZERO) RETURN
      VJ1(1)=VJ(4)/DETJ
      VJ1(2)=-VJ(2)/DETJ
      VJ1(3)=-VJ(3)/DETJ
      VJ1(4)=VJ(1)/DETJ
      RETURN

C----------------3 DIMENSIONAL INVERSION
60    DETJ=VJ(1)*(VJ(5)*VJ(9)-VJ(8)*VJ(6))+
     1     VJ(4)*(VJ(8)*VJ(3)-VJ(2)*VJ(9))+
     2     VJ(7)*(VJ(2)*VJ(6)-VJ(5)*VJ(3)) 
      IF (DETJ.LT.ZERO) THEN
	DETJ=-DETJ
      ENDIF
      IF (DETJ.EQ.ZERO) RETURN
      VJ1(1)=(VJ(5)*VJ(9)-VJ(6)*VJ(8))/DETJ
      VJ1(2)=(VJ(3)*VJ(8)-VJ(2)*VJ(9))/DETJ
      VJ1(3)=(VJ(2)*VJ(6)-VJ(3)*VJ(5))/DETJ
      VJ1(4)=(VJ(7)*VJ(6)-VJ(4)*VJ(9))/DETJ
      VJ1(5)=(VJ(1)*VJ(9)-VJ(7)*VJ(3))/DETJ
      VJ1(6)=(VJ(4)*VJ(3)-VJ(6)*VJ(1))/DETJ
      VJ1(7)=(VJ(4)*VJ(8)-VJ(7)*VJ(5))/DETJ
      VJ1(8)=(VJ(2)*VJ(7)-VJ(8)*VJ(1))/DETJ
      VJ1(9)=(VJ(1)*VJ(5)-VJ(4)*VJ(2))/DETJ
 
      RETURN

      END
