      SUBROUTINE COOR23D(VCORE,VCORE1)
C**********************************************************************C
C COOR23D.FOR                                                          C
C  3/10-1994  Kompaniet Numerisk lyd                                   C
C  Marc Dhainaut                                                       C
C**********************************************************************C
C  8 nodes quadratic element for 3 dimensional coordinates             C
C  Transform a 3D plate to 2D, using the coordinate system of the plateC
C                                                                      C
C  INPUT                                                               C
C     - VCORE : 3D coordinates of the 8 nodes                          C
C  OUTPUT                                                              C
C     - VCORE1: 2D coordinates of the 8 nodes                          C 
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'femak.cmn'
c      COMMON/COOR/NDIM
c      COMMON/RGDT/IEL,ITPE,ITPE1,IGRE,IDLE,ICE,IPRNE,IPREE,INEL,IDEG,
c     &            IPG,ICODE,IDLE0,INEL0,IPG0

c      COMMON/ES/M,MR,MP

      DIMENSION VCORE(1),VCORE1(1)

C-----------CHARACTERISTIC DIMENSIONS OF THE ELEMENT
C     DIMENSION VCORE1(NDIM*INEL),VECT(NDIM*NDIM)
      DIMENSION VECT(9),VECTINV(9)
C
      REAL*8 C,C0      
C---------------------------------------------------
      DATA ZERO/0.0D0/
      DATA EPS/1.D-6/

      SQRT(X)=DSQRT(X)
      
C--------initialize arrays
      DO I=1,9
         VECT(I) = ZERO
         VECTINV(I) = ZERO
      ENDDO
      DO I=1,3*INEL
         VCORE1(1) = ZERO
      ENDDO

               
      GO TO (100,200) NDIM

C*****************************************************
C************ 2D problem into 1D one *****************
C*****************************************************
100   CONTINUE    
C--------form the vector formed by nodes 1-3
      VECT(1)=VCORE(5)-VCORE(1)
      VECT(2)=VCORE(6)-VCORE(2)      
C--------normalized the vector 1-3 (1rst vector of new base)      
      VECT1 = SQRT(VECT(1)**2 + VECT(2)**2)
      VECT(1) = VECT(1)/VECT1
      VECT(2) = VECT(2)/VECT1
      
C--------second vector of the new base (orthogonal to the first)
      VECT(3) = -VECT(2)
      VECT(4) =  VECT(1)  
        
C--------inversion of matrix VECT
C----------------vectinv=(vect)t
      VECTINV(1)=VECT(1)
      VECTINV(2)=VECT(3)
      VECTINV(3)=VECT(2)
      VECTINV(4)=VECT(4) 
      
C--------origine of the new base at node 1
C--------writing in the new coordinate system
      DO JJ=1,INEL
         DO II=1,2
            C =ZERO
            C0=ZERO
            DO IJ=1,2
               C = C + VECTINV(II+2*IJ-2)*VCORE(IJ+2*JJ-2)
               C0= C0+ VECTINV(II+2*IJ-2)*VCORE(IJ)
            ENDDO
            VCORE1(II+2*JJ-2) = -C0+C
         ENDDO
      ENDDO 
      
C--------value of VCORE1
C      DO I=1,6
C         WRITE(MP,*) 'VCORE1=',VCORE1(I)      
C      ENDDO
C--------re-arrange VCORE1 to be compatible with ELEM08
      VCORE1(1)=VCORE1(1)
      VCORE1(2)=VCORE1(3)
      VCORE1(3)=VCORE1(5)

      RETURN
      
C*****************************************      
C*********** 3D problem into 2D one ******
C************ (8 or 9 node element) ******
C*****************************************
200   CONTINUE
C--------form the two vectors formed by nodes 1-3, 1-7
      I1=6
      DO JJ=1,2
         DO I=1,3
            VECT(I+3*JJ-3) = -VCORE(I)+VCORE(I+I1)
         ENDDO
         I1=I1+12
      ENDDO
C--------normalized the vector 1-3 (1rst vector of new base)      
      VECT1 = SQRT(VECT(1)**2 + VECT(2)**2 + VECT(3)**2)
      DO I=1,3
         VECT(I) = VECT(I)/VECT1
      ENDDO
      
C--------vectorial product between vector 1-3 and 1-7
      VECT(7) = VECT(2)*VECT(6)-VECT(3)*VECT(5)
      VECT(8) = VECT(3)*VECT(4)-VECT(1)*VECT(6)
      VECT(9) = VECT(1)*VECT(5)-VECT(2)*VECT(4)         
C--------normalized the vector obtained (3rd vector of new base)      
      VECT3 = SQRT(VECT(7)**2 + VECT(8)**2 + VECT(9)**2)
      DO I=1,3
         VECT(I+6) = VECT(I+6)/VECT3
      ENDDO
            
C--------vectorial product between 3rd and 1rst vectors of new base
      VECT(4) =  VECT(3)*VECT(8)-VECT(2)*VECT(9)
      VECT(5) =  VECT(1)*VECT(9)-VECT(3)*VECT(7)
      VECT(6) =  VECT(2)*VECT(7)-VECT(1)*VECT(8)       
C--------normalized the vector obtained (2rd vector of new base)      
      VECT2 = SQRT(VECT(4)**2 + VECT(5)**2 + VECT(6)**2)
      DO I=1,3
         VECT(I+3) = VECT(I+3)/VECT2
      ENDDO
                
C--------inversion of matrix VECT
C----------------3 DIMENSIONAL INVERSION
60    DETV=VECT(1)*(VECT(5)*VECT(9)-VECT(8)*VECT(6))+
     1     VECT(4)*(VECT(8)*VECT(3)-VECT(2)*VECT(9))+
     2     VECT(7)*(VECT(2)*VECT(6)-VECT(5)*VECT(3)) 
C      IF (DETV.LT.ZERO) THEN
C          DETV=-DETV
C      ENDIF
      IF (DETV.EQ.ZERO) RETURN
      VECTINV(1)=(VECT(5)*VECT(9)-VECT(6)*VECT(8))/DETV
      VECTINV(2)=(VECT(3)*VECT(8)-VECT(2)*VECT(9))/DETV
      VECTINV(3)=(VECT(2)*VECT(6)-VECT(3)*VECT(5))/DETV
      VECTINV(4)=(VECT(7)*VECT(6)-VECT(4)*VECT(9))/DETV
      VECTINV(5)=(VECT(1)*VECT(9)-VECT(7)*VECT(3))/DETV
      VECTINV(6)=(VECT(4)*VECT(3)-VECT(6)*VECT(1))/DETV
      VECTINV(7)=(VECT(4)*VECT(8)-VECT(7)*VECT(5))/DETV
      VECTINV(8)=(VECT(2)*VECT(7)-VECT(8)*VECT(1))/DETV
      VECTINV(9)=(VECT(1)*VECT(5)-VECT(4)*VECT(2))/DETV

C--------origine of the new base at node 1
C--------writing in the new coordinate system
      DO JJ=1,INEL
         DO II=1,3
            C =ZERO
            C0=ZERO
            DO IJ=1,3
               C = C + VECTINV(II+3*IJ-3)*VCORE(IJ+3*JJ-3)
               C0= C0+ VECTINV(II+3*IJ-3)*VCORE(IJ)
            ENDDO
            VCORE1(II+3*JJ-3) = -C0+C
         ENDDO
      ENDDO
C--------value of VCORE1
C      DO I=1,3*INEL
C         WRITE(MP,*) 'VCORE1=',VCORE1(I)      
C      ENDDO
C--------re-arrange VCORE1 to be compatible with ELEM08
      II=1
      JJ=1
      DO I=1,INEL
          VCORE1(II)   = VCORE1(JJ)
          VCORE1(II+1) = VCORE1(JJ+1)
          II=II+2
          JJ=JJ+3
      ENDDO

      RETURN
      END
