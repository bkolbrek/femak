      SUBROUTINE GAUSS(IPGKED,NDIM,VKPG,VCPG,IPG)
C******************************************************************************C
C    GAUSS.FOR                                                                 C
C     28/2-1992    Tonni F. Johansen                                           C
c    Kompaniet numerisk lyd                                                    C 
C******************************************************************************C
C   To form arrays of coordinates and weights at gauss points                  C 
C   (1, 2 and 3 dimensions) (1, 2, 3 or 4g.p. per dimension)                   c
C    input                                                                     c 
C       IPGKED  number of points in ksi,eta zeta directions                    C 
C       NDIM    number of dimensions (1, 2 or 3)                               C 
C     output                                                                   C 
C        VKPG    coordinates of gauss points                                   C 
C        VCPG    weights at gauss points                                       C 
C        IPG     total numbers of gauss points                                 C 
C******************************************************************************C
C                                  programmert etter Dhatt & Touzot s.265-266  C
C******************************************************************************C

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IPGKED(1),VKPG(1),VCPG(1),G(10),P(10),INDIC(4)
      DATA INDIC/1,2,4,7/
      DATA G/0.0D0,-.577350269189626D0,.577350269189626D0,
     1     -.774596669241483D0,0.0D0,.774596669241483D0,
     2     -.861136311594050D0,-.339981043584860D0,
     3     .339981043584860D0,.861136311594050D0/
      DATA P/2.0D0,1.0D0,1.0D0,
     1     0.555555555555556D0,0.888888888888889D0,0.555555555555556D0,
     2     0.347854845137450D0,0.652145154862550D0,
     3     0.652145154862550D0,0.347854845137450D0/ 

C--------------------------------------------------------------------------
      II=IPGKED(1)
      IMIN=INDIC(II)
      IMAX=IMIN+II-1
      IF (NDIM-2) 10,20,30

C----------1 DIMENSION
10    IPG=0
      DO 15 I=IMIN,IMAX
        IPG=IPG+1
        VKPG(IPG)=G(I)
        VCPG(IPG)=P(I)
15    CONTINUE
      RETURN

C-----------2 DIMENSIONS
20    II=IPGKED(2)
      JMIN=INDIC(II)
      JMAX=JMIN+II-1
      IPG=0
      L=1
      DO 25 I=IMIN,IMAX
        DO 25 J=JMIN,JMAX
          IPG=IPG+1
          VKPG(L)=G(I)
          VKPG(L+1)=G(J)
          L=L+2
          VCPG(IPG)=P(I)*P(J)
25    CONTINUE
      RETURN
C-----------3 DIMENSIONS
30    II=IPGKED(2)
      JMIN=INDIC(II)
      JMAX=JMIN+II-1
      II=IPGKED(3)
      KMIN=INDIC(II)
      KMAX=KMIN+II-1
      IPG=0
      L=1
      DO 35 I=IMIN,IMAX
        DO 35 J=JMIN,JMAX
          DO 35 K=KMIN,KMAX  
            IPG=IPG+1
            VKPG(L)=G(I)
            VKPG(L+1)=G(J)
            VKPG(L+2)=G(K)
            L=L+3
            VCPG(IPG)=P(I)*P(J)*P(K)
35    CONTINUE
      RETURN

      END
