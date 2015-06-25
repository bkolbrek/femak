      SUBROUTINE EXSOLC(VFG,KDLNC,KNEQ)
C****************************************************************C
C     EXSOLC.FOR                                                 C
C****************************************************************C
C     2/3 1992                                                   C
C     Ulf R. Kristiansen, NTH                                    C
C****************************************************************C
C     TO EXECUTE BLOCK 'SOLC'                                    C
C     READ CONCENTRATED LOADS                                    C
C     Programmed after Dhatt og Touzot p.400                     C
C================================================================C
      IMPLICIT REAL*8(A-H,O-Z)

      COMMON/TRVL/KV(16),V(14)

      include 'femak.cmn'

      DIMENSION VFG(1),KDLNC(1),KNEQ(1)
      DATA L16/16/
C----------------------------------------------------------------
C------  READ DATA
      IF(M.GE.0)WRITE(MP,2000)
2000  FORMAT(//' CARDS OF NODAL LOADS'//)
      I0=MIN0(7,NDLN)
C10    READ(M1,1000) IG,(V(I),I=1,I0)
10    READ(MR,1000) IG,(V(I),I=1,I0)
1000  FORMAT(I5,7F10.0)
C      IF(NDLN.GT.7) READ(M1,1005) (V(I),I=8,NDLN)
      IF(NDLN.GT.7) READ(MR,1005) (V(I),I=8,NDLN)
1005  FORMAT(5X,7F10.0)
      IF(M.GE.0)WRITE(MP,2010)IG,(V(I),I=1,NDLN)
2010  FORMAT(' >>>>>',I5,7E12.5/(' >>>>>',5X,7E12.5))
      IF(IG.LE.0) GO TO 60
C20    READ(M1,1010)(KV(I),I=1,L16)
20    READ(MR,1010)(KV(I),I=1,L16)
1010  FORMAT(16I5)
      IF(M.GE.0)WRITE(MP,2020)(KV(I),I=1,L16)
2020  FORMAT(' >>>>>',16I5)
C------  DECODE NODAL DATA
      DO 50 IN=1,L16
      I1=KV(IN)
      IF(I1.GT.NNT) CALL ERREUR(61,I1,NNT,1)
      IF(I1)10,10,30
30    ID1=KDLNC(I1)+1
      ID2=KDLNC(I1+1)
      J=0
      DO 50 ID=ID1,ID2
        J=J+1
        IEQ=KNEQ(ID)
        IF (IEQ)50,50,40
40      VFG(IEQ)=VFG(IEQ)+V(J)
50    CONTINUE
      GO TO 20

C------------OUTPUT
60    IF (M.GE.1) WRITE(MP,2030)(VFG(I),I=1,NEQ)
2030  FORMAT(//'TOTAL LOAD VECTOR'/(10X,10E12.5))
      
      RETURN
      END
