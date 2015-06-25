      SUBROUTINE EXCOND(VCORG,KDLNC,KNEQ,VDIMP)
C**********************************************************************C
C EXCOND.FOR                                                           C
C  1/2-1992   Kompaniet Numerisk lyd                                   C
C  Tonni F. Johansen                                                   C
C**********************************************************************C
C  To execute  Block 'COND'                                            C
C  Read boundary conditions and generate table (NEQ)                   C
C  Programmed after Dhatt & Touzot, p.389-390                          C
C**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

c      COMMON/COOR/NDIM,NNT,NDLN,NDLT,FAC(3)
c      COMMON/COND/NCLT,NCLZ,NCLNZ
c      COMMON/RESO/NEQ
c      COMMON/ES/M,MR,MP,M1
      COMMON/TRVL/KV(16),V(10),H(20),ICOD(10)
      include 'femak.cmn'


      DIMENSION VCORG(1),KDLNC(1),KNEQ(1),VDIMP(1)
      DATA L7/7/,L8/8/,L16/16/, X1/0.D0/,X2/0.D0/,X3/0.0D0/ZERO/0.D0/

C------------------------------------------------------------------------
C---------CUMULATIVE TABLE KDLNC
      DO 10 IN=1,NNT
	KDLNC(IN+1)=KDLNC(IN)+KDLNC(IN+1)
10    CONTINUE
      I1=NNT+1
      IF (M.GE.2) WRITE(MP,2000)(KDLNC(IN),IN=1,I1)
2000  FORMAT(//'NUMBER OF D.O.F. PRECEDING EACH NODE    (DLNC)'/
     1	  (1X,10I10))

C-----------INITIALIZE
      NCLT=0
      NCLNZ=0
      NCLZ=0
      IF (M.GE.0) WRITE(MP,2010)
2010  FORMAT(// 'BOUNDARY CONDITION CARDS '/)

C----------READ A B.C. GROUP CARD : 10 CODES + PRESCRIBED VALUE
20    READ(MR,1000)ICOD,(V(I),I=1,L7)
1000  FORMAT(10I1,7F10.0)
      IF (M.GE.0) WRITE(MP,2020)ICOD,(V(I),I=1,L7)
2020  FORMAT('>>>>>>',10I1,7E12.5)

C----------CHECK FOR A BLANK CARD
      J=0
      DO 30 I=1,10
	J=J+ICOD(I)
30    CONTINUE
      IF (J.EQ.0) GO TO 110

C---------READ ADDITIONAL CARD IF REQUIRED
      I2=0
      DO 40 ID=1,NDLN
	IF (ICOD(ID).LT.2) GO TO 40
	I2=I2+1
	IF (I2.NE.L8) GO TO 40
	READ(MR,1010)(V(I),I=L8,NDLN)
1010	FORMAT(10X,7F10.0)
	IF (M.GE.0) WRITE(MP,2030)(V(I),I=L8,NDLN)
2030	FORMAT('>>>>>>>>>',10X,7E12.5)
40    CONTINUE

C---------READ CODIFICATION OF THE NODE NUMBER
45    READ(MR,1040)IHOW
1040  FORMAT(I1)    
      IF (IHOW.EQ.0) GOTO 20

C---------READ NODE CARDS
50    READ(MR,1020)(KV(IN),IN=1,L16)
1020  FORMAT(16I6)
      IF (M.GE.0) WRITE(MP,2040)(KV(IN),IN=1,L16)
2040  FORMAT('>>>>>>',10X,16I6)

C-----------DIFFERENT BEHAVIOUR AS FUNCTION OF IHOW
      IF (IHOW.EQ.1) THEN
C-----------FORM NEQ
      DO 100 IN=1,L16
	I2=KV(IN)

C-----------END OF GROUP OF B.C. OR END OF NODES OR ANALYSIS OF A NODE
	IF (I2) 45,45,60
60	IF (I2.GT.NNT) CALL ERREUR(32,I2,NNT,1)
	I1=KDLNC(I2)
	IDN=KDLNC(I2+1)-I1

C-----------GENERATE VDIMP, PUT IT IN KNEQ (THE PRESCIBED D.O.F. ADDRESS)
	IV=0
	DO 90 ID=1,IDN
	  I1=I1+1
	  IC=ICOD(ID)-1
	  IF (IC) 90,70,80
70	  NCLT=NCLT+1
	  VDIMP(NCLT)=ZERO
	  NCLZ=NCLZ+1
	  KNEQ(I1)=-NCLT
	  GO TO 90

80	  NCLT=NCLT+1
	  IV=IV+1
	  VDIMP(NCLT)=V(IV)
	  NCLNZ=NCLNZ+1
	  KNEQ(I1)=-NCLT
90	CONTINUE
100   CONTINUE

C------------ADDITIONAL CARD OF NODE NUMBERS
      GO TO 50

      ELSE

      IF (KV(3).EQ.0)KV(3)=1

      DO 155,I2=KV(1),KV(2),KV(3)
      IF (I2) 45,45,210
210   IF (I2.GT.NNT) CALL ERREUR(32,I2,NNT,1)

      I1=KDLNC(I2)
      IDN=KDLNC(I2+1)-I1
      IV=0
	DO 290 ID=1,IDN
	  I1=I1+1
	  IC=ICOD(ID)-1
	  IF (IC) 290,270,280
270	  NCLT=NCLT+1
	  VDIMP(NCLT)=ZERO
	  NCLZ=NCLZ+1
	  KNEQ(I1)=-NCLT
	  GO TO 290

280	  NCLT=NCLT+1
	  IV=IV+1
	  VDIMP(NCLT)=V(IV)
	  NCLNZ=NCLNZ+1
	  KNEQ(I1)=-NCLT
290	CONTINUE
155	CONTINUE

      ENDIF
C---------ADDITIONNAL INITIAL CONDITIONS
      GOTO 45

C-------------GENERATE EQUATION NUMBERS IN NEQ
110   I1=0
      DO 150 IN=1,NNT
	ID=KDLNC(IN)
120	ID=ID+1
	IF (ID.GT.KDLNC(IN+1)) GO TO 150
	IF (KNEQ(ID)) 120,130,120
130	I1=I1+1
	KNEQ(ID)=I1
	GO TO 120
150   CONTINUE
      NEQ=I1

C-------------OUTPUT
      IF (M.LT.0) GO TO 170
      WRITE(MP,2050)NNT,NDLT,NEQ,NCLNZ,NCLZ,NCLT
2050  FORMAT(//
     1	15X,'TOTAL NUMBER OF NODES                     (NNT)=',I5/
     2	15X,'TOTAL NUMBER OF D.O.F                    (NDLT)=',I5/
     3	15X,'NUMBER OF EQUATIONS TO BE SOLVED          (NEQ)=',I5/
     4	15X,'NUMBER OF PRESCRIBED NON ZERO D.O.F.    (NCLNZ)=',I5/
     5	15X,'NUMBER OF PRESCRIBED ZERO D.O.F.         (NCLZ)=',I5/
     6	15X,'TOTAL NUMBER OF PRESCRIBED D.O.F.        (NCLT)=',I5/)

      IF ((M.GE.2).AND.(NCLT.GT.0)) WRITE(MP,2060)(VDIMP(I),I=1,NCLT)
2060  FORMAT(//'PRESCRIBED VALUES  (VDIMP)'//(10X,10E12.5))

      WRITE(MP,2070)
2070  FORMAT(//' NODAL COORDINATES ARRAY'//
     1	'   NO  D.L.',5X,'X',12X,'Y',12X,'Z','EQUATION NUMBER (NEQ)'/)

      I2=0
      DO 160 IN=1,NNT
	I1=I2+1
	I2=I2+NDIM
	ID1=KDLNC(IN)+1
	ID2=KDLNC(IN+1)
	ID=ID2-ID1+1
	IF (ID2.LT.ID1) ID2=ID1
	X1=VCORG(I1)
	IF (NDIM.GE.2) X2=VCORG(I1+1)
	IF (NDIM.GE.3) X3=VCORG(I1+2)
	WRITE(MP,2080) IN,ID,X1,X2,X3,(KNEQ(I),I=ID1,ID2)
160   CONTINUE
2080  FORMAT(1X,2I5,3E12.5,10X,10I6)


170   RETURN
      END

