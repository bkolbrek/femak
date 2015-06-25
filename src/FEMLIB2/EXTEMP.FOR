      SUBROUTINE EXTEMP(VCORG,KDLNC,VDIMP,KNEQ,KLD,KLOCE,VCORE,VPRNE,
     1  VPREE,KNE,VKE,VME,VFE,VDLE,VKGS,VKGD,VKGI,VFG,VRES,VDLG,
     2  VDLE0,VDLG0,VFG0)
C*********************************************************************C
C  EXTEMP.FOR                                                         C
C  02/12-1994                                                         C
C  Bruno Brouard NTH, LAUM                                            C
C*********************************************************************C
C  To execute  Block 'TEMP'                                           C
C  Programmed after Dhatt & Touzot, p.431-432                         C
C  Modified to read initial conditions in FEMIN.DAT                   C
C*********************************************************************C
C     TO SOLVE AN UNSTEADY PROBLEM (LINEAR OR NOT)
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      include 'femak.cmn'

      DIMENSION VCORG(1),KDLNC(1),VDIMP(1),KNEQ(1),KLD(1),KLOCE(1),
     1  VCORE(1),VPRNE(1),VPREE(1),KNE(1),VKE(1),VME(1),VFE(1),VDLE(1),
     2  VKGS(1),VKGD(1),VKGI(1),VFG(1),VRES(1),VDLG(1),VDLE0(1),
     3  VDLG0(1),VFG0(1)

      COMMON/TRVL/KV(16),V(10),ICOD(10)

      DATA L7/7/,L8/8/,L16/16/,ZERO/0.D0/,UN/1.D0/
C----------------------------------------------------------------------- 
      DPAS0=ZERO
      XPAS=ZERO
      IPAS=0
C-------  READ INITIAL D.O.F. ON FILE M3
c      IF(M3.EQ.0) GO TO 5
c      REWIND M3
c      READ(M3) (VDLG(I),I=1,NEQ)
C------- DEFAULT INITIAL CONDITION FOR ALL D.O.F.
      DO 5 I=1,NEQ
        VDLG(I)=ZERO
5     CONTINUE

      IF (M.GE.0) WRITE(MP,1000)
1000  FORMAT(// 'INITIAL CONDITION CARDS '/)

C----------READ AN INITIAL CONDITION GROUP CARD : 10 CODES + PRESCRIBED VALUE
20    READ(MR,1010)ICOD,(V(I),I=1,L7)
1010  FORMAT(10I1,7F10.0)
      IF (M.GE.0) WRITE(MP,2000)ICOD,(V(I),I=1,L7)
2000  FORMAT('##### ',10I1,7E12.5)

C----------CHECK FOR A BLANK CARD
      J=0
      DO 30 I=1,10
	J=J+ICOD(I)
30    CONTINUE
      IF (J.EQ.0) GO TO 160

C---------READ ADDITIONAL CARD IF REQUIRED
      I2=0
      DO 40 ID=1,NDLN
	IF (ICOD(ID).LT.1) GO TO 40
	I2=I2+1
	IF (I2.NE.L8) GO TO 40
	READ(MR,1030)(V(I),I=L8,NDLN)
1030	FORMAT(10X,7F10.0)
	IF (M.GE.0) WRITE(MP,2010)(V(I),I=L8,NDLN)
2010	FORMAT('#####',10X,7E12.5)
40    CONTINUE

C---------READ CODIFICATION OF THE NODE NUMBER
50    READ(MR,1040)IHOW
1040  FORMAT(I1)    
      IF (IHOW.EQ.0) GOTO 20

C---------READ NODE CARDS
60    READ(MR,1050)(KV(IN),IN=1,L16)
1050  FORMAT(16I6)
      IF (M.GE.0) WRITE(MP,2020)(KV(IN),IN=1,L16)
2020  FORMAT('>>>>>>',10X,16I6)

C-----------DIFFERENT BEHAVIOUR AS FUNCTION OF IHOW
      IF (IHOW.EQ.1) THEN
      DO 105 IN=1,L16
      I2=KV(IN)
C-----------END OF NODES OR ANALYSIS OF A NODE
      IF (I2) 50,50,70
70    IF (I2.GT.NNT) CALL ERREUR(32,I2,NNT,1)

      I1=KDLNC(I2)
      IDN=KDLNC(I2+1)-I1
      IV=0
      DO 100 ID=1,IDN
        I1=I1+1
        IC=ICOD(ID)
        IF (IC.EQ.1) THEN
        IF(KNEQ(I1).LE.0) THEN
          CALL ERREUR(61,ID,I2,1)
        ELSE
          IV=IV+1
          VDLG(KNEQ(I1))=V(IV)
        ENDIF
        ENDIF  
100   CONTINUE
105   CONTINUE
C------------ADDITIONAL CARD OF NODE NUMBERS
      GOTO 60

      ELSE
      IF (KV(3).EQ.0)KV(3)=1

      DO 150,I2=KV(1),KV(2),KV(3)
      IF (I2) 50,50,110
110   IF (I2.GT.NNT) CALL ERREUR(32,I2,NNT,1)

      I1=KDLNC(I2)
      IDN=KDLNC(I2+1)-I1
      IV=0
      DO 140 ID=1,IDN
        I1=I1+1
        IC=ICOD(ID)
        IF (IC.EQ.1) THEN
        IF(KNEQ(I1).LE.0) THEN
          CALL ERREUR(61,ID,I2,1)
        ELSE
          IV=IV+1
          VDLG(KNEQ(I1))=V(IV)
        ENDIF
        ENDIF  
140   CONTINUE
150   CONTINUE

      ENDIF
C---------ADDITIONNAL INITIAL CONDITIONS
      GOTO 50

C-------  SAVE THE REFERENCE INITIAL CONDITIONS
 160  CALL MAJ(UN,ZERO,NEQ,VDLG,VDLG0)
C-------  SAVE THE REFERENCE LOAD CONDITIONS
      CALL MAJ(UN,ZERO,NEQ,VFG,VFG0)
C-------   READ A CARD DEFINING A SET OF IDENTICAL STEPS
170   READ(MR,1060) DPAS,I1,I2,I3,X1,X2 
1060  FORMAT(F10.0,3I5,2F10.0)
      IF(DPAS.EQ.ZERO) GO TO 230
      IF(I1.GT.0) NPAS=I1
      IF(I2.GT.0) NITER=I2
      IF(I3.GT.0) IMETH=I3
      IF(X1.GT.ZERO) EPSDL=X1
      IF(X2.NE.ZERO) OMEGA=X2
C
C------- LOOP OVER THE STEPS
C
      DO 220 IP=1,NPAS
      CALL INIT(ZERO,NEQ,VFG)
      IPAS=IPAS+1
      XPAS=XPAS+DPAS
      WRITE(MP,2030) IPAS,DPAS,XPAS,NITER,IMETH,EPSDL,OMEGA
 2030 FORMAT(/1X,13('-'),'STEP NUMBER (IPAS):',I5//
     1           14X,'INCREMENT                   (DPAS)=',E12.5/
     2           14X,'TOTAL LEVEL                 (XPAS)=',E12.5/
     3           14X,'NUMBER OF ITERATIONS       (NITER)=',I12/
     4           14X,'METHOD NUMBER              (IMETH)=',I12/
     5           14X,'TOLERANCE                  (EPSDL)=',E12.5/
     6           14X,'COEFFICIENT ALPHA          (OMEGA)=',E12.5/)
C
C-------   LOOP OVER EQUILIBRIUM ITERATIONS
C
      DO 200 ITER=1,NITER
C-------   CHOOSE THE METHOD
      IF(IMETH.GT.3) GO TO 180
C-------   EULER TYPE METHODS
      CALL EULER(VCORG,KDLNC,VDIMP,KNEQ,KLD,KLOCE,VCORE,VPRNE,VPREE,
     1  KNE,VKE,VME,VFE,VDLE,VKGS,VKGD,VKGI,VFG,VRES,VDLG,
     2  VDLE0,VDLG0,VFG0)
      GO TO 190
C-------   OTHER METHODS .......
 180   CONTINUE
      WRITE(MP,2040) IMETH
 2040 FORMAT(' ** ERROR, METHOD:',I3,' UNKNOWN')
      STOP
C-------    COMPUTE THE NORM 
 190  CALL NORME(NEQ,VRES,VDLG,XNORM)      
      IF(M.GT.0) WRITE(MP,2050) ITER,XNORM
 2050 FORMAT(5X,'ITERATION (ITER):',I3,' NORM (XNORM)=',E12.5)      
      IF(M.GE.2) CALL PRSOL(KDLNC,VCORG,VDIMP,KNEQ,VDLG)
      IF(XNORM.LE.EPSDL) GO TO 210
 200  CONTINUE
C-------     END OF STEP
 210  DPAS0=DPAS
      CALL MAJ(UN,ZERO,NEQ,VDLG,VDLG0)
      CALL PRSOL(KDLNC,VCORG,VDIMP,KNEQ,VDLG)
 220  CONTINUE
      GO TO 170
C-------  SAVE THE SOLUTION ON FILE M4
 230  IF(M4.NE.0) WRITE(M4) (VDLG(I),I=1,NEQ)
      RETURN
      END
                  