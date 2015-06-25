C   FEMAK.FOR                                                          C
C   Tonni F. Johansen. 9/12-1991                                       C
C   Kompaniet Numerisk lyd/NTH                                         C
C**********************************************************************C
C   This is the main program for the finite element program: FEMAK     C
C   It is based on the program F.E.M. in Dhatt & Touzot; The Finite    C
C   Element Method Displayed.                                          C
C   Pr.date it is planned to install elements for plates, and acoustic C
C   problems. The last will be developed by us.                        C
C                                                                      C
C   The basis of the code can be found at p.378-381 in the above       C
C   mentioned book.                                                    C
C**********************************************************************C


      IMPLICIT REAL*8(A-H,O-Z)

      INTEGER ARRAYMAX
      PARAMETER (ARRAYMAX=1000000)

      REAL*4 BLOC,BLOCS
      COMMON /ALLOC/NVA,IVA,IVAMAX,NREEL,NTBL
      COMMON /ES/M,MR,MP,MLUN(10)
      COMMON VA(ARRAYMAX)
c      INTEGER UNIT(10)
c      CHARACTER*9 NAME(10)
      DIMENSION BLOCS(21)
      DATA BLOCS/4HIMAG,4HCOMT,4HCOOR,4HDLPN,4HCOND,4HPRND,4HPREL,
     &           4HELEM,4HSOLC,4HSOLR,4HLINM,4HLIND,4HNLIN,4HTEMP,
     &           4HVALP,4HCIMP,4HDYNM,4HSRCE,4H....,4H....,4HSTOP/
      DATA NB/21/

c      DATA UNIT/10*0/
c      DATA NAME/'FEM1.DAT','FEM2.DAT','FEM3.DAT','FEM4.DAT',
c     1    'FEM5.DAT','FEM6.DAT','FEM7.DAT','FEM8.DAT','FEM9.DAT',
c     2    'FEM10.DAT'/
C---------------------------------------------------------

      OPEN(MP,FILE='FEMOUT.DAT')
      OPEN(MR,FILE='FEMIN.DAT')
      OPEN(MLUN(1),FILE='FEMAK1.DAT',FORM='UNFORMATTED')
      OPEN(MLUN(2),FILE='FEMAK2.DAT',FORM='UNFORMATTED')
      OPEN(MLUN(3),FILE='FEMAK3.DAT',FORM='UNFORMATTED')
      OPEN(MLUN(4),FILE='FEMAK4.DAT',FORM='UNFORMATTED')
      OPEN(MLUN(5),FILE='FEMAK5.DAT',FORM='UNFORMATTED')
      OPEN(MLUN(6),FILE='FEMAK6.DAT',FORM='UNFORMATTED')
      OPEN(MLUN(7),FILE='FEMAK7.DAT',FORM='UNFORMATTED')
      OPEN(MLUN(8),FILE='FEMAK8.DAT',FORM='UNFORMATTED')
      OPEN(MLUN(9),FILE='FEMAK9.DAT',FORM='UNFORMATTED')
      OPEN(MLUN(10),FILE='FEMAK10.DAT',FORM='UNFORMATTED')

C.........Length of blank common in real words (table va)
      NVA=ARRAYMAX

C.........Heading
      WRITE(MP,2000)
2000  FORMAT(20X,'Acoustic Finite Element Package'/23X,
     1       'Kompaniet Numerisk Lyd'/23X' after G.Touzot, G.Dhatt'/23x,
     2       27('-')//)

C........Read BLOCK title
c10    READ(MR,1000)BLOC,M,MLUN
c1000  FORMAT(A4,I6,10I5)
10    READ(MR,1000)BLOC,M
1000  FORMAT(A4,I6)

C........Check whether MLUN includes new file to be opened.
C........PC-version only
c      DO I=1,10
c       IF(MLUN(I).NE.0) THEN
c         IF (UNIT(I).EQ.0) THEN
c           OPEN(MLUN(I),FILE=NAME(MLUN(I)),FORM='UNFORMATTED')
c           UNIT(I)=MLUN(I)
c         ENDIf
c       ENDIF
c      ENDDO

C........Search for Block to executed
      DO 20 I=1,NB
	IF (BLOC.EQ.BLOCS(I)) GO TO 30
20    CONTINUE
      WRITE(MP,2010)
2010  FORMAT('** ERROR, Missing block calling card')
      GOTO 10

30    GO TO (110,120,130,140,150,160,170,180,190,200,
     &       210,220,230,240,250,260,270,280,290,300,
     &       999),I

C..........Block to print images of data cards (I=1)           'IMAG'
110   CALL BLIMAG
      GO TO 10

C..........Block to read and print comments (I=2)              'COMT'
120   CALL BLCOMT
      GO TO 10

C..........Block to read nodal points coordinates (I=3)        'COOR'
130   CALL BLCOOR
      GO TO 10

C..........Block to read degrees of freedom per node (I=4)     'DLPN'
140   CALL BLDLPN
      GO TO 10

C..........Block to read boundary conditions (I=5)             'COND'
150   CALL BLCOND
      GO TO 10

C..........Block to read nodal properties (I=6)                'PRND'
160   CALL BLPRND
      GO TO 10

C..........Block to read element properties (I=7)              'PREL'
170   CALL BLPREL
      GO TO 10

C..........Block to read element data (I=8)                    'ELEM'
180   CALL BLELEM
      GO TO 10

C..........Block to read concentrated loads (I=9)              'SOLC'
190   CALL BLSOLC
      GO TO 10

C..........Block to read distributed loads (I=10)              'SOLR'
200   CALL BLSOLR
      GO TO 10

C..........Block for in core assembling and linear solution
C..........                                         (I=11)     'LINM'
210   continue
      CALL BLLINM
      GO TO 10

C..........Block for on disk assembling and linear solution
C..........                                         (I=12)     'LIND'
220   continue
c     CALL BLLIND
      GO TO 10

C..........Block for nonliear problem solution (I=13)          'NLIN'
230   continue
c     CALL BLNLIN
      GO TO 10

C..........Block for unsteady problem (I=14)                   'TEMP'
240   CALL BLTEMP
      GO TO 10

C..........Block to compute eigenvalues (subspace) (I=15)      'VALP'
250   CALL BLVALP
      GO TO 10

C..........Block to read impedance conditions (I=16)           'CIMP'
260   CALL BLCIMP
      GO TO 10

C..........Block to assemble and solve dynamic problems (ACOUSTIC)
C..........                               (I=17)               'DYNM'
270   CALL BLDYNM
      GO TO 10

C..........Block to read impedance conditions (I=16)           'SRCE'
280   CALL BLSRCE
      GO TO 10

C...........Undefined blocks
290   CONTINUE
300   CONTINUE
      GO TO 10

C...........End of program                         (I=21)      'STOP'
999   WRITE(MP,2020)IVAMAX,NVA
2020  FORMAT(//' End of problem,',I10,'Utilized real words over ',I10)

c      if (m.gt.2) then
c       write(mp,*)' VA HELT PAA SLUTTEN'
c       do 10000 ik=1,ivamax/10+1
c          write(mp,10010)(ik-1)*10,(va((ik-1)*10+j),j=1,10)
c10010      format(i4,10e10.3)
c10000  continue
c      endif

      STOP
      END



C**********************************************************************C
C   FEMMAIN.FOR                                                        C
C   Tonni F. Johansen. 13/12-1991                                      C
C   Kompaniet Numerisk lyd/NTH                                         C
C**********************************************************************C
C   Initialize labelled COMMON blokker.                                C
C                                                                      C
C   The code can be found at p.378-381 in the above mentioned book.    C
C**********************************************************************C

      BLOCK DATA

      IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'FEMAK.CMN'


C.......This is a test for definition of maximum size of the TRVL-block
C.......Max.size found in EXCOND (problems when using cray, Trondheim)
      INTEGER KV(16),ICOD(10)
      REAL*8 V(10),H(20)
C      COMMON/TRVL/KV(16),V(10),H(20),ICOD(10)
      COMMON/TRVL/KV,V,H,ICOD


      DIMENSION LXX(61)
      EQUIVALENCE (LXX(1),LCORG)

C.........COMMON /COOR/
      DATA NNT/20/,NDLN/2/,NDIM/2/,FAC/3*1.D0/,IAX/0/

C.........COMMON /PRND/
      DATA NPRN/0/

C.........COMMON /PREL/
      DATA NGPE/0/,NPRE/0/

C.........COMMON /ELEM/
      DATA NELT/20/,NNEL/8/,NTPE/1/,NGRE/1/,ME/1/,NIDENT/0/

C.........COMMON /ASSE/
      DATA NSYM/0/

C.........COMMON /RESO/
      DATA NRES/0/,MRES/2/

C.........COMMON /RGDT/
      DATA ITPE1/0/

C.........COMMON /LIND/
      DATA MKG1/4/,MKG2/7/

C.........COMMON /NLIN/
      DATA EPSDL/1.D-2/,OMEGA/1.D0/,DPAS/.2D0/,NPAS/1/,NITER/5/,IMETH/1/

C.........COMMON /VALP/
      DATA NITER1/10/,NMDIAG/0/,EPSLB/1.D-3/,SHIFT/0.D0/,NSS/5/,
     &     NSWM/12/,TOLJAC/1.D-12/,NVALP/3/

C.........COMMON /ES/
      DATA  MR/11/,MP/12/,M1/21/,M2/22/,M3/23/,M4/24/,M5/25/,M6/26/
      DATA M7/27/,M8/28/,M9/29/,M10/30/

C.........COMMON /ALLOC/
      DATA IVA/1/,IVAMAX/1/,NTBL/53/


C.......Define here the number of integers contained in a real     
C       for the computer employed
C        examples   : IBM SIMPLE PRECISION NREEL.EQ.1
C                     IBM DOUBLE PRECISION NREEL.EQ.2
C                     CDC                  NREEL.EQ.1
      DATA NREEL/2/
C......................

C...........COMMON /LOC/
      DATA LXX/61*1/

C...........COMMON /AKUS/
      DATA NGIMP/0/,NGSRC/0/,NFREQ/0/,NSIMP/0/,NSRC/0/
      DATA PRPAR/0/,VELPAR/0/,INTPAR/0/,NPEXT/0/
      DATA UPDTFRQ/.FALSE./,CMAT/.FALSE./,INWVE/.FALSE./

C............COMMON/IMPEL/
      DATA IFL/0/,IGIMP/0/,IFREQ/0/,IGSRC/0/,IFL1/0/
      DATA NI100NOD/0/,NI100MOD/0/,NI100FL/0/,I100RAD/0.0/
      DATA W/0.0/,W2/0.0/
      END


C**********************************************************************C
C   FEMMAIN.FOR                                                        C
C   Tonni F. Johansen. 9/12-1991                                       C
C   Kompaniet Numerisk lyd/NTH                                         C
C   SUBROUTINE ERREUR                                                  C
C**********************************************************************C
C   PRINT ERROR MESSAGES FOR BLOCKS READING DATA                       C
C                                                                      C
C   It is based on the program F.E.M. in Dhatt & Touzot; The Finite    C
C   Element Method Displayed.                                          C
C                                                                      C
C   The code can be found at p.378-381 in the above mentioned book.    C
C   Input:                                                             C
C                                                                      C
C   Output:                                                            C
C**********************************************************************C

      SUBROUTINE ERREUR(IERR,I1,I2,INIV)

      COMMON/ES/M,MR,MP,MLUN(10)
C.................................................................
C......BLOCK 'COOR'

      IF (IERR.GT.19) GO TO 200
      IE=IERR-10
      GO TO (110,120,130,140,150,160,160,180),IE

110   WRITE(MP,2110)I1,I2
2110  FORMAT('***  ERROR, FIRST NODE NUMBER(',I4,') IS GREATER THAN',
     &'NNT=',I4)
      GO TO 900

120   WRITE(MP,2120)I1,I2
2120  FORMAT('***  ERROR, SECOND NODE NUMBER(',I4,') IS GREATER THAN',
     &' NNT=',I4)
      GO TO 900

130   WRITE(MP,2130)I1,I2
2130  FORMAT('**  ERROR, NODAL NUMBER OF D.O.F.(',I4,') IS GREATER',
     &' THAN NDLN=',I4)
      GO TO 900

140   WRITE(MP,2140)
2140  FORMAT('***  ERROR, FIRST AND SECOND NODE NUMBERS ARE INCOMPATIBLE
     & WITH THE GENERATION PARAMETER')
      GO TO 900

150   WRITE(MP,2150)I1
2150  FORMAT('***  ERROR, NODE ',I4,' IS DEFINED MORE THAN ONCE')
      GO TO 900

160   WRITE(MP,2160)I1
2160  FORMAT('***  ERROR, NODE ',I4,' IS NOT DEFINED ')
      GO TO 900

180   WRITE(MP,2180)I2,I1
2180  FORMAT('**  ERROR, GENERATED NODES NUMBER(',I4,') IS LESS THAN NNT
     &=',I4)
      GO TO 900

C..........BLOCK 'DLPN'
200   IF (IERR.GT.29) GO TO 300
      IE=IERR-20
      GO TO (210,220)IE

210   WRITE(MP,2210)I1,I2
2210  FORMAT('** ERROR, NUMBER OF D.O.F (',I2,') IS GREATER THAN NDLN='
     &,I2)
      GO TO 900

220   WRITE(MP,2220)I1,I2
2220  FORMAT('** ERROR, NODE NUMBER(',I4,') IS GREATER THAN NNT='
     &,I4)
      GO TO 900


C..........BLOCK 'COND'
300   IF (IERR.GT.39) GO TO 400
      IE=IERR-30
      GO TO (900,320,900)IE
320   GO TO 220


C..........BLOCK 'PREL'
400   IF (IERR.GT.49) GO TO 500
      IE=IERR-40
      GO TO (410,900)IE

410   WRITE(MP,2410)I1,I2
2410  FORMAT('** ERROR, GROUP NUMBER (',I3,') IS GREATER THAN NGPE='
     &,I3)
      GO TO 900

C......BLOCK 'ELEM'

500   IF (IERR.GT.59) GO TO 600
      IE=IERR-50
      GO TO (510,900,530,540,550,560,570),IE

510   WRITE(MP,2510)I1,I2
2510  FORMAT('***  ERROR, NUMBER OF NODES (',I3,') IS GREATER THAN NNEL=
     &',I3)
      GO TO 900

530   WRITE(MP,2530)I1,I2
2530  FORMAT('**  ERROR, PROPERTY NUMBER (',I3,') IS GREATER THAN NGPE='
     &,I3)
      GO TO 900

540   WRITE(MP,2540)I1,I2
2540  FORMAT('***  ERROR, GROUP NUMBER (',I3,') IS GREATER THAN NGRE=',
     &I3)
      GO TO 900

550   WRITE(MP,2550)I1,I2
2550  FORMAT('***  ERROR, ELEMENT NUMBER (',I4,') IS GREATER THAN',
     &' NELT=',I4 )
      GO TO 900

560   GO TO 220

570   WRITE(MP,2570)I1,I2
2570  FORMAT('**  ERROR, NUMBER OF ELEMENTS (',I4,') IS GREATER',
     &' THAN NELT=',I4)

      GO TO 900


C..........BLOCK 'TEMP'
600   IF (IERR.GT.69) GO TO 900
      IE=IERR-60
      GO TO (610)IE

610   WRITE(MP,2610)I1,I2
2610  FORMAT('** WARNING, THE D.O.F (',I2,') OF THE NODE ',I3,'IS NOT ',
     &       'FREE ANYMORE')

C......... END
900   I1=I2
      IF (INIV.GE.2) STOP

      RETURN
      END
