      SUBROUTINE EXSRCE(VPRSG,KTYPRSG,VNDSRC,KLOCE,VCORE,KNE,
     1              VPRNE,VPREE,VPRNG)
C************************************************************C
C     EXSRCE.FOR                                             C
C     Tonni F. Johansen, NTH                                 C
C     24/8-1992                                              C
C************************************************************C
C     TO EXECUTE BLOCK 'SRCE'                                C
C     TO READ SOURCE CONDITIONS                              C
C************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)
      include 'femak.cmn'
c   $DECLARE 
 
      COMPLEX*16 VPRSG(1),VNDSRC(1)
      REAL*8 VCORE(1),VPRNE(1),VPREE(1),VPRNG(1)
      INTEGER KLOCE(1),KNE(1),KTYPRSG(1),NODGRP
      INTEGER SRCPAR,I,J,I2,J2,II,JJ,IK,IEL0
      INTEGER init
c      integer NSEL
c      LOGICAL ELFOUND
      
      REAL*8 ABSK,MAXNODGRP

      II=0
      JJ=0
      init = 0
C-----------------------------------------------------------
C------- READ SOURCE GROUP, SOURCE PARAMETER AND
C        SOURCE DESCRIPTION
      DO IK=1,NGSRC
C        READ(M1,1000)IGSRC,SRCPAR
        READ(MR,1000)IGSRC,SRCPAR
1000    FORMAT(2I5)
        IF(M.GT.0) THEN
          CALL PRSRCE(1,SRCPAR,VPRSG,VNDSRC)
        ENDIF
1020    FORMAT(6E12.5)
        KTYPRSG(IGSRC)=SRCPAR
C.........Her m† det komme inn muligheter for † se p† spesieller typer kilder
C.........etterhvert.
C..........pr.30/8-92 konsentrerer jeg meg om den mest generelle, dvs. den som
C.........angir nodehastihet p†randa i detalj. (SRCPAR=20)
C..........Den trenger ingen videre innlesning av data her.       
        IF ((SRCPAR.EQ.0).OR.(SRCPAR.EQ.3)) THEN
          IF (MOD(NFREQ,4).EQ.0) THEN
            I2=NFREQ/4
          ELSE
            I2=NFREQ/4+1
          ENDIF
          DO 10 I=1,I2
            J2=MIN(4,NFREQ-(I-1)*4)
C            READ(M1,1010)(VPRSG((IGSRC-1)*NFREQ+(I-1)*4+J),J=1,J2)
            READ(MR,1010)(VPRSG((IGSRC-1)*NFREQ+(I-1)*4+J),J=1,J2)
10        CONTINUE
          IF(M.GT.0) THEN
            CALL PRSRCE(1,SRCPAR,VPRSG,VNDSRC)
          ENDIF
1010      FORMAT(8E12.5)
CMARC
C========== type source 1
        ELSEIF (SRCPAR.EQ.1) THEN
          READ(MR,'(4F8.4)') vit_son1,rhof1,vit_infi,del_cl
          READ(MR,'(2I5)') dir_xy,excit
          READ(MR,'(5I5)') idsp1,idsp2,idsp3,idsp4,idsp5        
        ENDIF
C========== END type source 1
CMARC        


      ENDDO
C........Write text for next input
      IF(M.GT.0) THEN
        CALL PRSRCE(0,0,VPRSG,VNDSRC)
      ENDIF

      DO IK=1,NSRC
        READ(MR,1012)IEL0,IFL,IGSRC
1012    FORMAT(3I5)
        IF ((KTYPRSG(IGSRC).EQ.0).OR.(KTYPRSG(IGSRC).EQ.3)) THEN

C.......write source specifications to file 
C.......(dummy element for source type 3)
          CALL FNDELM(VPRSG,KTYPRSG,VNDSRC,KLOCE,VCORE,KNE,
     1              VPRNE,VPREE,IEL0)

          IF (M.GT.0) THEN
            CALL PRSRCE(2,KTYPRSG(IGSRC),VPRSG,VNDSRC)
          ENDIF

        ELSEIF (KTYPRSG(IGSRC).EQ.1) THEN
C...........SRCPAR=1        
          IF (init.EQ.0) idlinf = ifl
          init   = 1
          idlsup = ifl
          
          CALL FNDELM(VPRSG,KTYPRSG,VNDSRC,KLOCE,VCORE,KNE,
     1              VPRNE,VPREE,IEL0)

          IF (M.GT.0) THEN
            CALL PRSRCE(2,KTYPRSG(IGSRC),VPRSG,VNDSRC)
          ENDIF        
        
        ELSEIF (KTYPRSG(IGSRC).EQ.20) THEN
C...........SRCPAR=20
          IF (NDIM.EQ.2) THEN 
            DO I=1,NFREQ
              READ(MR,1050)(VNDSRC((I-1)*4+J),J=1,3)
            ENDDO
          ELSEIF (NDIM.EQ.3) THEN 
            DO I=1,NFREQ
              READ(MR,1055)(VNDSRC((I-1)*10+J),J=1,8)
            ENDDO
          ENDIF
1050      FORMAT(6E12.5)
1055      FORMAT(16E12.5)

C.......write source specifications to file
          CALL FNDELM(VPRSG,KTYPRSG,VNDSRC,KLOCE,VCORE,KNE,
     1              VPRNE,VPREE,IEL0)

        ELSEIF (KTYPRSG(IGSRC).EQ.21) THEN
C...........SRCPAR=21
          IF (NDIM.EQ.2) THEN 

C..........Read input lines
            READ(MR,1050)(VNDSRC(J),J=1,3)
            
C...........Duplicate for equivalence with SRCPAR=20
            DO I=2,NFREQ
              DO J=1,3
                VNDSRC((I-1)*4+J)=VNDSRC(J)
                WRITE(MP,2099)VNDSRC((I-1)*4+J)
2099  format(2e12.5)                                            
              ENDDO
            ENDDO
          ELSEIF (NDIM.EQ.3) THEN 
            READ(MR,1055)(VNDSRC(J),J=1,8)
            DO I=2,NFREQ
              DO J=1,8
                VNDSRC((I-1)*10+J)=VNDSRC(J)
              ENDDO
            ENDDO
          ENDIF
C.......write source specifications to file
          CALL FNDELM(VPRSG,KTYPRSG,VNDSRC,KLOCE,VCORE,KNE,
     1              VPRNE,VPREE,IEL0)


        ELSEIF (KTYPRSG(IGSRC).EQ.25) THEN
C...........SRCPAR=25 
          IF (NDIM.EQ.2) THEN 
            DO I=1,NFREQ
              READ(MR,1051)VNDSRC((I-1)*4+4)
C......Check for directional cosines .ne. 1
              IF (ABS(VNDSRC((I-1)*4+4)).NE.1.0) THEN
C.........Scaling if not equal to 1.0
                VNDSRC((I-1)*4+4)=VNDSRC((I-1)*4+4)/
     1                        ABS(VNDSRC((I-1)*4+4))
              ENDIF
              READ(MR,1050)(VNDSRC((I-1)*4+J),J=1,3)
            ENDDO
          ELSEIF (NDIM.EQ.3) THEN 
            DO I=1,NFREQ
              READ(MR,1056)(VNDSRC((I-1)*10+J),J=9,10)
              READ(MR,1055)(VNDSRC((I-1)*10+J),J=1,8)
            ENDDO
          ENDIF
1051      FORMAT(2E12.5)
1056      FORMAT(4E12.5)
C.......write source specifications to file
          CALL FNDELM(VPRSG,KTYPRSG,VNDSRC,KLOCE,VCORE,KNE,
     1              VPRNE,VPREE,IEL0)

        ELSEIF (KTYPRSG(IGSRC).EQ.26) THEN
C...........SRCPAR=26
          IF (NDIM.EQ.2) THEN 

C..........Read input lines
            READ(MR,1051)VNDSRC(4)
C......Check for directional cosines .ne. 1
            IF (ABS(VNDSRC(4)).NE.1.0) THEN
C.........Scaling if not equal to 1.0
              VNDSRC(4)=VNDSRC(4)/ABS(VNDSRC(4))
            ENDIF
            READ(MR,1050)(VNDSRC(J),J=1,3)
            
C...........Duplicate for equivalence with SRCPAR=25
            DO I=2,NFREQ
              DO J=1,4
                VNDSRC((I-1)*4+J)=VNDSRC(J)
              ENDDO
            ENDDO
          ELSEIF (NDIM.EQ.3) THEN 
            READ(MR,1056)(VNDSRC(J),J=9,10)
            READ(MR,1055)(VNDSRC(J),J=1,8)
            DO I=2,NFREQ
              DO J=1,10
                VNDSRC((I-1)*10+J)=VNDSRC(J)
              ENDDO
            ENDDO
          ENDIF
C.......write source specifications to file
          CALL FNDELM(VPRSG,KTYPRSG,VNDSRC,KLOCE,VCORE,KNE,
     1              VPRNE,VPREE,IEL0)

        ELSEIF (KTYPRSG(IGSRC).EQ.50) THEN
C...........SRCPAR=50, incoming plane wave

C...........Set Incoming wave flag
          INWVE=.TRUE.
          
          IF (NDIM.EQ.2) THEN 

C..........Read input lines
            
            READ(MR,1052)NODGRP
            VPRSG((IGSRC-1)*NFREQ+1)=FLOAT(NODGRP)
            READ(MR,1050)(VNDSRC(J),J=1,3)
C......Check for directional cosines .ne. 1
            ABSK=ABS(SQRT(VNDSRC(2)*VNDSRC(2)+VNDSRC(3)*VNDSRC(3)))
            IF (ABSK.NE.1.0) THEN
C.........Scaling if not equal to 1.0
              VNDSRC(2)=VNDSRC(2)/ABSK
              VNDSRC(3)=VNDSRC(3)/ABSK
            ENDIF
            
C...........Duplicate for equivalence with SRCPAR=20
            DO I=2,NFREQ
              VPRSG((IGSRC-1)*NFREQ+I)=VPRSG((IGSRC-1)*NFREQ+1)
              DO J=1,3
                VNDSRC((I-1)*4+J)=VNDSRC(J)
              ENDDO
            ENDDO
          ENDIF
1052  FORMAT(I5)          
C.......write source specifications to file
          CALL FNDELM(VPRSG,KTYPRSG,VNDSRC,KLOCE,VCORE,KNE,
     1              VPRNE,VPREE,IEL0)


        ELSEIF (KTYPRSG(IGSRC).EQ.55) THEN
C...........SRCPAR=55, incoming plane wave

C...........Set Incoming wave flag
          INWVE=.TRUE.
          
          IF (NDIM.EQ.2) THEN 

C..........Read input lines
            
            READ(MR,1052)NODGRP
            VPRSG((IGSRC-1)*NFREQ+1)=FLOAT(NODGRP)
            READ(MR,1050)(VNDSRC(J),J=1,3)
C......Check for directional cosines .ne. 1
            ABSK=ABS(SQRT(VNDSRC(2)*VNDSRC(2)+VNDSRC(3)*VNDSRC(3)))
            IF (ABSK.NE.1.0) THEN
C.........Scaling if not equal to 1.0
              VNDSRC(2)=VNDSRC(2)/ABSK
              VNDSRC(3)=VNDSRC(3)/ABSK
            ENDIF
            
C...........Duplicate for equivalence with SRCPAR=20
            DO I=2,NFREQ
              VPRSG((IGSRC-1)*NFREQ+I)=VPRSG((IGSRC-1)*NFREQ+1)
              DO J=1,3
                VNDSRC((I-1)*4+J)=VNDSRC(J)
              ENDDO
            ENDDO
          ENDIF

C.........Find nodes in infinite elements and assign them to a new nodegroup
C.........First find max.group number
          MAXNODGRP=0.0
          DO I=1,NNT
            MAXNODGRP=MAX(MAXNODGRP,VPRNG((I-1)*NPRN+1))
          ENDDO
          VPRSG((IGSRC-1)*NFREQ+1)=
     1             VPRSG((IGSRC-1)*NFREQ+1)+CMPLX(0.0,MAXNODGRP)

C.........Find element of type 11 (always complex elements)          
          REWIND(M2)
          DO I=1,NELTC
            CALL RDELEM(M2,KLOCE,VCORE,VPRNE,VPREE,KNE)
            IF (ITPE.EQ.11) THEN
              IF (IDINT(REAL(VPRSG((IGSRC-1)*NFREQ+1))).EQ.
     1              (IDINT(VPRNG(KNE(1)-1)*NPRN+1))) THEN         
C.................Checking if element in contact with incoming wave
C.................if yes all outer nodes will be assigned to the new
C.................nodegroup (MAXNODGRP)
                IF (NDIM.EQ.2) THEN
                  DO J=2,6
                    VPRNG(KNE(J))=IDINT(MAXNODGRP)
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          
C.......write source specifications to file
          CALL FNDELM(VPRSG,KTYPRSG,VNDSRC,KLOCE,VCORE,KNE,
     1              VPRNE,VPREE,IEL0)

          
C...........2. change boundary condition arrays and variables.
C...........3. more to be changed ????!!!!
C

        ENDIF
        
        
C....Next section is replaced by the subroutine FNDELM
C...........Search for information about ELEMNR in elementfile
        
c          REWIND(M1)
c          IEL=0
c--------------Search for "real-elements"
c          NSEL=0

c          DO WHILE ((IEL0.NE.IEL).AND.(NSEL.LT.NELTR))
c            CALL RDELEM(M1,KLOCE,VCORE,VPRNE,VPREE,KNE)
c            NSEL=NSEL+1
c          ENDDO     

C............If not among "real-elements" searc "complex-elements"
c          IF (IEL.NE.IEL0) THEN
c            REWIND(M2)
c            IEL=0
C...........Serch "complex-elements"
c            NSEL=0
c            DO WHILE ((IEL0.NE.IEL).AND.(NSEL.LT.NELTC))
c              CALL RDELEM(M2,KLOCE,VCORE,VPRNE,VPREE,KNE)
c              NSEL=NSEL+1
c            ENDDO     
          
c          ENDIF          



C...........If found, write element properties to source file
c          IF (IEL.EQ.IEL0) THEN
c            CALL WRSRCE(M7,KLOCE,VCORE,VPRNE,VPREE,KNE,VNDSRC,KTYPRSG)
c            IF (M.GT.0) THEN
c              CALL PRSRCE(2,KTYPRSG(IGSRC),VPRSG,VNDSRC)
c            ENDIF
c          ELSE
c            WRITE(*,1060)IEL0
c1060        FORMAT('ERROR; IEL0 ,'I4,' DOES NOT EXIST !')
c          ENDIF
C!!!NEEd appropriate error messages
      ENDDO
      
      RETURN
      END


