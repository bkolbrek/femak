      SUBROUTINE FNDELM(VPRSG,KTYPRSG,VNDSRC,KLOCE,VCORE,KNE,
     1              VPRNE,VPREE,IEL0)
C************************************************************
C     FNDELM.FOR
C     Tonni F. Johansen
C     Kompaniet Numerisk Lyd
C     27/3-1995
C
C************************************************************
C     Search for an element and write source specifications
C************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include 'femak.cmn'
c  $DECLARE 
 
      COMPLEX*16 VPRSG(1),VNDSRC(1)
      REAL*8 VCORE(1),VPRNE(1),VPREE(1)
      INTEGER KLOCE(1),KNE(1),KTYPRSG(1)
      INTEGER IEL0,NSEL

C
C...........Search for information about ELEMNR in elementfile
        
      REWIND(M1)
      IEL=0
c--------------Search for "real-elements"
      NSEL=0

      DO WHILE ((IEL0.NE.IEL).AND.(NSEL.LT.NELTR))
        CALL RDELEM(M1,KLOCE,VCORE,VPRNE,VPREE,KNE)
        NSEL=NSEL+1
      ENDDO     

C............If not among "real-elements" searc "complex-elements"
      IF (IEL.NE.IEL0) THEN
        REWIND(M2)
        IEL=0
C...........Serch "complex-elements"
        NSEL=0
        DO WHILE ((IEL0.NE.IEL).AND.(NSEL.LT.NELTC))
          CALL RDELEM(M2,KLOCE,VCORE,VPRNE,VPREE,KNE)
          NSEL=NSEL+1
        ENDDO     
          
      ENDIF          



C...........If found, write element properties to source file
      IF (IEL.EQ.IEL0) THEN
        CALL WRSRCE(M7,KLOCE,VCORE,VPRNE,VPREE,KNE,VNDSRC,KTYPRSG)
        IF (M.GT.0) THEN
          CALL PRSRCE(2,KTYPRSG(IGSRC),VPRSG,VNDSRC)
        ENDIF
      ELSE
        WRITE(*,1060)IEL0
1060    FORMAT('ERROR; IEL0 ,'I4,' DOES NOT EXIST !')
      ENDIF
C!!!NEEd appropriate error messages
      
      RETURN
      END


