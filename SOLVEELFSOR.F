C############################################################################
C                           SOR method for equations solving                #
c                            VERSION 1.0 (24/05/2012)                       #
C                            AUTHORIZED BY ZHANG JINGXIN                    #
C                            SHANGHAI JIAO TONG UNIVERSITY                  #
C                                 SHANGHAI, CHINA                           #
c                                                                           #
c############################################################################
	Subroutine SOLVEELFSOR
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-12,EPSI=1.E-6,OMIGA=1.2,INUM=500)
	Common/ELFBLK/CS(IJM,IPOLYGEN),CB(IJM),CP(IJM),X(IJM)

	Common/ELFPRE/ASL(IJM,IPOLYGEN),APL(IJM),
     &	          ASR(IJM,IPOLYGEN),APR(IJM),BBP(IJM),XP(IJM)

	Dimension R1(IJM),U1(IJM),V1(IJM),P1(IJM)
	Dimension TRACE(IJM)
C===========================================================================C
C                optimazing the Matrix                                      c
c===========================================================================c
	Call PREMATRIXSORELF                           ! Preconditioning Matrix
c-----  INITIAL VALUES	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      R1(I) = 0.0
	      X(I) = 0.0
CC	         U1(I,K) = 0.0
	      V1(I) = 0.0
	      P1(I) = 0.0
	      U1(I) = OMIGA * CB(I)
	   Endif
	Enddo
C---------------------------------------------------------------------------C
      KNUM = 0
CC	ALPHA = 1.0
CC	ROU = 1.0
10    Continue
      KNUM = KNUM + 1
C	ROU = 0.0
CC	ROUT = 0.0
C-----   PART1
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      BBP(I) = 0.0
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	            BBP(I) = BBP(I) + OMIGA * 
     &				  ASR(I,J) * R1(CELL_SIDE(I,J,2))
	         Endif
	      Enddo
	      BBP(I) = BBP(I) + (1. - OMIGA) * CP(I) * R1(I)
	   Endif
	Enddo
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
            BBP(I) = BBP(I) + U1(I)
	   Endif
	Enddo 
	Call SOLVETRISORELF(OMIGA)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      V1(I) = XP(I)
	   Endif
	Enddo
C-----   CONVERGENCE
	ERROR = 0.0
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      AR = Abs(V1(I)-X(I)) / (1. + Abs(V1(I)))
	      If(AR .GT. ERROR) ERROR = AR	         
	   Endif
	Enddo
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      X(I) = V1(I)
	      R1(I) = X(I)
	   Endif
	Enddo
	PRINT*, ERROR, KNUM
	If(ERROR .LT. EPSI .OR. KNUM .GE. INUM) Return
	Goto 10
	Return
	End
c===========================================================================c
c         preconditioning the matrix for SOR method                         c
c===========================================================================c
	Subroutine PREMATRIXSORELF
      Include './Include/OCERM_INF'
	Common/ELFBLK/CS(IJM,IPOLYGEN),CB(IJM),CP(IJM),X(IJM)

	Common/ELFPRE/ASL(IJM,IPOLYGEN),APL(IJM),
     &	          ASR(IJM,IPOLYGEN),APR(IJM),BBP(IJM),XP(IJM)

	Do I = 1, IJM
	   Do J = 1, CELL_POLYGEN(I)
	      ASL(I,J) = 0.0
	      ASR(I,J) = 0.0
	   Enddo
	   APL(I) = 0.0
	   APR(I) = 0.0
	Enddo
C-----  MATRIX L AND U
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	            INDEXCELL = CELL_SIDE(I,J,2)
	            If(INDEXCELL .LT. I) Then
		           ASL(I,J) = CS(I,J)
		        Endif
	            If(INDEXCELL .GT. I) Then
	               ASR(I,J) = CS(I,J)
	            Endif
	         Endif
	      Enddo
	   Endif
	Enddo
	Return
	End

C############################################################################
C                     method for tri-dignal equations solving               #
c                            VERSION 1.0 (21/05/2012)                       #
C                            AUTHORIZED BY ZHANG JINGXIN                    #
C                            SHANGHAI JIAO TONG UNIVERSITY                  #
C                                 SHANGHAI, CHINA                           #
c                                                                           #
c############################################################################
	Subroutine SOLVETRISORELF(OMIGA)
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-12,EPSI=1.E-5,INUM=200)
	Common/ELFBLK/CS(IJM,IPOLYGEN),CB(IJM),CP(IJM),X(IJM)

	Common/ELFPRE/ASL(IJM,IPOLYGEN),APL(IJM),
     &	          ASR(IJM,IPOLYGEN),APR(IJM),BBP(IJM),XP(IJM)
	Dimension R1(IJM),V1(IJM),P1(IJM)
C-----   INITIAL VALUES
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      R1(I) = 0.0
	      V1(I) = 0.0
	      P1(I) = 0.0
	   Endif
	Enddo
C-----   FORWARD
      Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      V1(I) = 0.0
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
                  V1(I) = V1(I) + OMIGA * 
     &				 ASL(I,J) * P1(CELL_SIDE(I,J,2))
	         Endif
	      Enddo
	      V1(I)= BBP(I) + V1(I)
	      P1(I) = V1(I) / CP(I)
	   Endif
	Enddo
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      XP(I) = P1(I)
	   Endif
	Enddo
	Return
	End