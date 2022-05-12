C############################################################################
C                           SOR method for equations solving                #
c                            VERSION 1.0 (24/05/2012)                       #
C                            AUTHORIZED BY ZHANG JINGXIN                    #
C                            SHANGHAI JIAO TONG UNIVERSITY                  #
C                                 SHANGHAI, CHINA                           #
c                                                                           #
c############################################################################
	Subroutine SOLVEDYN
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-12,EPSI=1.E-6,OMIGA=1.25,INUM=100)
	Common/DYNBLK/AS(IJM,KB,IPOLYGEN),AB(IJM,KB),AT(IJM,KB),
     &	          AP(IJM,KB),BB(IJM,KB),X(IJM,KB)
	Common/DYNPRE/ASL(IJM,KB,IPOLYGEN),ABL(IJM,KB),ATL(IJM,KB),
     &	          APL(IJM,KB),ASR(IJM,KB,IPOLYGEN),ABR(IJM,KB),
     &              ATR(IJM,KB),APR(IJM,KB),BBP(IJM,KB),XP(IJM,KB)
	Dimension R1(IJM,0:KB),U1(IJM,0:KB),V1(IJM,KB),P1(IJM,0:KB)
	Dimension TRACE(IJM,KB)
C===========================================================================C
C                optimazing the Matrix                                      c
c===========================================================================c
	Call PREMATRIX                                 ! Preconditioning Matrix
c-----  INITIAL VALUES	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
CCC	         R1(I,K) = BB(I,K)
	         R1(I,K) = PN(I,K)
	         X(I,K) = PN(I,K)
CC	         U1(I,K) = 0.0
	         V1(I,K) = 0.0
	         P1(I,K) = 0.0
	         U1(I,K) = OMIGA * BB(I,K)
	      Enddo
	   Endif
	Enddo
c	Call SOLVETRI(OMIGA)
c	Do I = 1, IJM
c	   If(CCM(I) .EQ. 1.0) Then
c	      Do K = 1, KBM
c	         U1(I,K) = XP(I,K)
c	      Enddo
c	   Endif
c	Enddo
CC	Do I = 1, IJM
CC	   If(CCM(I) .EQ. 1.0) Then
CC	      Do K = 1, KBM
CC	         U1(I,K) = XP(I,K)
CC	         P1(I,K) = U1(I,K)
ccc	         ROU = ROU + R1(I,K) * U1(I,K)
CC	      Enddo
CC	   Endif
CC	Enddo
c	ROUINI = ROU
cc	CALL CPU_TIME(END1)
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
	     Do K = 1, KBM
	        BBP(I,K) = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	              BBP(I,K) = BBP(I,K) + OMIGA * 
     &				  ASR(I,K,J) * R1(CELL_SIDE(I,J,2),K)
	           Endif
	        Enddo
	        BBP(I,K) = BBP(I,K) + OMIGA * 
     & 			(ATR(I,K) * R1(I,K-1) + ABR(I,K) * R1(I,K+1))
	        BBP(I,K) = BBP(I,K) + (1. - OMIGA) * AP(I,K) * R1(I,K)
	     Enddo
	   Endif
	Enddo
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
               BBP(I,K) = BBP(I,K) + U1(I,K)
	      Enddo
	   Endif
	Enddo 
	Call SOLVETRI(OMIGA)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         V1(I,K) = XP(I,K)
	      Enddo
	   Endif
	Enddo
C-----   CONVERGENCE
	ERROR = 0.0
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         AR = Abs(V1(I,K)-X(I,K)) / (1. + Abs(V1(I,K)))
	         If(AR .GT. ERROR) ERROR = AR	         
	      Enddo
	   Endif
	Enddo
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         X(I,K) = V1(I,K)
	         R1(I,K) = X(I,K)
	      Enddo
	   Endif
	Enddo
CCC	PRINT*, ERROR, KNUM
	If(ERROR .LT. EPSI .OR. KNUM .GE. INUM) Return
	Goto 10
	Return
	End
c===========================================================================c
c         preconditioning the matrix for SOR method                         c
c===========================================================================c
	Subroutine PREMATRIX
      Include './Include/OCERM_INF'
	Common/DYNBLK/AS(IJM,KB,IPOLYGEN),AB(IJM,KB),AT(IJM,KB),
     &	          AP(IJM,KB),BB(IJM,KB),X(IJM,KB)
	Common/DYNPRE/ASL(IJM,KB,IPOLYGEN),ABL(IJM,KB),ATL(IJM,KB),
     &	          APL(IJM,KB),ASR(IJM,KB,IPOLYGEN),ABR(IJM,KB),
     &              ATR(IJM,KB),APR(IJM,KB),BBP(IJM,KB),XP(IJM,KB)
	Do K = 1, KB
	   Do I = 1, IJM
	      Do J = 1, CELL_POLYGEN(I)
	         ASL(I,K,J) = 0.0
	         ASR(I,K,J) = 0.0
	      Enddo
		  ABL(I,K) = 0.0
		  ATL(I,K) = 0.0
		  APL(I,K) = 0.0
		  ABR(I,K) = 0.0
		  ATR(I,K) = 0.0
		  APR(I,K) = 0.0
	   Enddo
	Enddo
C-----  MATRIX L AND U
	Do K = 1, KBM
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         Do J = 1, CELL_POLYGEN(I)
	            If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	               INDEXCELL = CELL_SIDE(I,J,2)
	               If(INDEXCELL .LT. I) Then
		              ASL(I,K,J) = AS(I,K,J)
		           Endif
	               If(INDEXCELL .GT. I) Then
	                  ASR(I,K,J) = AS(I,K,J)
	               Endif
	            Endif
	         Enddo
	         ATL(I,K) = AT(I,K)
	         ABR(I,K) = AB(I,K) 
	      Endif
	   Enddo
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
	Subroutine SOLVETRI(OMIGA)
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-12,EPSI=1.E-4,INUM=200)
	Common/DYNBLK/AS(IJM,KB,IPOLYGEN),AB(IJM,KB),AT(IJM,KB),
     &	          AP(IJM,KB),BB(IJM,KB),X(IJM,KB)
	Common/DYNPRE/ASL(IJM,KB,IPOLYGEN),ABL(IJM,KB),ATL(IJM,KB),
     &	          APL(IJM,KB),ASR(IJM,KB,IPOLYGEN),ABR(IJM,KB),
     &              ATR(IJM,KB),APR(IJM,KB),BBP(IJM,KB),XP(IJM,KB)
	Dimension R1(IJM,0:KB),V1(IJM,KB),P1(IJM,0:KB)
C-----   INITIAL VALUES
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      P1(I,0) = 0.0
	      R1(I,0) = 0.0
	      Do K = 1, KB
	         R1(I,K) = 0.0
	         V1(I,K) = 0.0
	         P1(I,K) = 0.0
	      Enddo
	   Endif
	Enddo
C-----   FORWARD
	Do K = 1, KBM
         Do I = 1, IJM
		  If(CCM(I) .EQ. 1.0) Then
	        V1(I,K) = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
                   V1(I,K) = V1(I,K) + OMIGA * 
     &				 ASL(I,K,J) * P1(CELL_SIDE(I,J,2),K)
	           Endif
	        Enddo
	        V1(I,K)= V1(I,K) + OMIGA * 
     & 			(ATL(I,K) * P1(I,K-1) + ABL(I,K) * P1(I,K+1))
	        V1(I,K)= BBP(I,K) + V1(I,K)
	        P1(I,K) = V1(I,K) / AP(I,K)
	     Endif
	   Enddo
	Enddo
	Do K = 1, KBM
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         XP(I,K) = P1(I,K)
	      Endif
	   Enddo
	Enddo
	Return
	End