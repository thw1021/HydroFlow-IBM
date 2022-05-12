C############################################################################
C                           SOR method for equations solving                #
c                            VERSION 1.0 (24/05/2012)                       #
C                            AUTHORIZED BY ZHANG JINGXIN                    #
C                            SHANGHAI JIAO TONG UNIVERSITY                  #
C                                 SHANGHAI, CHINA                           #
c                                                                           #
c############################################################################
	Subroutine SOLVEDYNSOR
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-12,EPSI=1.E-6,OMIGA=1.25,INUM=100)
	Common/DYNBLK/AS(IJM,KB,IPOLYGEN),AB(IJM,KB),AT(IJM,KB),
     &	          AP(IJM,KB),BB(IJM,KB),X(IJM,KB),XINI(IJM,KB)
	Common/DYNPRE/ASL(IJM,KB,IPOLYGEN),ABL(IJM,KB),ATL(IJM,KB),
     &	          APL(IJM,KB),ASR(IJM,KB,IPOLYGEN),ABR(IJM,KB),
     &              ATR(IJM,KB),APR(IJM,KB),BBP(IJM,KB),XP(IJM,KB)
	Dimension R1(IJM,0:KB),U1(IJM,0:KB),V1(IJM,KB),P1(IJM,0:KB)
	Dimension TRACE(IJM,KB)

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
	PRINT*, ERROR, KNUM
	If(ERROR .LT. EPSI .OR. KNUM .GE. INUM) Return
	Goto 10
	Return
	End
c===========================================================================c
c         preconditioning the matrix for SOR method                         c
c===========================================================================c
	Subroutine PREMATRIXSOR
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
	Subroutine SOLVETRISOR
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-12,EPSI=1.E-5,INUM=200)
	Common/DYNPRE/ASL(IJM,KB,IPOLYGEN),ABL(IJM,KB),ATL(IJM,KB),
     &	          APL(IJM,KB),ASR(IJM,KB,IPOLYGEN),ABR(IJM,KB),
     &              ATR(IJM,KB),APR(IJM,KB),BBP(IJM,KB),XP(IJM,KB)
	Dimension R1(IJM,0:KB),V1(IJM,KB),P1(IJM,0:KB)
C-----   INITIAL VALUES
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K)
!$OMP DO
!$OMP END DO
!$OMP BARRIER
C-----   FORWARD

!$OMP BARRIER
C-----   BACKWARD
C	R1(IJM,KBM) = P1(IJM,KBM) / P1(IJM,KBM)
	Do K = KBM, 1, -1
!$OMP DO	
	   Do I = IJM, 1, -1
	      If(CCM(I) .EQ. 1.0) Then
	        V1(I,K) = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	             V1(I,K) = V1(I,K) + ASR(I,K,J)*R1(CELL_SIDE(I,J,2),K)
	           Endif
	        Enddo
	        V1(I,K)= V1(I,K) + ATR(I,K)*R1(I,K-1) + ABR(I,K)*R1(I,K+1)
	        V1(I,K)= P1(I,K) + V1(I,K)
C	        R1(I,K) = V1(I,K) / APR(I,K)
	      Endif
	   Enddo
!$OMP END DO
	Enddo
	Do K = KBM, 1, -1
!$OMP DO	
	   Do I = IJM, 1, -1
	      If(CCM(I) .EQ. 1.0) Then
C	        V1(I,K)= V1(I,K) + ATR(I,K)*R1(I,K-1) + ABR(I,K)*R1(I,K+1)
C	        V1(I,K)= P1(I,K) + V1(I,K)
	        R1(I,K) = V1(I,K) / APR(I,K)
	      Endif
	   Enddo
!$OMP END DO
!$OMP BARRIER
	Enddo
	Do K = 1, KBM
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         XP(I,K) = R1(I,K)
	      Endif
	   Enddo
!$OMP END DO
	Enddo
!$OMP END PARALLEL
	Return
	End