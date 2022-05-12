C############################################################################
C                           SOR method for equations solving                #
c                            VERSION 1.0 (24/05/2012)                       #
C                            AUTHORIZED BY ZHANG JINGXIN                    #
C                            SHANGHAI JIAO TONG UNIVERSITY                  #
C                                 SHANGHAI, CHINA                           #
c                                                                           #
c############################################################################
	Subroutine SOLVEDYNPRESOR
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-12,EPSI=1.E-6,OMIGA=1.25,INUM=100)
	Common/DYNBLK/AS(IJM,KB,IPOLYGEN),AB(IJM,KB),AT(IJM,KB),
     &	          AP(IJM,KB),BB(IJM,KB),X(IJM,KB),XINI(IJM,KB)
	Common/DYNPRE/ASL(IJM,KB,IPOLYGEN),ABL(IJM,KB),ATL(IJM,KB),
     &	          APL(IJM,KB),ASR(IJM,KB,IPOLYGEN),ABR(IJM,KB),
     &              ATR(IJM,KB),APR(IJM,KB),BBP(IJM,KB),XP(IJM,KB)
	Dimension R1(IJM,0:KB),U1(IJM,0:KB),V1(IJM,KB),P1(IJM,0:KB)
	Dimension TRACE(IJM,KB)
C===========================================================================C
C                optimazing the Matrix                                      c
c===========================================================================c
C-----  NORMALZING THE MATRIX
C	Do I = 1, IJM
C	   If(CCM(I) .EQ. 1.0) Then
C	      Do K = 1, KBM
C		     Do J = 1, CELL_POLYGEN(I)
C	            AS(I,K,J) = AS(I,K,J) / AP(I,K)
C		     Enddo
C	         AT(I,K) = AT(I,K) / AP(I,K)
C	         AB(I,K) = AB(I,K) / AP(I,K)
C	         BB(I,K) = BB(I,K) / AP(I,K)
C	         AP(I,K) = AP(I,K) / AP(I,K)
C	      Enddo
C	   Endif
C	Enddo
	Call PREMATRIXSOR(OMIGA)                       ! Preconditioning Matrix
c-----  INITIAL VALUES	
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
CCC	         R1(I,K) = BB(I,K)
	         R1(I,K) = 0.0
	         X(I,K) = 0.0
CC	         U1(I,K) = 0.0
	         V1(I,K) = 0.0
	         P1(I,K) = 0.0
	         U1(I,K) = OMIGA * BB(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
C---------------------------------------------------------------------------C
      KNUM = 0
	AR = 1.E6
CC	ALPHA = 1.0
CC	ROU = 1.0
10    Continue
      Do While (AR .GT. EPSI .AND. KNUM .LT. INUM)
         KNUM = KNUM + 1
C	ROU = 0.0
CC	ROUT = 0.0
C-----   PART1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K)
	   Do K = 1, KBM
!$OMP DO
	      Do I = 1, IJM
	         If(CCM(I) .EQ. 1.0) Then
	           BBP(I,K) = 0.0
	           Do J = 1, CELL_POLYGEN(I)
	              If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	                 BBP(I,K) = BBP(I,K) + OMIGA * 
     &				  ASR(I,K,J) * R1(CELL_SIDE(I,J,2),K)
	              Endif
	           Enddo
	           BBP(I,K) = BBP(I,K) + OMIGA * 
     & 			   (ATR(I,K) * R1(I,K-1) + ABR(I,K) * R1(I,K+1))
	           BBP(I,K) = BBP(I,K) + (1. - OMIGA) * AP(I,K) * R1(I,K)
	         Endif
	      Enddo
!$OMP END DO
	   Enddo
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         Do K = 1, KBM
                  BBP(I,K) = BBP(I,K) + U1(I,K)
	         Enddo
	      Endif
	   Enddo
!$OMP END DO
!$OMP END PARALLEL	    
	   Call SOLVETRISOR
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         Do K = 1, KBM
	            V1(I,K) = XP(I,K)
	         Enddo
	      Endif
	   Enddo
!$OMP END PARALLEL DO
C-----   CONVERGENCE
	   AR = 0.0
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         Do K = 1, KBM
	            AR = Dmax1(AR,Abs(V1(I,K)-X(I,K)) / (1.+Abs(V1(I,K))))
	         Enddo
	      Endif
	   Enddo
	   PRINT*, AR, KNUM
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         Do K = 1, KBM
	            X(I,K) = V1(I,K)
	            R1(I,K) = X(I,K)
	         Enddo
	      Endif
	   Enddo
!$OMP END PARALLEL DO
	Enddo
	Return
	End
c===========================================================================c
c         preconditioning the matrix for SOR method                         c
c===========================================================================c
	Subroutine PREMATRIXSOR(OMIGA)
      Include './Include/OCERM_INF'
	Common/DYNBLK/AS(IJM,KB,IPOLYGEN),AB(IJM,KB),AT(IJM,KB),
     &	          AP(IJM,KB),BB(IJM,KB),X(IJM,KB),XINI(IJM,KB)
	Common/DYNPRE/ASL(IJM,KB,IPOLYGEN),ABL(IJM,KB),ATL(IJM,KB),
     &	          APL(IJM,KB),ASR(IJM,KB,IPOLYGEN),ABR(IJM,KB),
     &              ATR(IJM,KB),APR(IJM,KB),BBP(IJM,KB),XP(IJM,KB)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K)
	Do K = 1, KB
!$OMP DO	
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
!$OMP END DO
	Enddo
!$OMP END PARALLEL
C-----  L and LT
	GSCOE = Sqrt(OMIGA * (2. - OMIGA))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K)
	Do K = 1, KBM
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         APL(I,K) = Sqrt(AP(I,K)) / GSCOE
	         APR(I,K) = Sqrt(AP(I,K)) / GSCOE
	         Do J = 1, CELL_POLYGEN(I)
	            If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	               INDEXCELL = CELL_SIDE(I,J,2)
	               If(INDEXCELL .LT. I) Then
	                  ASL(I,K,J) = 
     &					  AS(I,K,J) * OMIGA / GSCOE / 
     &                      Sqrt(AP(INDEXCELL,K))
	               Endif
	            Endif
	         Enddo
	         If(Abs(AT(I,K)) .GT. 0.0) Then
		        ATL(I,K) = AT(I,K) * OMIGA / GSCOE / Sqrt(AP(I,K-1))
	            ABR(I,K-1) = ATL(I,K)
	         Endif
	      Endif
	   Enddo
!$OMP END DO NOWAIT
	Enddo
	Do K = 1, KBM
!$OMP DO	
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         Do J = 1, CELL_POLYGEN(I)
	            If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	               INDEXCELL = CELL_SIDE(I,J,2)
	               If(INDEXCELL .GT. I) Then
	                  ASR(I,K,J) = 
     &					  AS(I,K,J) * OMIGA / GSCOE / Sqrt(AP(I,K))
	               Endif
	            Endif
	         Enddo
	      Endif
	   Enddo
!$OMP END DO NOWAIT	   
	Enddo
!$OMP END PARALLEL	
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
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      P1(I,0) = 0.0
	      R1(I,0) = 0.0
	      Do K = 1, KB
	         R1(I,K) = 0.0
CCC	         V1(I,K) = 0.0
	         P1(I,K) = 0.0
	      Enddo
	   Endif
	Enddo
!$OMP END DO
!$OMP BARRIER
C-----   FORWARD
	Do K = 1, KBM
!$OMP DO	
         Do I = 1, IJM
		  If(CCM(I) .EQ. 1.0) Then
	        V1(I,K) = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
                   V1(I,K) = V1(I,K) + ASL(I,K,J)*P1(CELL_SIDE(I,J,2),K)
	           Endif
	        Enddo
	        V1(I,K)= V1(I,K) + ATL(I,K)*P1(I,K-1) + ABL(I,K)*P1(I,K+1)
	        V1(I,K)= BBP(I,K) + V1(I,K)
C	        P1(I,K) = V1(I,K) / APL(I,K)
	      Endif
	    Enddo
!$OMP END DO
	Enddo
	Do K = 1, KBM
!$OMP DO	
         Do I = 1, IJM
		  If(CCM(I) .EQ. 1.0) Then
C	        V1(I,K)= V1(I,K) + ATL(I,K)*P1(I,K-1) + ABL(I,K)*P1(I,K+1)
C	        V1(I,K)= BBP(I,K) + V1(I,K)
	        P1(I,K) = V1(I,K) / APL(I,K)
	     Endif
	   Enddo
!$OMP END DO
	Enddo
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