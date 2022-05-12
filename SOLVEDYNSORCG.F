C############################################################################
C                     BI-CGSTAB method for equations solving                #
c                            VERSION 1.0 (28/07/2009)                       #
C                            AUTHORIZED BY ZHANG JINGXIN                    #
C                            SHANGHAI JIAO TONG UNIVERSITY                  #
C                                 SHANGHAI, CHINA                           #
c                                                                           #
c############################################################################
	Subroutine SOLVEDYNSORCG
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-20,EPSI=1.E-10,INUM=200,OMIGA = 1.25)
	Common/DYNBLK/AS(IJM,KB,IPOLYGEN),AB(IJM,KB),AT(IJM,KB),
     &	          AP(IJM,KB),BB(IJM,KB),X(IJM,KB),XINI(IJM,KB)
	Common/DYNPRE/ASL(IJM,KB,IPOLYGEN),ABL(IJM,KB),ATL(IJM,KB),
     &	          APL(IJM,KB),ASR(IJM,KB,IPOLYGEN),ABR(IJM,KB),
     &              ATR(IJM,KB),APR(IJM,KB),BBP(IJM,KB),XP(IJM,KB)
	Dimension R1(IJM,0:KB),U1(IJM,0:KB),V1(IJM,KB),P1(IJM,0:KB)
	Dimension TRACE(IJM,KB)
C===========================================================================C
C                Normalizing the Matrix                                     c
c===========================================================================c
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K)
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
		     Do J = 1, CELL_POLYGEN(I)
	            If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	              AS(I,K,J) = AS(I,K,J) / 
     &				   Sqrt(AP(I,K)) / Sqrt(AP(CELL_SIDE(I,J,2),K))
	            Endif
		     Enddo
	         If(K .GT. 1) Then
	            AT(I,K) = AT(I,K) / 
     &                 Sqrt(AP(I,K)) / Sqrt(AP(I,K-1))
	         Endif
	         If(K .LT. KBM) Then
	            AB(I,K) = AB(I,K) / 
     &                 Sqrt(AP(I,K)) / Sqrt(AP(I,K+1))
	         Endif
	         BB(I,K) = BB(I,K) / Sqrt(AP(I,K))
c	         AT(I,K) = AT(I,K) / AP(I,K)
c	         AB(I,K) = AB(I,K) / AP(I,K)
c	         BB(I,K) = BB(I,K) / AP(I,K)
c	         AP(I,K) = AP(I,K) / AP(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END DO NOWAIT
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         TRACE(I,K) = AP(I,K)
	         AP(I,K) = 1.0
	      Enddo
	   Endif
	Enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
C===========================================================================C
C                optimazing the Matrix                                      c
c===========================================================================c
	Call PREMATRIXSORCG(OMIGA)                     ! Preconditioning Matrix
c-----  INITIAL VALUES	
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,K,J)
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         X(I,K)  = XINI(I,K) * Sqrt(TRACE(I,K))
	      Enddo
	   Endif
	Enddo
!$OMP END DO 
!$OMP BARRIER
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        R1(I,K) = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	              R1(I,K) = R1(I,K) + 
     &				        AS(I,K,J) * X(CELL_SIDE(I,J,2),K)
	           Endif
             Enddo
             If(K .EQ. 1)
     &	        R1(I,K) = R1(I,K) + AB(I,K) * X(I,K+1)
             If(K .EQ. KBM)
     &	        R1(I,K) = R1(I,K) + AT(I,K) * X(I,K-1) 
             If(K .GT. 1 .AND. K .LT. KBM)
     &	        R1(I,K) = R1(I,K) + AT(I,K) * X(I,K-1) + 
     &			          AB(I,K) * X(I,K+1)
	       R1(I,K) = BB(I,K) + R1(I,K) - AP(I,K) * X(I,K)
	     Enddo
	   Endif
	Enddo
!$OMP END DO NOWAIT	
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         U1(I,K) = 0.0
	         V1(I,K) = 0.0
	         BBP(I,K) = R1(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END DO
!$OMP END PARALLEL
	Call SOLVETRISORCG
	ROU = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K) REDUCTION(+:ROU)	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         U1(I,K) = XP(I,K)
	         P1(I,K) = U1(I,K)
	         ROU = ROU + R1(I,K) * U1(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
c	ROUINI = ROU
C	CALL CPU_TIME(BEGIN1)
C---------------------------------------------------------------------------C
      KNUM = 0
	ALPHA = 1.0
	AR = 1.0E6
CCC	ROU = 1.0
 10   Continue
      Do While (AR .GT. EPSI .AND. KNUM .LT. INUM)
      KNUM = KNUM + 1
C-----   AP = (A,P)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K)
	Do K = 1, KBM
!$OMP DO	
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	        V1(I,K) = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	              V1(I,K) = V1(I,K) + AS(I,K,J)*P1(CELL_SIDE(I,J,2),K)
	           Endif
	        Enddo
	        V1(I,K) = V1(I,K) + AT(I,K)*P1(I,K-1) + AB(I,K)*P1(I,K+1)
	        V1(I,K) = -V1(I,K) + AP(I,K) * P1(I,K)
	      Endif
	   Enddo
!$OMP END DO 
	Enddo
!$OMP END PARALLEL
C-----    BV=(P,AP)
	BV = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K) REDUCTION(+:BV)	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        BV = BV + P1(I,K) * V1(I,K)
	     Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO	
C-----    ALPHA
	ALPHA = ROU / (BV + EPSON)
C-----    NEW X
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)
	Do I = 1, IJM
     	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         X(I,K) = X(I,K) + ALPHA * P1(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)	
	Do I = 1, IJM
     	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         R1(I,K) = R1(I,K) - ALPHA * V1(I,K)
	         BBP(I,K) = R1(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO	
C-----    PRE-MATRIX
	Call SOLVETRISORCG
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         U1(I,K) = XP(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO	
	ROUT = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K) REDUCTION(+:ROUT)	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         ROUT = ROUT + R1(I,K) * U1(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO	
	BETA = ROUT / (ROU + EPSON)
	ROU = ROUT
C-----   NEW P
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)
	Do I = 1, IJM
     	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         P1(I,K) = U1(I,K) + BETA * P1(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
	AR = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K) REDUCTION(+:AR)	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        AR = AR + R1(I,K) * R1(I,K)
C			Dmax1(AR,Abs(R1(I,K))/AP(I,K))
	     Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
	AR = Sqrt(AR)
	Enddo
c	print*, ar, knum
C      If(AR .LT. EPSI .OR. KNUM .GE. INUM) Goto 1000
C	Goto 10
C---------------------------------------------------------------------------C
C1000	Continue
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         X(I,K) = X(I,K) / Sqrt(TRACE(I,K))
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO	
c-----   End ---------------------------------------------------------------c
	Return
	End

c===========================================================================c
c         preconditioning the matrix for CG method                          c
c===========================================================================c
	Subroutine PREMATRIXSORCG(OMIGA)
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
C	               Do JJ = 1, CELL_POLYGEN(INDEXCELL)
C	                  If(CFM(CELL_SIDE(INDEXCELL,JJ,1)) .EQ. 1.0) Then
C	                     INDEXCELL2 = CELL_SIDE(INDEXCELL,JJ,2)
C						 If(INDEXCELL2 .EQ. I) Then
C	                        ASR(INDEXCELL,K,JJ) = 
C     &							AS(INDEXCELL,K,JJ) * OMIGA / GSCOE /
C     &                            Sqrt(AP(INDEXCELL,K))
C	                     Endif
C	                  Endif
C	               Enddo
	            Endif
	         Enddo
	         If(K .NE. 1) Then
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
	Subroutine SOLVETRISORCG
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