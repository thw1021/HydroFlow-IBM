C############################################################################
C                     BI-CGSTAB method for equations solving                #
c                            VERSION 1.0 (28/07/2009)                       #
C                            AUTHORIZED BY ZHANG JINGXIN                    #
C                            SHANGHAI JIAO TONG UNIVERSITY                  #
C                                 SHANGHAI, CHINA                           #
c                                                                           #
c############################################################################
	Subroutine SOLVEDYNPOLCGCG
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-12,EPSI=1.E-6,INUM=50,OMIGA = 1.25)
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
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
		     Do J = 1, CELL_POLYGEN(I)
	            If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
c	            AS(I,K,J) = AS(I,K,J) / AP(I,K)
	              
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
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         TRACE(I,K) = AP(I,K)
	         AP(I,K) = 1.0
	      Enddo
	   Endif
	Enddo
C---------------------------------------------------------------------------C
	Call PREMATRIX2CG(OMIGA)                       ! Preconditioning Matrix
c-----  INITIAL VALUES	
	Do I=1,IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         R1(I,K) = BB(I,K)
	         X(I,K) = 0.0
	         U1(I,K) = 0.0
	         V1(I,K) = 0.0
ccc	         P1(I,K) = 0.0
	         BBP(I,K) = R1(I,K)
	      Enddo
	   Endif
	Enddo
	Call SOLVETRI2CG
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         U1(I,K) = XP(I,K)
	         P1(I,K) = U1(I,K)
	         ROU = ROU + R1(I,K) * U1(I,K)
	      Enddo
	   Endif
	Enddo
c	ROUINI = ROU
cc	CALL CPU_TIME(END1)
C---------------------------------------------------------------------------C
      KNUM = 0
	ALPHA = 1.0
CCC	ROU = 1.0
 10   Continue
      KNUM = KNUM+1
C	ROUT = ROU
C	   ROU = (R,S)
CC	ROU = 0.0
CC	Do I = 1, IJM
CC	   If(CCM(I) .EQ. 1.0) Then
CC	      Do K = 1, KBM
CC	         ROU = ROU + R1(I,K) * U1(I,K)
CC	      Enddo
CC	   Endif
CC	Enddo
C-----   AP = (A,P)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        V1(I,K) = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	              V1(I,K) = V1(I,K) + AS(I,K,J)*P1(CELL_SIDE(I,J,2),K)
	           Endif
	        Enddo
	        V1(I,K) = V1(I,K) + AT(I,K)*P1(I,K-1) + AB(I,K)*P1(I,K+1)
	        V1(I,K) = -V1(I,K) + AP(I,K) * P1(I,K)
	     Enddo
	   Endif
	Enddo
C-----    BV=(P,AP)
	BV = 0.0
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        BV = BV + P1(I,K) * V1(I,K)
	     Enddo
	   Endif
	Enddo
C-----    ALPHA
	ALPHA = ROU / (BV + EPSON)
C-----    NEW X
	Do I = 1, IJM
     	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         X(I,K) = X(I,K) + ALPHA * P1(I,K)
	         R1(I,K) = R1(I,K) - ALPHA * V1(I,K)
	         BBP(I,K) = R1(I,K)
	      Enddo
	   Endif
	Enddo
C-----    PRE-MATRIX
	Call SOLVETRI2CG
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         U1(I,K) = XP(I,K)
	      Enddo
	   Endif
	Enddo
	ROUT = 0.0
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         ROUT = ROUT + R1(I,K) * U1(I,K)
	      Enddo
	   Endif
	Enddo
	BETA = ROUT / (ROU + EPSON)
	ROU = ROUT
C-----   NEW P
	Do I = 1, IJM
     	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         P1(I,K) = U1(I,K) + BETA * P1(I,K)
	      Enddo
	   Endif
	Enddo
	AR = 0.0
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        AR = AR + Abs(R1(I,K))
	     Enddo
	   Endif
	Enddo
	AR = Sqrt(AR)
	print*, ar, knum
      If(AR .LT. EPSI .OR. KNUM .GE. INUM) Goto 1000
	Goto 10
C---------------------------------------------------------------------------C
1000	Continue
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         X(I,K) = X(I,K) / Sqrt(TRACE(I,K))
	      Enddo
	   Endif
	Enddo
C---------------------------------------------------------------------------C
	Return
	End

c===========================================================================c
c         preconditioning the matrix for CG method                          c
c===========================================================================c
	Subroutine PREMATRIX2CG(OMIGA)
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
C-----  L and LT
	GSCOE = Sqrt(OMIGA * (2. - OMIGA))
	Do K = 1, KBM
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
	               Do JJ = 1, CELL_POLYGEN(INDEXCELL)
	                  If(CFM(CELL_SIDE(INDEXCELL,JJ,1)) .EQ. 1.0) Then
	                     INDEXCELL2 = CELL_SIDE(INDEXCELL,JJ,2)
						 If(INDEXCELL2 .EQ. I) Then
	                        ASR(INDEXCELL,K,JJ) = 
     &							AS(INDEXCELL,K,JJ) * OMIGA / GSCOE /
     &                            Sqrt(AP(INDEXCELL,K))
	                     Endif
	                  Endif
	               Enddo
	            Endif
	         Enddo
	         If(Abs(AT(I,K)) .GT. 0.0) Then
		        ATL(I,K) = AT(I,K) * OMIGA / GSCOE / Sqrt(AP(I,K-1))
	            ABR(I,K-1) = ATL(I,K)
	         Endif
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
	Subroutine SOLVETRI2CG
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-12,EPSI=1.E-5,INUM=200)
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
                   V1(I,K) = V1(I,K) + ASL(I,K,J)*P1(CELL_SIDE(I,J,2),K)
	           Endif
	        Enddo
	        V1(I,K)= V1(I,K) + ATL(I,K)*P1(I,K-1) + ABL(I,K)*P1(I,K+1)
	        V1(I,K)= BBP(I,K) + V1(I,K)
	        P1(I,K) = V1(I,K) / APL(I,K)
	     Endif
	   Enddo
	Enddo
C-----   BACKWARD
C	R1(IJM,KBM) = P1(IJM,KBM) / P1(IJM,KBM)
	Do K = KBM, 1, -1
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
	        R1(I,K) = V1(I,K) / APR(I,K)
	      Endif
	   Enddo
	Enddo
	Do K = 1, KBM
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         XP(I,K) = R1(I,K)
	      Endif
	   Enddo
	Enddo
	Return
	End