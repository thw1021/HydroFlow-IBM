C############################################################################
C                     BI-CGSTAB method for equations solving                #
c                            VERSION 1.0 (28/07/2009)                       #
C                            AUTHORIZED BY ZHANG JINGXIN                    #
C                            SHANGHAI JIAO TONG UNIVERSITY                  #
C                                 SHANGHAI, CHINA                           #
c                                                                           #
c############################################################################
	Subroutine SOLVEDYN3DPOLBICG
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-20,EPSI=1.E-6,INUM=100)
C	Common/DYNBLK/AS(IJM,KB,IPOLYGEN),AB(IJM,KB),AT(IJM,KB),
C     &	          AP(IJM,KB),BB(IJM,KB),X(IJM,KB)
     
     	Common/DYNBLK/AS(IJM,KB,IPOLYGEN),AB(IJM,KB),AT(IJM,KB),
     &	          AP(IJM,KB),BB(IJM,KB),X(IJM,KB),XINI(IJM,KB),
     &            ASU(IJM,KB,IPOLYGEN),ASD(IJM,KB,IPOLYGEN),
     &            ABB(IJM,KB),ATT(IJM,KB)

	Dimension R1(IJM,KB),U1(IJM,KB),V1(IJM,KB),P1(IJM,0:KB)
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
	Do I=1,IJM
	   If(CCM(I) .EQ. 1.0) Then
	      P1(I,0) = 0.0
	      Do K = 1, KBM
C	         R1(I,K) = BB(I,K)
C	         X(I,K) = 0.0
	         U1(I,K) = 0.0
	         V1(I,K) = 0.0
	         P1(I,K) = 0.0
	      Enddo
	   Endif
	Enddo
C-----   INITIAL VALUES
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         X(I,K)  = PT(I,K) * Sqrt(TRACE(I,K))
	      Enddo
	   Endif
	Enddo
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
C---------------------------------------------------------------------------C
	ROU = 1.0
	ALPHA = 1.0
	OMEGA = 1.0
      KNUM = 0
 10   Continue
      KNUM = KNUM+1
      BETA = ROU
	ROU = 0.0
CD    ROU=(B,R)
      Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        ROU = ROU + BB(I,K) * R1(I,K)
	     Enddo
	   Endif
	Enddo
	BETA = ROU / (BETA + EPSON) * ALPHA / OMEGA
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        P1(I,K) = R1(I,K) + BETA * (P1(I,K) - OMEGA * V1(I,K))
	     Enddo
	   Endif
	Enddo
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
C	     Do K = 2, KBM - 1
C	        V1(I,K) = V1(I,K) + AT(I,K)*P1(I,K-1) + AB(I,K)*P1(I,K+1)
C	     Enddo
C	     V1(I,1) = V1(I,1) + AB(I,1) * P1(I,2)
C	     V1(I,KBM) = V1(I,KBM) + AT(I,KBM) * P1(I,KBM-1)
C	     Do K = 1, KBM
C		    V1(I,K) = -V1(I,K) + AP(I,K) * P1(I,K)
C	     Enddo
	   Endif
	Enddo
CD    ALPHA=ROU/(B,V)
	BV = 0.0
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        BV = BV + BB(I,K) * V1(I,K)
	     Enddo
	   Endif
	Enddo
	ALPHA = ROU / (BV + EPSON)
	Do I = 1,IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         R1(I,K) = R1(I,K) - ALPHA * V1(I,K)
	      Enddo
	   Endif
	Enddo
	AR=0.0
	Do I = 1, IJM
   	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         AR = AR + R1(I,K) * R1(I,K)
	      Enddo
	   Endif
      Enddo
	   AR = Sqrt(AR)
C	PRINT*, AR
      If(AR .LT. EPSI) Then
	  Do I = 1, IJM
     	     If(CCM(I) .EQ. 1.0) Then
	        Do K = 1, KBM
	           X(I,K) = X(I,K) + ALPHA * P1(I,K)
	        Enddo
	     Endif
	  Enddo
	  Goto 1000
	Endif
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        U1(I,K) = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	             U1(I,K) = U1(I,K) + AS(I,K,J)*R1(CELL_SIDE(I,J,2),K)
	           Endif
	        Enddo
	        U1(I,K) = U1(I,K) + AT(I,K)*R1(I,K-1) + AB(I,K)*R1(I,K+1)
	        U1(I,K) = -U1(I,K) + AP(I,K) * R1(I,K)
	     Enddo
C	     Do K = 2, KBM - 1
C	        U1(I,K) = U1(I,K) + AT(I,K)*R1(I,K-1) + AB(I,K)*R1(I,K+1)
C	     Enddo
C	     U1(I,1) = U1(I,1) + AB(I,1) * R1(I,2)
C	     U1(I,KBM) = U1(I,KBM) + AT(I,KBM) * R1(I,KBM-1)
C	     Do K = 1, KBM
C	        U1(I,K) = -U1(I,K) + AP(I,K) * R1(I,K)
C	     Enddo
	   Endif
	Enddo
CD    OMEGA=(U,R)/(U,U)
	UR1 = 0.0
	UU1 = 0.0
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         UR1 = UR1 + U1(I,K) * R1(I,K)
	         UU1 = UU1 + U1(I,K) * U1(I,K)
	      Enddo
	   Endif
	Enddo
	OMEGA = UR1 / (UU1 + EPSON)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         X(I,K) = X(I,K) + ALPHA * P1(I,K) + OMEGA * R1(I,K)
	      Enddo
	   Endif
	Enddo
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         R1(I,K) = R1(I,K) - OMEGA * U1(I,K)
	      Enddo
	   Endif
	Enddo
	AR = 0.0
	Do I = 1, IJM
   	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         AR = AR + R1(I,K) * R1(I,K)
	      Enddo
	   Endif
	Enddo
	AR = Sqrt(AR)
c       PRINT*, AR,KNUM
      If(AR .LT. EPSI .OR. KNUM .GE. INUM) Then
	  Goto 1000
	Else
C	PRINT*, AR,KNUM
C	PAUSE
	  Goto 10
	Endif
C-------------------------------------------------------------------------C
1000	Continue
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         X(I,K) = X(I,K) / Sqrt(TRACE(I,K))
	      Enddo
	   Endif
	Enddo
	Return
	End