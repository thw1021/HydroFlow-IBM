C############################################################################
C                     BI-CGSTAB method for equations solving                #
c                            VERSION 1.0 (28/07/2009)                       #
C                            AUTHORIZED BY ZHANG JINGXIN                    #
C                            SHANGHAI JIAO TONG UNIVERSITY                  #
C                                 SHANGHAI, CHINA                           #
c                                                                           #
c############################################################################
	Subroutine SOLVEDYNPOLCG
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-10,EPSI=1.E-6,INUM=100)
	Common/DYNBLK/AS(IJM,KB,IPOLYGEN),AB(IJM,KB),AT(IJM,KB),
     &	          AP(IJM,KB),BB(IJM,KB),X(IJM,KB)
	Dimension R1(IJM,KB),U1(IJM,KB),V1(IJM,KB),P1(IJM,0:KB)
	Dimension TRACE(IJM,KB)
C===========================================================================C
C                optimazing the Matrix                                      c
c===========================================================================c
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K)
!$OMP DO
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
!$OMP END DO
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         TRACE(I,K) = AP(I,K)
	         AP(I,K) = 1.0
	      Enddo
	   Endif
	Enddo
!$OMP END DO
!$OMP DO
	Do I = 1, IJM
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
!$OMP END DO
C-----   INITIAL VALUES
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         X(I,K)  = PN(I,K) * Sqrt(TRACE(I,K))
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
	        P1(I,K) = R1(I,K)
	     Enddo
	   Endif
	Enddo
!$OMP END DO
!$OMP END PARALLEL
CCC	CALL CPU_TIME(begin1)
C---------------------------------------------------------------------------C
	ROU = 1.0
	ALPHA = 1.0
	OMEGA = 1.0
      KNUM = 0
 10   Continue
      KNUM = KNUM+1
      BETA = ROU
	ROU = 0.0
CD    ROU=(P,B)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K) REDUCTION(+:ROU)
      Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        ROU = ROU + R1(I,K) * R1(I,K)
	     Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
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
CD    BV=(P,AP)
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
	ALPHA = ROU / (BV + EPSON)
CD    X(N+1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K) 
      Do I = 1, IJM
     	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         X(I,K) = X(I,K) + ALPHA * P1(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
CD    R(N+1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        R1(I,K) = R1(I,K) - ALPHA * V1(I,K)
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
	     Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
      AR = Sqrt(AR)
CCC	print*, ar, knum
      If(AR .LT. EPSI .OR. KNUM .GE. INUM) Goto 1000
	UR1 = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K) REDUCTION(+:UR1)	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         UR1 = UR1 + R1(I,K) * R1(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO	
CD    BETA=(R,AP)/(P,AP)
	BETA = UR1 / (ROU + EPSON)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         P1(I,K) = R1(I,K) + BETA * P1(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO	
	Goto 10
C-------------------------------------------------------------------------C
1000	Continue
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         X(I,K) = X(I,K) / Sqrt(TRACE(I,K))
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
	Return
	End