C############################################################################
C                     BI-CGSTAB method for equations solving                #
c                            VERSION 1.0 (28/07/2009)                       #
C                            AUTHORIZED BY ZHANG JINGXIN                    #
C                            SHANGHAI JIAO TONG UNIVERSITY                  #
C                                 SHANGHAI, CHINA                           #
c                                                                           #
c############################################################################
	Subroutine SOLVEDYNICBICG
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-20,EPSI=1.E-6,INUM=100)
	Common/DYNBLK/AS(IJM,KB,IPOLYGEN),AB(IJM,KB),AT(IJM,KB),
     &	          AP(IJM,KB),BB(IJM,KB),X(IJM,KB),XINI(IJM,KB)
	Dimension R1(IJM,0:KB),U1(IJM,0:KB),V1(IJM,0:KB),P1(IJM,0:KB),
     &            Z1(IJM,0:KB),Y1(IJM,0:KB),S1(IJM,0:KB)	
	Dimension ASC(IJM,KB,IPOLYGEN),ABC(IJM,KB),ATC(IJM,KB),
     &	          APC(IJM,KB),ASCT(IJM,KB,IPOLYGEN),ABCT(IJM,KB),
     &            ATCT(IJM,KB),APCT(IJM,KB)
     	Dimension TRACE(IJM,KB)
      IIII = 0
      IJM_B = IJM_DYN_B + IIII
      IJM_E = IJM_DYN_E + IIII

C===========================================================================C
C                Preconditioning Matrix A                                   c
c---------------------------------------------------------------------------c
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K) 
	Do K = 1, KB
!$OMP DO
	   Do I = IJM_B, IJM_E
	      Do J = 1, CELL_POLYGEN(I)
	         ASC(I,K,J) = 0.0
	         ASCT(I,K,J) = 0.0
	      Enddo
		    ABC(I,K) = 0.0
		    ATC(I,K) = 0.0
		    APC(I,K) = 0.0
		    ABCT(I,K) = 0.0
		    ATCT(I,K) = 0.0
		    APCT(I,K) = 0.0
	   Enddo
!$OMP END DO	   
	Enddo
!$OMP END PARALLEL
	
C	APC(1,1) = Sqrt(AP(1,1))
C	APCT(1,1) = Sqrt(AP(1,1))

      Do I = IJM_B, IJM_E
         If(CCM(I) .EQ. 1.0) Then
            Do K = 1, KBM
               Do J = 1, CELL_POLYGEN(I)
                  IL = I
	            IR = CELL_SIDE(I,J,2)
	            If(IR .LT. I .AND. IR .NE. -999)
     &            APC(I,K) = APC(I,K) + (AS(I,K,J) / APC(IR,K)) ** 2.
               Enddo
               If(K .GT. 1)
     &              APC(I,K) = APC(I,K) + (AT(I,K) / APC(I,K-1)) ** 2.
               If(AP(I,K) .GT. APC(I,K)) Then
                  APC(I,K) = Sqrt(AP(I,K) - APC(I,K))
                  APCT(I,K) = APC(I,K)
               Else
                 APC(I,K) = Sqrt(AP(I,K))
                 APCT(I,K) = APC(I,K)
               Endif 
            Enddo
         Endif
      Enddo
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K) 
!$OMP DO      
      Do I = IJM_B, IJM_E
         If(CCM(I) .EQ. 1.0) Then
            Do K = 1, KBM
               Do J = 1, CELL_POLYGEN(I)
                  IL = I
	            IR = CELL_SIDE(I,J,2)
	            If(IR .LT. I .AND. IR .NE. -999) Then
                     ASC(I,K,J) = AS(I,K,J) / APC(IR,K) 
                  Else
                     ASCT(I,K,J) = AS(I,K,J) / APCT(I,K)
                  Endif  
               Enddo
               If(K .GT. 1)
     &            ATC(I,K) =  AT(I,K) / APC(I,K-1)
               If(K .LT. KBM)
     &            ABCT(I,K) = AB(I,K) / APCT(I,K)
            Enddo
         Endif   
      Enddo
!$OMP END DO
!$OMP END PARALLEL   
C---------------------------------------------------------------------------C
C           ENDING MATRIX PREDICTIONING                                     C   
C===========================================================================C
C-----   INITIAL VALUES
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K) 
!$OMP DO
	Do I=1,IJM
	   If(CCM(I) .EQ. 1.0) Then
	      P1(I,0) = 0.0
	      R1(I,0) = 0.0
	      U1(I,0) = 0.0
	      V1(I,0) = 0.0
	      Z1(I,0) = 0.0
	      Y1(I,0) = 0.0
	      S1(I,0) = 0.0
	      Do K = 1, KB
	         R1(I,K) = 0.0
C	         X(I,K) = 0.0
	         U1(I,K) = 0.0
	         V1(I,K) = 0.0
	         P1(I,K) = 0.0
	         Z1(I,K) = 0.0
	         Y1(I,K) = 0.0
	         S1(I,K) = 0.0
	         X(I,K) = 0.0
	      Enddo
	   Endif
	Enddo
!$OMP END DO
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
	ROU = 1.0
	ALPHA = 1.0
	OMEGA = 1.0
      KNUM = 0
C===========================================================================C      
C----      ITERATION PROGRESS  ---------------------------------------------C      
 10   Continue
      KNUM = KNUM+1
      BETA = ROU
C---  Y = INV(M) * P      
      Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         U1(I,K) = 0.0
	         Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	              U1(I,K) = U1(I,K) + ASC(I,K,J)*Y1(CELL_SIDE(I,J,2),K) 
	           Endif
	         Enddo
	         U1(I,K) = U1(I,K) + ATC(I,K) * Y1(I,K-1)                   
	         Y1(I,K) = (P1(I,K) + U1(I,K)) / APC(I,K)                   
	      Enddo
	   Endif
	Enddo   
	Do I = IJM, 1, -1
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = KBM, 1, -1
	         U1(I,K) = 0.0
	         Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	              U1(I,K) = U1(I,K)+ASCT(I,K,J)*Y1(CELL_SIDE(I,J,2),K)  
	           Endif
	         Enddo
	         U1(I,K) = U1(I,K) + ABCT(I,K) * Y1(I,K+1)                  
	         Y1(I,K) = (Y1(I,K) + U1(I,K)) / APCT(I,K)                  
	      Enddo
	   Endif
	Enddo
C---  Q*Y	
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,K)
      Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        V1(I,K) = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	              V1(I,K) = V1(I,K) + AS(I,K,J)*Y1(CELL_SIDE(I,J,2),K)  
	           Endif
	        Enddo
	        V1(I,K) = V1(I,K) + AT(I,K)*Y1(I,K-1) + AB(I,K)*Y1(I,K+1)   
	        V1(I,K) = -V1(I,K) + AP(I,K) * Y1(I,K)                      
	     Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
C      BETA = ROU
C	ROU = 0.0
CD    ROU=(B,R)
C      Do I = 1, IJM
C	   If(CCM(I) .EQ. 1.0) Then
C	     Do K = 1, KBM
C	        ROU = ROU + BB(I,K) * R1(I,K)                              
C	     Enddo
C	   Endif
C	Enddo
	BV = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K) REDUCTION(+:BV)	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        BV = BV + BB(I,K) * V1(I,K)
	     Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
C---  ALPHA = ROU / R0 * QY	
	ALPHA = ROU / (BV + EPSON)                                          
C---  S(K) = R(K) - ALPHA * QY
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)
	Do I = 1,IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         S1(I,K) = R1(I,K) - ALPHA * V1(I,K)                        
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
C---  Z = INV(M) * S
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         U1(I,K) = 0.0
	         Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	              U1(I,K) = U1(I,K) + ASC(I,K,J)*Z1(CELL_SIDE(I,J,2),K) 
	           Endif
	         Enddo
	         U1(I,K) = U1(I,K) + ATC(I,K) * Z1(I,K-1)                   
	         Z1(I,K) = (S1(I,K) + U1(I,K)) / APC(I,K)                   
	      Enddo
	   Endif
	Enddo   
	Do I = IJM, 1, -1
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = KBM, 1, -1
	         U1(I,K) = 0.0
	         Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	              U1(I,K) = U1(I,K)+ASCT(I,K,J)*Z1(CELL_SIDE(I,J,2),K)  
	           Endif
	         Enddo
	         U1(I,K) = U1(I,K) + ABCT(I,K) * Z1(I,K+1)                  
	         Z1(I,K) = (Z1(I,K) + U1(I,K)) / APCT(I,K)                  
	      Enddo
	   Endif
	Enddo
C---  Q*Z
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        U1(I,K) = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	             U1(I,K) = U1(I,K) + AS(I,K,J)*Z1(CELL_SIDE(I,J,2),K)  
	           Endif
	        Enddo
	        U1(I,K) = U1(I,K) + AT(I,K)*Z1(I,K-1) + AB(I,K)*Z1(I,K+1)  
	        U1(I,K) = -U1(I,K) + AP(I,K) * Z1(I,K)                     
	     Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
C---  OMEGA = (QZ,S) / (QZ,QZ)
	UR1 = 0.0
	UU1 = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K) REDUCTION(+:UR1,UU1)	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         UR1 = UR1 + U1(I,K) * S1(I,K)
	         UU1 = UU1 + U1(I,K) * U1(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
	OMEGA = UR1 / (UU1 + EPSON)                                        
C---  X(K+1) = X(K) + ALPHA * Y + OMEGA * Z
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,K)
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         X(I,K) = X(I,K) + ALPHA * Y1(I,K) + OMEGA * Z1(I,K)       
	      Enddo
	   Endif
	Enddo
!$OMP END DO	
C---  R(K+1) = S - OMEGA * QZ
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         R1(I,K) = S1(I,K) - OMEGA * U1(I,K)                       
	      Enddo
	   Endif
	Enddo
!$OMP END DO
!$OMP END PARALLEL	
C---  BETA = ALPHA / OMEGA * (R0,R(K+1)) / (R0,R(K))
      ROU = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K) REDUCTION(+:ROU)	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         ROU = ROU + BB(I,K) * R1(I,K)
	      Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
	BETA = ROU / (BETA + EPSON) * ALPHA / (OMEGA + EPSON)              
C---  P(K+1) = R(K+1) + BETA * (P(K) - OMEGA * QZ)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        P1(I,K) = R1(I,K) + BETA * (P1(I,K) - OMEGA * V1(I,K))      
	     Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
	AR = 0.0
	B1 = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K) REDUCTION(+:AR,B1)	
	Do I = IJM_B, IJM_E
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        AR = AR + R1(I,K) * R1(I,K)
	        B1 = B1 + BB(I,K) * BB(I,K)
	     Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
      AR = Sqrt(AR)
      B1 = Sqrt(B1)
      AR = AR / B1
C	print*, ar, knum
      If(AR .LT. EPSI .OR. KNUM .GE. INUM) Goto 1000
      Goto 10
C-------------------------------------------------------------------------C
1000	Continue
	Return
	End