C############################################################################
C                     BI-CGSTAB method for equations solving                #
c                            VERSION 1.0 (25/05/2009)                       #
C                            AUTHORIZED BY ZHANG JINGXIN                    #
C                            SHANGHAI JIAO TONG UNIVERSITY                  #
C                                 SHANGHAI, CHINA                           #
c                                                                           #
c############################################################################
      Subroutine SOLVEELFPOLCG
      Include './Include/OCERM_INF'
	Parameter(EPSON=1.E-20,EPSI=1.E-6)
	Common/ELFBLK/CS(IJM,IPOLYGEN),CB(IJM),CP(IJM),X(IJM)
	Dimension R1(IJM),U1(IJM),V1(IJM),P1(IJM),R0(IJM),S1(IJM),T1(IJM)
	Dimension TRACE(IJM)
C===========================================================================C
C                optimazing the Matrix                                      c
c===========================================================================c
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J)
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
		  Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	            CS(I,J) = CS(I,J) / 
     &				   Sqrt(CP(I)) / Sqrt(CP(CELL_SIDE(I,J,2)))
	         Endif
		  Enddo
	      CB(I) = CB(I) / Sqrt(CP(I))
	   Endif
	Enddo
!$OMP END DO
!$OMP DO	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      TRACE(I) = CP(I)
	      CP(I) = 1.0
	   Endif
	Enddo
!$OMP END DO
!$OMP DO	
	Do I=1,IJM
	   If(CCM(I) .EQ. 1.0) Then
	      U1(I) = 0.0
	      V1(I) = 0.0
	      P1(I) = 0.0
	      S1(I) = 0.0
	      T1(I) = 0.0
	      R0(I) = 0.0
	      X(I) = 0.0
	   Endif
	Enddo
!$OMP END DO	
C-----   INITIAL VALUES
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      X(I)  = ELF(I) * Sqrt(TRACE(I))
	   Endif
	Enddo
!$OMP END DO
!$OMP BARRIER
!$OMP DO	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      R1(I) = 0.0
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	            R1(I) = R1(I) + 
     &			        CS(I,J) * X(CELL_SIDE(I,J,2))
	         Endif
            Enddo
	      R1(I) = CB(I) + R1(I) - CP(I) * X(I)
	      R0(I) = CB(I)
	      P1(I) = 0.0
	   Endif
	Enddo
!$OMP END DO	
!$OMP END PARALLEL
C---------------------------------------------------------------------------C
	ROU = 1.0
	ALPHA = 1.0
	OMEGA = 1.0
      K = 0
 10   Continue
      K = K+1
      print*, k
      BETA = ROU
	ROU = 0.0
CD    ROU=(B,R)
      Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     ROU = ROU + CB(I) * R1(I)
	   Endif
	Enddo
	If(ROU .EQ. 0.0) Then
	   BETA = 0.0
	Else
         BETA = ROU / (BETA) * ALPHA / OMEGA
	Endif
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     P1(I) = R1(I) + BETA * (P1(I) - OMEGA * V1(I))
	   Endif
	Enddo
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     V1(I) = 0.0
	     Do J = 1, CELL_POLYGEN(I)
	        If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	           V1(I) = V1(I) + CS(I,J) * P1(CELL_SIDE(I,J,2))
	        Endif
	     Enddo
	     V1(I) = -V1(I) + CP(I) * P1(I)
	   Endif
	Enddo
CD    ALPHA=ROU/(B,V)
	BV = 0.0
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      BV = BV + CB(I) * V1(I)
	   Endif
	Enddo
	If(ROU .EQ. 0.0) Then
	   ALPHA = 0.0
	Else
	   ALPHA = ROU / (BV)
	Endif
	Do I = 1,IJM
	   If(CCM(I) .EQ. 1.0) Then
	      R1(I) = R1(I) - ALPHA * V1(I)
	   Endif
	Enddo
	AR=0.0
	Do I = 1, IJM
   	   If(CCM(I) .EQ. 1.0) Then
	      AR = AR + R1(I) * R1(I)
	   Endif
      Enddo
	AR = Sqrt(AR)
c	ART = AR / AR0
      If(AR .LT. EPSI) Then
	   Do I = 1, IJM
     	      If(CCM(I) .EQ. 1.0) Then
	         X(I) = X(I) + ALPHA * P1(I)
	      Endif
	   Enddo
	   goto 1000
	Endif
	AR0 = AR
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      U1(I) = 0.0
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	            U1(I) = U1(I) + CS(I,J) * R1(CELL_SIDE(I,J,2))
	         Endif
	      Enddo
	      U1(I) = -U1(I) + CP(I) * R1(I)
	   Endif
	Enddo
CD    OMEGA=(U,R)/(U,U)
	UR1 = 0.0
	UU1 = 0.0
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      UR1 = UR1 + U1(I) * R1(I)
	      UU1 = UU1 + U1(I) * U1(I)
	   Endif
	Enddo
	If(UR1 .EQ. 0.0) Then
	   OMEGA = 0.0
	Else
	   OMEGA = UR1 / (UU1)
	Endif
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      X(I) = X(I) + ALPHA * P1(I) + OMEGA * R1(I)
	   Endif
	Enddo
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      R1(I) = R1(I) - OMEGA * U1(I)
	   Endif
	Enddo
	AR = 0.0
	Do I = 1, IJM
   	   If(CCM(I) .EQ. 1.0) Then
	      AR = AR + R1(I) * R1(I)
	   Endif
	Enddo
	AR = Sqrt(AR) 
c	ART = AR / AR0
C		PRINT*, AR,K
      If(AR .LT. EPSI .OR. K .GE. 200) Then
	   goto 1000
	Else
c	PRINT*, AR,K
	  Goto 10
	Endif
C-------------------------------------------------------------------------C
1000	Continue
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      X(I) = X(I) / Sqrt(TRACE(I))
C	      X(I) = Anint(X(I) * 1.E8)/1.E8
	   Endif
	Enddo
!$OMP END PARALLEL DO	
	Return
	End