C==============================================================================C
C               Gauss-Seidle iteration method for linear equations             c
C==============================================================================C
	Subroutine GSDL
	Include './Include/OCERM_INF'
	Parameter(EPSON = 1.E-2, ALF = 0.8)
	Common/DYNBLK/AS(IJM,KB,IPOLYGEN),AB(IJM,KB),AT(IJM,KB),
     &	          AP(IJM,KB),BB(IJM,KB),X(IJM,KB)
	Dimension R1(IJM,KB),U1(IJM,KB),V1(IJM,KB),P1(IJM,0:KB)
C===========================================================================C
C                optimazing the Matrix                                      c
c===========================================================================c
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
		     Do J = 1, CELL_POLYGEN(I)
	            AS(I,K,J) = AS(I,K,J) / AP(I,K)
		     Enddo
	         AT(I,K) = AT(I,K) / AP(I,K)
	         AB(I,K) = AB(I,K) / AP(I,K)
	         BB(I,K) = BB(I,K) / AP(I,K)
	         AP(I,K) = AP(I,K) / AP(I,K)
	      Enddo
	   Endif
	Enddo
C---------------------------------------------------------------------------C
	Do I=1,IJM
	   If(CCM(I) .NE. 0.0) Then
	      P1(I,0) = 0.0
	      Do K = 1, KBM
	         R1(I,K) = BB(I,K)
	         X(I,K) = 0.0
	         U1(I,K) = 0.0
	         V1(I,K) = 0.0
	         P1(I,K) = 0.0
	      Enddo
	   Endif
	Enddo
C----------------------------------------------------------------------------C
	KNUM = 0	
10	Continue
	KNUM = KNUM + 1
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        U1(I,K) = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	              U1(I,K) = U1(I,K) + AS(I,K,J)*P1(CELL_SIDE(I,J,2),K)
	           Endif
	        Enddo
	        U1(I,K) = U1(I,K) + AT(I,K)*P1(I,K-1) + AB(I,K)*P1(I,K+1)
	        U1(I,K) = BB(I,K) + U1(I,K)
	     Enddo
	   Endif
	Enddo
C-----   UPDATING
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         U1(I,K) = ALF * U1(I,K) + (1. - ALF) * P1(I,K)
	      Enddo
	   Endif
	Enddo
C-----   DECIDING WHETHER TO GO ON
	ERRORMAX = 1.E-20
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 1, KBM
	         ERROR = Abs(U1(I,K) - P1(I,K)) / (1. + Abs(U1(I,K)))
	         If(ERROR .GT. ERRORMAX) ERRORMAX = ERROR
	      Enddo
	   Endif
	Enddo
	PRINT*, KNUM, ERRORMAX
	If(ERRORMAX .LT. EPSON) Then
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         Do K = 1, KBM
	            X(I,K) = U1(I,K)
	         Enddo
	      Endif
	   Enddo
	   Goto 100
	Else
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         Do K = 1, KBM
	            P1(I,K) = U1(I,K)
	         Enddo
	      Endif
	   Enddo
	   Goto 10
	Endif
100	Continue
	Return
	End