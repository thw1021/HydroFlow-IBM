C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (02/08/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c      wave generation by means of adding source terms into the             #
c      governing equations                                                  #
c                                                                           #
c############################################################################

	Subroutine WAVEGEN(KKK)
      Include './Include/OCERM_INF'
	Parameter (IWT = 1,I_D = 1)
	Dimension FRE_IRR(101), PHA(100), FRE_CHAR(100),WAVE_K(100) 
c	Dimension IGEN(N_SOURCE)
c===========================================================================c
c                 define cells as sources for wave generation               c
c===========================================================================c
	K = 0
	Do I = 1, IJM
	   If(CXY(I,1) .GT. 0.0 .AND. CXY(I,1) .LT. 0.005) Then
	      K = K + 1
		  IGEN(K) = I
	   Endif
	Enddo
c===========================================================================c
c                 wave parameters                                           c
c===========================================================================c
	WAVENUM = 2.*Acos(-1.0) / 0.01
	WAVEA = 0.003
	HGEN = 0.05
	WAVEFRE = Sqrt(9.8*WAVENUM*Tanh(WAVENUM*HGEN))
	WPEROID = 2.*Acos(-1.0) / WAVEFRE
	W_DIR_X = 1.0
	W_DIR_Y = 0.0
c===========================================================================c
c                 calculation of sources                                    c
c===========================================================================c
	If(KKK.EQ.0) Then
	   WTIME = (NSTEP - 1) * DTI
	Else
	   WTIME=NSTEP * DTI
	Endif
	Do I = 1, N_SOURCE
	   II = IGEN(I)
	   Do K = 1, KBM
	      WGEN(I,K) = 0.0
            If(IWT .EQ. 1) Then
	         Do J = 1, CELL_POLYGEN(II)
C-----     Shallow water
                  If(I_D .EQ. 1) Then
	  	        WGEN(I,K) = WGEN(I,K) + 
     *     			 Sqrt(9.8 * Tanh(WAVENUM * HC(II)) / WAVENUM) * 
     *			     WAVEA * Sin(WAVEFRE * WTIME) * CELL_CUV(II,J,6) *
     *                 Abs(CELL_CUV(II,J,7) * W_DIR_X +
     &                     CELL_CUV(II,J,8) * W_DIR_Y) /
     &                 AREA(II) 
                  Endif     
c-----     Middle water depth
                  If(I_D .EQ. 2) Then
	  	        WGEN(I,K) = WGEN(I,K) + 
     *     			 WAVEFRE * Cosh(WAVENUM*(1.+ZZ(K))*HC(II)) /
     &                         Sinh(WAVENUM*HC(II)) *  
     *			     WAVEA * Sin(WAVEFRE * WTIME) * CELL_CUV(II,J,6) *
     *                 Abs(CELL_CUV(II,J,7) * W_DIR_X +
     &                     CELL_CUV(II,J,8) * W_DIR_Y) /
     &                 AREA(II) * HC(II)
                  Endif
c-----     Deep water depth  
                  If(I_D .EQ. 3) Then
	  	        WGEN(I,K) = WGEN(I,K) + 
     *     			 WAVEFRE * Exp(WAVENUM*ZZ(K) * HC(II)) *  
     *			     WAVEA * Sin(WAVEFRE * WTIME) * CELL_CUV(II,J,6) *
     *                 Abs(CELL_CUV(II,J,7) * W_DIR_X +
     &                     CELL_CUV(II,J,8) * W_DIR_Y) /
     &                 AREA(II) * HC(II)
                  Endif
	         Enddo
	      Endif
            If(IWT.EQ.2) Then
               Do J = 1, CELL_POLYGEN(II)
		         WGEN(I,K) = WGEN(I,K) +
     *		          Sqrt(9.8 * Tanh(WAVENUM * HC(II)) / WAVENUM) *
     *		         (WAVEA * WAVENUM *
     *		     Cosh(WAVENUM*(1.+ZZ(K))*HC(II)) /
     *            Sinh(WAVENUM * HC(II)) * Cos(WAVEFRE * WTIME) +
     *            3.0 / 4.0 * WAVEA ** 2.0 * WAVENUM ** 2. *
     *            Cosh(2.0 * WAVENUM * (1.+ZZ(K)) * HC(II)) /
     *            (Sinh(WAVENUM * HC(II))) ** 4.0 * 
     *            Cos(2 * WAVEFRE * WTIME)) *  CELL_CUV(II,J,6) *
     *                 Abs(CELL_CUV(II,J,7) * W_DIR_X +
     &                     CELL_CUV(II,J,8) * W_DIR_Y) /
     &                 AREA(II) * HC(II)
	         Enddo
	      Endif
            If(IWT.EQ.3) Then
		       
		       WVEL= Sqrt(9.8 * Tanh(WAVENUM * HC(II)) / WAVENUM *
     &	                 (1. + (WAVENUM * WAVEA) ** 2. *
     &                    (8. + Cosh(4. * WAVENUM * HC(II))) /
     &                    8. / (Cosh(WAVENUM * HC(II))) ** 4.)) 
	         F1  = WAVENUM * WAVEA / Sinh(WAVENUM * HC(II))
	         F2  = 3. * (WAVENUM * WAVEA) ** 2. / 4. /
     &	            (Sinh(WAVENUM * HC(II))) ** 4.
	         F3  = 3. * (WAVENUM * WAVEA) ** 3. *
     &			    (11. - 2. * Cosh(2. * WAVENUM * HC(II))) / 64. /
     &              (Sinh(WAVENUM * HC(II))) ** 7.

               Do J = 1, CELL_POLYGEN(II)
		         WGEN(I,K) = WGEN(I,K) + WVEL * (
     &			             F1*Cosh(WAVENUM * (1. + ZZ(K)) * HC(II)) *
     &                       Cos(WAVEFRE * WTIME) +
     &                       F2*Cosh(2. * WAVENUM * (1.+ZZ(K))*HC(II)) *
     &                       Cos(2. * WAVEFRE * WTIME) +
     &                       F3*Cosh(3. * WAVENUM * (1.+ZZ(K))*HC(II)) *
     &                       Cos(3. * WAVEFRE * WTIME)) *
     &                       CELL_CUV(II,J,6) *
     *                       Abs(CELL_CUV(II,J,7) * W_DIR_X +
     &                           CELL_CUV(II,J,8) * W_DIR_Y) / AREA(II) 
     &                       * HC(II)
	         Enddo
	      Endif
c-------  Cnoidal wave generation  ----------------------------------------            
		  If(IWT .EQ. 5) Then
		       WVEL= Sqrt(9.8 * Tanh(WAVENUM * HC(II)) / WAVENUM *
     &	                 (1. + (WAVENUM * WAVEA) ** 2. *
     &                    (8. + Cosh(4. * WAVENUM * HC(II))) /
     &                    8. / (Cosh(WAVENUM * HC(II))) ** 4.)) 
	         F1  = WAVENUM * WAVEA / Sinh(WAVENUM * HC(II))
	         F2  = 3. * (WAVENUM * WAVEA) ** 2. / 4. /
     &	            (Sinh(WAVENUM * HC(II))) ** 4.
	         F3  = 3. * (WAVENUM * WAVEA) ** 3. *
     &			    (11. - 2. * Cosh(2. * WAVENUM * HC(II))) / 64. /
     &              (Sinh(WAVENUM * HC(II))) ** 7.

               Do J = 1, CELL_POLYGEN(II)
		         WGEN(I,K) = WGEN(I,K) + WVEL * (
     &			             F1*Cosh(WAVENUM * (1. + ZZ(K)) * HC(II)) *
     &                       Cos(WAVEFRE * WTIME) +
     &                       F2*Cosh(2. * WAVENUM * (1.+ZZ(K))*HC(II)) *
     &                       Cos(2. * WAVEFRE * WTIME) +
     &                       F3*Cosh(3. * WAVENUM * (1.+ZZ(K))*HC(II)) *
     &                       Cos(3. * WAVEFRE * WTIME)) *
     &                       CELL_CUV(II,J,6) *
     *                       Abs(CELL_CUV(II,J,7) * W_DIR_X +
     &                           CELL_CUV(II,J,8) * W_DIR_Y) / AREA(II) 
     &                       * HC(II)
	         Enddo
	    Endif
c-------  irregular wave generation  ----------------------------------------            
		  If(IWT.EQ.6) Then
	        FRE_M = 20
	        FRE_L = 3.
	        FRE_R = 30.
	        FRE_IRR(1) = 3.
	        GAMA = 3.3
	        
	        COE1=0.0624 / (0.23 + 0.0336 * GAMA - 0.185 / (1.9+GAMA))
	        
	        Do I6 = 2, 101
	           FRE_IRR(I6) = FRE_IRR(I6-1) + (FRE_R - FRE_L) / 100.
	        Enddo
		      Call RANDOM_NUMBER(FRE_CHAR)
	        Do I6 = 1, 100
	           FRE_CHAR(I6) =
     &			           ((1. - FRE_CHAR(I6)) * FRE_IRR(I6) +
     &                        FRE_CHAR(I6) * FRE_IRR(I6+1))
	        Enddo
		      Call RANDOM_NUMBER(PHA)
	        Do I6 = 1, 100
	           PHA(I6) = PHA(I6) * 2. * 3.1415926
	        Enddo
	       	Do M = 1, 100
	           X1 = 1.
	           X2 = 600.
	           F1 = FRE_CHAR(M) **2. / 9.8 / X1 - Tanh(X1*HC(II))
	           F2 = FRE_CHAR(M) **2. / 9.8 / X2 - Tanh(x2*HC(II))
	           If(F1 * F2 .GT. 0.0) Then
	             Print*, 'error, input new initial datas !'
	           Stop
                Endif
100	           Continue
	           XNEW = (X1 + X2) / 2.
	           FNEW = FRE_CHAR(M) **2. / 9.8 / XNEW - Tanh(XNEW*HC(II))
	           If(Abs(FNEW) .LE. 1.e-5) Then
	             WAVE_K(M) = XNEW
	           Else
	             If(FNEW * F1 .LT. 0.0) X2 = XNEW
	             If(FNEW * F2 .LT. 0.0) X1 = XNEW
	             Goto 100
	           Endif
	        Enddo
	        
	        WGEN(I,K)=0.0
	        Do K123 = 1,100
	           If(FRE_CHAR(K123) .LE. FRE_M) Then
	              SIGMA = 0.07
	           Else
	              SIGMA = 0.09
	           Endif
			   FRE_MID = .5 * (FRE_IRR(K123 + 1) + FRE_IRR(K123))
			   WGEN(I,K) = WGEN(I,K) +
     *			    Sqrt(2. * (FRE_IRR(K123+1) - FRE_IRR(K123))*
     *                     COE1 * 9.8 **2. / FRE_MID ** 5. * 
     *                     Exp(-1.25 * (FRE_M / FRE_MID) **4.) *
     *                     3.3 ** (Exp(-(FRE_MID - FRE_M) **2. /
     *                     2. / SIGMA ** 2. / FRE_M ** 2.))) * 
     *                     FRE_CHAR(K123) * 
     *                     Cosh(WAVE_K(K123) * (1.+ZZ(K)) * HC(II))/
     *                     Sinh(WAVE_K(K123) * HC(II)) *
     *                     Cos(FRE_CHAR(K123) * WTIME + PHA(K123))
	         Enddo
	         WGEN(I,K) = 2. * WGEN(I,K) / Sqrt(HC(II))
	      Endif	      
c----------------------------------------------------------------------------
		    If(WTIME .LE. 3.0 * WPEROID) 
     *         WGEN(II,K) = WGEN(II,K) *
     *                      (1.0 - Exp(-2.0 * WTIME / WPEROID))
	   Enddo
	Enddo
C	   WRITE(IUTSR,'(100F12.4)') WTIME, (WGEN(1,K), K=1, KBM)
	Return
	End