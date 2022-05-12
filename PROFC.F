C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (25/05/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                     computes the velocity                                 #
c                                                                           #
c############################################################################
      Subroutine PROFC
      Include './Include/OCERM_INF'
	Parameter(SCHEME = 2,PRNUM=0.78)
	Dimension AAAA(KBM,KBM), BBBB(KBM)
	Dimension AA(3*KBM-2), BB(KBM)
c===========================================================================c
C              initializing the arrays                                      c
c===========================================================================c
	Do K1 = 1, KBM
         Do K2 = 1, KBM
	      AAAA(K1,K2) = 0.0
	   Enddo
	Enddo
	Do K1 = 1, 3 * KBM - 2
	   AA(K1) = 0.0
	Enddo
	Do K1 = 1, KBM
	   BB(K1) = 0.0
	   BBBB(K1) = 0.0
	Enddo
c===========================================================================c
c                  arrays for the turbulent kinetic governing equation      c
c===========================================================================c
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,K,K1,K2,AAAA,BBBB,AA,BB,L)
	Do I = 1, IJM
	   IF(CCM(I) .EQ. 1.0) Then
	 	   Do K = 2, KBM 
	        AAAA(K,K-1) = -DTI * PRNUM* (UMOL+(KM(I,K)+KM(I,K-1)) / 2.) / 
     &			            DC(I) ** 2. / DZZ(K-1)
	        AAAA(K-1,K) = AAAA(K,K-1) 
	     Enddo
	     Do K = 2, KBM - 1
	        AAAA(K,K) = DZ(K) - AAAA(K,K-1) - AAAA(K,K+1)  
	     Enddo
	     Do K = 1, KBM
	        BBBB(K) = CF(I,K) / AREA(I)
		   Enddo
c----      Boundary condition   --------------------------------------------c      
	     If(KBM .GT. 1) Then
	        AAAA(1,1) = DZ(1) - AAAA(1,2)
	     Else
	        AAAA(1,1) = DZ(1)
	     Endif
c-----     at the bottom  ----------------------------------------c
	     If(KBM .GT. 1) Then
	        AAAA(KBM,KBM) = DZ(KBM) - AAAA(KBM,KBM-1) 
	     Else
	        AAAA(KBM,KBM) = DZ(KBM) 
	     Endif
C===========================================================================C
C            solve by forward elimination and back-substitution             c
C===========================================================================C
	     AA(1) = AAAA(1,1)
	     If(KBM .GT. 1) AA(2) = AAAA(1,2)
		   Do K = 2, KBM - 1
		      AA(2*(K-1)+K-1) = AAAA(K,K-1)
	        AA(2*(K-1)+K)   = AAAA(K,K)
		      AA(2*(K-1)+K+1) = AAAA(K,K+1)
		   Enddo 
	     If(KBM .GT. 1) AA(3*KBM-3) = AAAA(KBM,KBM-1)
	     AA(3*KBM-2) = AAAA(KBM,KBM)
		   Do K = 1, KBM
	        BB(K) = BBBB(K)
	     Enddo
           Call ATRDE(AA,KBM,3*KBM-2,BB,L)
	     Do K = 1, KBM
	        C(I,K) = BB(K)
	     Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
C---------------------------------------------------------------------------C
      Return
      End
