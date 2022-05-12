C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (16/02/2011)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                     computes the velocity                                 #
c                                                                           #
c############################################################################
      Subroutine PROFVIS(SIGMA)
      Include './Include/OCERM_INF'
	Common/VISBLK/VISSOURCE(IJM,KB),VISMIN(IJM,KB)
	Parameter(SCHEME = 2)
	Dimension AAAA(KBM,KBM), BBBB(KBM)
	Dimension AA(3*KBM-2), BB(KBM)
C===========================================================================C
C                       matrix of the equations                             c 
C===========================================================================C
c============================================================================c
c                    open boundary treatments                                c
c============================================================================c	   
C-----   elevation boundary condition
	   If(NUMEBC .NE. 0) Then
		 Do N = 1, NUMEBC
	        ID = IEBC(N)
	        Do K = 1, KBM
	           VIS(ID,K) = 0.0
	        Enddo
	     Enddo
	   Endif
C-----   astrotidal boundary condition
  	   If(NUMAST .NE. 0) Then
	     Do N = 1, NUMAST
	        ID = IABC(N)
	        Do K = 1, KBM
	           VIS(ID,K) = 0.0
	        Enddo
	     Enddo
	   Endif
c-----   discharge boundary condition
	   If(NUMQBC .NE. 0) Then
	      Do N = 1, NUMQBC
	         ID = IQBC(N)
	         Do K = 1, KBM
			    VIS(ID,K) = 0.0
	         Enddo
	      Enddo
	   Endif
C-----   velocity boundary condition
	   If(NUMVBC .NE. 0) Then
	      Do N = 1, NUMVBC
	         ID = IVBC(N)
	         Do K = 1, KBM
			    VIS(ID,K) = 0.0
	         Enddo
	      Enddo
	   Endif
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
C--------------------------------------------------------------------------C
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,K,K1,K2,AAAA,BBBB,AA,BB,L)
C!$OMP DO
	Do I = 1, IJM
	   IF(CCM(I) .EQ. 1.0) Then
C	     Do K1 = 1, KBM
C	        Do K2 = 1, KBM
C	           AAAA(K1,K2) = 0.0
C	        Enddo
C	     Enddo
	 	   Do K = 2, KBM 
	        AAAA(K,K-1) = -DTI * ((VIS(I,K-1)+VIS(I,K))/2. + UMOL) / 
     &			          DC(I) ** 2. / DZZ(K-1) / SIGMA
	        AAAA(K-1,K) = AAAA(K,K-1) 
	     Enddo
	     Do K = 2, KBM - 1
	        AAAA(K,K) = DZ(K) - AAAA(K,K-1) - AAAA(K,K+1)  
     &			          - DTI * VISSOURCE(I,K) 
	     Enddo
	     Do K = 1, KBM
	        BBBB(K) = VISF(I,K) / AREA(I) 
           Enddo
c----      Boundary condition   --------------------------------------------c      
	     If(KBM .GT. 1) Then
c	        AAAA(1,1) = DZ(1) - 2. * AAAA(1,2) - 
c     &		  DTI * (VIS(I,1)+UMOL) /      !  zero at the free surface 
c     &			        (.5 * DZ(1)) / DC(I) ** 2. / SIGMA
C-----   3MOL at the free surface     
c	        AAAA(1,1) = 1.0
c	        AAAA(1,2) = 0.0
C              BBBB(1) = 5. * UMOL
c	        BBBB(1) = 1.0*(AREA(I)*DC(I)*DZ(1))**2./3.*VORTICITY(I,1)
c	        BBBB(1) = (0.15*Dmin1(D2D(I),DC(I)))**2.*TENSOR(I,1)
c              AAAA(1,2) = 0.0
              AAAA(1,1) = DZ(1) -  2. * AAAA(1,2) - DTI * VISSOURCE(I,1)
              BBBB(1) = BBBB(1) - AAAA(1,2) *   
     &                  (0.1*Dmin1(D2D(I),DC(I)))**2. * TENSOR(I,1) 

ccc              AAAA(1,2) = 0.0
	     Else
	        AAAA(1,1) = DZ(1) +
     &		  DTI * (VIS(I,1)+UMOL) /      !  zero at the free surface
     &			        (.5 * DZ(1)) / DC(I) ** 2. / SIGMA 
     &                    -DTI * VISSOURCE(I,1)                
	     Endif
c-----     vis = 0.0 at the bottom  ----------------------------------------c
	     If(KBM .GT. 1) Then
	        AAAA(KBM,KBM) = DZ(KBM) - AAAA(KBM,KBM-1) + DTI * 
     &                      (UMOL+VIS(I,KBM)) / (.5*DZ(KBM)*DC(I))/
     &                      DC(I) / SIGMA - DTI * VISSOURCE(I,KBM)
	     Else
	        AAAA(KBM,KBM) = DZ(KBM) + DTI * 
     &                      (UMOL+VIS(I,KBM)) / (.5*DZ(KBM)*DC(I))/
     &                      DC(I) / SIGMA - DTI * VISSOURCE(I,KBM)
	     Endif
c---------------------------------------------------------------------------c
		   Goto (1,2) SCHEME
C===========================================================================C
C                  solve by operation of matrixs                            c
C===========================================================================C
1	     Continue
		   Call BRINV(AAAA,KBM,L)
		   Do K = 1, KBM
	        VIS(I,K) = 0.0
	        Do J = 1, KBM
	           VIS(I,K) = VIS(I,K) + AAAA(K,J) * BBBB(J) 
	        Enddo
	     Enddo
	     Goto 100
C===========================================================================C
C            solve by forward elimination and back-substitution             c
C===========================================================================C
2          Continue
C-----     Mass concentration  
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
	        VIS(I,K) = BB(K)
C	        IF(I .EQ. 30) PRINT*, VIS(I,K),K
	        VIS(I,K) = Dmax1(BB(K), 0.1*UMOL) 
C	        If(DES .NE. 'NEGLECT   ') Then
C	           D1 = DC(I) * (1. + ZZ(K)) 
C	           D2 = CDES * Sqrt(4. * AREA(I) / ACOS(-1.0))
C                 DTIDLE = Dmin1(D1,D2)
C                 If(DTIDLE .EQ. D2) 
C     &              VIS(I,K) = Dmax1(VIS(I,K),VISMIN(I,K))
C              Endif
	     Enddo
	     Goto 100
100	   Continue
	   Endif
	Enddo
!$OMP END PARALLEL DO
CC!$OMP END PARALLEL	
C---------------------------------------------------------------------------C
      Return
      End
