C###########################################################################C
c                 SUBROUNTINE PROGRAM                                       #
C                 SOLID BOUNDARY CONDITION BASED ON DES                     #
c                               VERSION 1.0 (16/02/2011)                    #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C############################################################################
	Subroutine DESMODULE
      Include './Include/OCERM_INF'
	Parameter(CB1=0.1355, CB2=0.622, VK=0.41, SIGMA=2./3.)
      Parameter(CW2=0.3, CW3=2.0, CV1=7.1, CDES=0.65, CT3=1.2, CT4=0.5)
	Dimension SGS(IJM,KB)
c===========================================================================c
c              transport equation of working variable                       c
c===========================================================================c
	  CW1 = CB1 / VK ** 2. + (1. + CB2) / SIGMA 
        Call ADVVIS(CB1,CB2,CW1,CW2,CW3,SIGMA,VK,CV1,CDES,CT3,CT4)
	  Call PROFVIS(SIGMA)
c===========================================================================c
c               natural eddy viscosity                                      c
c===========================================================================c
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(COELES,ZSTAR,XX,FV1,D1,D2,
!$OMP&            DTIDLE,SLES)
!$OMP DO
	  Do I = 1, IJM
	     If(CCM(I) .EQ. 1.0) Then
              Do K = 1, KBM
	           COELES = 0.0
		         ZSTAR = DC(I) * (1. + ZZ(K)) / 1.E-6 * 
     &			       Sqrt(Sqrt(TBX(I)**2.+TBY(I)**2.)/RMEAN(I,KBM)) 
c	           VIS(I,K) = VIS(I,K) * Sqrt(1. - Exp(-(ZSTAR/25.)**3.))
  			     XX = VIS(I,K) / UMOL
CC	           XX = (VIS(I,K) + UMOL) / UMOL + 0.5 * 30. * Z01(I) / 
CC     &			    (DC(I) * (1. + ZZ(K)) + Z01(I))
			     FV1 = XX ** 3. / (XX ** 3. + CV1 ** 3.)
C	           D1 = DC(I) * (1. + ZZ(K)) 
C	           D2 = CDES * Sqrt(4. * AREA(I) / ACOS(-1.0))
C                DTIDLE = Dmin1(D1,D2)
C	           If(DTIDLE .EQ. D2) COELES = 1.0
C	           SLES = (0.1*4.* AREA(I)/ACOS(-1.0))**2. * TKE(I,K)
C			     SGS(I,K) = Dmax1(VIS(I,K)*FV1, COELES*SLES)
				 SGS(I,K) = VIS(I,K) * FV1
CCC     &			          * Sqrt(1. - Exp(-(ZSTAR/25.)**3.))
	        Enddo
           Endif
	  Enddo
!$OMP END DO	  
c===========================================================================c
c               horizontal and vertical eddy viscosity                      c
c===========================================================================c
	Do K = 1, KBM
!$OMP DO	
	   Do I = 1, IJM
	      If(VERTMIX .EQ. 'CONSTANT  ') Then
	         AAM(I,K) = SGS(I,K) + UMOL
	         KM(I,K)  = SGS(I,K) + UMOL
	      Else
	         KM(I,K)  = SGS(I,K) + UMOL
	         If(HORZMIX .NE. 'CLOSURE   ') AAM(I,K) = SGS(I,K) + UMOL
	      Endif
	   Enddo
!$OMP END DO
	Enddo
!$OMP END PARALLEL
c============================================================================c
c                    open boundary treatments                                c
c============================================================================c	   
C-----   elevation boundary condition
	   If(NUMEBC .NE. 0) Then
		 Do N = 1, NUMEBC
	        ID = IEBC(N)
	        Do K = 1, KBM
	           AAM(ID,K) = 0.0
	        Enddo
	     Enddo
	   Endif
C-----   astrotidal boundary condition
  	   If(NUMAST .NE. 0) Then
	     Do N = 1, NUMAST
	        ID = IABC(N)
	        Do K = 1, KBM
	           AAM(ID,K) = 0.0
	        Enddo
	     Enddo
	   Endif
c-----   discharge boundary condition
	   If(NUMQBC .NE. 0) Then
	      Do N = 1, NUMQBC
	         ID = IQBC(N)
	         Do K = 1, KBM
			    AAM(ID,K) = 0.0
	         Enddo
	      Enddo
	   Endif
C-----   velocity boundary condition
	   If(NUMVBC .NE. 0) Then
	      Do N = 1, NUMVBC
	         ID = IVBC(N)
	         Do K = 1, KBM
			    AAM(ID,K) = 0.0
	         Enddo
	      Enddo
	   Endif
c---------------------------------------------------------------------------c
	Return
	End