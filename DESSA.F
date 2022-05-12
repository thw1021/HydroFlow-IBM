C###########################################################################C
c                 SUBROUNTINE PROGRAM                                       #
C                 SOLID BOUNDARY CONDITION BASED ON DES                     #
c                               VERSION 1.0 (16/02/2011)                    #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C############################################################################
	Subroutine DESSA
      Include './Include/OCERM_INF'
      Include './Include/VORGEN_INF'

	Parameter(CB1=0.2, CB2=0.622, VK=0.41, SIGMA=2./3.,CDES2 = 6.)
      Parameter(CW2=0.3, CW3=2.0, CV1=7.1, CDES=0.3, CT3=1.2, CT4=0.5)
	Dimension SGS(IJM,KB)
c===========================================================================c
c              transport equation of working variable                       c
c===========================================================================c
	  CW1 = CB1 / VK ** 2. + (1. + CB2) / SIGMA 
        Call ADVVIS(CB1,CB2,CW1,CW2,CW3,SIGMA,VK,CV1,CDES,CT3,CT4,CDES2)
	  Call PROFVIS(SIGMA)
c===========================================================================c
c               natural eddy viscosity                                      c
c===========================================================================c
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(COELES,ZSTAR,XX,FV1,D1,D2,
!$OMP&            DTIDLE,SLES,ZDES,D2MAX,I,J,K,SLOWERPART1,
!$OMP&            SLOWERPART2,SLOWER)
!$OMP DO
	  Do I = 1, IJM
	     If(CCM(I) .EQ. 1.0) Then
	     	D2MAX = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CELL_CUV(I,J,6) .GT. D2MAX) D2MAX = CELL_CUV(I,J,6) 
	        Enddo
                Do K = 1, KBM
	           COELES = 0.0
		         ZSTAR = DC(I) * (1. + ZZ(K)) / 1.E-6 * 
     &		           Sqrt(Sqrt(TBX(I)**2.+TBY(I)**2.)/RMEAN(I,KBM)) 
c	           VIS(I,K) = VIS(I,K) * Sqrt(1. - Exp(-(ZSTAR/25.)**3.))
  		   XX = VIS(I,K) / UMOL
CC	           XX = (VIS(I,K) + UMOL) / UMOL + 0.5 * 30. * Z01(I) / 
CC     &			    (DC(I) * (1. + ZZ(K)) + Z01(I))
		 FV1 = XX ** 3. / (XX ** 3. + CV1 ** 3.)
                 FV2 = 1.0 - XX / (1. + XX * FV1)
	           FT2 = CT3 * Exp(-CT4 * XX ** 2.) 
	           SLOWERPART1 = 100.
		         SLOWERPART2 = (1. - CB1 / CW1 / VK ** 2. / 0.424 * 
     &                         (FT2+(1. - FT2)* FV2))/Dmax1(FV1,1.E-10)/
     &                          Dmax1(1.-FT2, 1.E-10)
	           SLOWER = Sqrt(Dmin1(SLOWERPART1,SLOWERPART2))
C                   SLOWER = 1.0
C-----           ZDEX
                   If(DES .EQ. 'SAZDES ') Then
                     D1 = DC(I) * (1. + ZZ(K)) + Z01(I)
                     D1 = Dmin1(D1,D2D(I))
      	             D2 = CDES * Dmax1(D2MAX,DC(I)*DZ(K))
CCC	              D2 = CDES * Sqrt(4. * AREA(I) / ACOS(-1.0))
	              !LES length-scale DES97
                     If(I_DES_MODE(I,K) .EQ. 1) Then
C                        D2 = CDES*(AREA(I)*DC(I)*DZ(K))**(1./3.) 
C                        DTIDLE = Dmin1(D1,D2*SLOWER)
                        DTIDLE = D1
                        If(DTIDLE .GE. DWINTERFACE) FV1 = 1.0
                     Endif
                     If(I_DES_MODE(I,K) .EQ. 3) Then
C                        D2 = (AREA(I)*DC(I)*DZ(K))**(1./3.)
C                        D2 = Dmax1(3**2.*D2/D1*CDES*D2,D2)
                        DTIDLE = D1
C                        DTIDLE = Dmin1(D1,CDES2*D2)
                        If(DTIDLE .GE. DWINTERFACE) FV1 = 1.0
                     Endif
                   Endif   
		   SGS(I,K) = Dmax1(VIS(I,K) * FV1, 0.0 * UMOL)
	        Enddo
           Endif
        Enddo
C       stop
!$OMP END DO	  
c===========================================================================c
c               horizontal and vertical eddy viscosity                      c
c===========================================================================c
	Do K = 1, KBM
!$OMP DO	
	   Do I = 1, IJM
	      If(VERTMIX .EQ. 'CONSTANT  ') Then
	         AAM(I,K) = SGS(I,K) 
	         KM(I,K)  = SGS(I,K) 
              Else
	         KM(I,K)  = SGS(I,K) 
	         If(HORZMIX .NE. 'CLOSURE   ') AAM(I,K) = SGS(I,K) 
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
CC	           KM(ID,K) = 0.0
	        Enddo
	     Enddo
	   Endif
C-----   astrotidal boundary condition
  	   If(NUMAST .NE. 0) Then
	     Do N = 1, NUMAST
	        ID = IABC(N)
	        Do K = 1, KBM
	           AAM(ID,K) = 0.0
CC	           KM(ID,K) = 0.0
	        Enddo
	     Enddo
	   Endif
c-----   discharge boundary condition
	   If(NUMQBC .NE. 0) Then
	      Do N = 1, NUMQBC
	         ID = IQBC(N)
	         Do K = 1, KBM
			    AAM(ID,K) = 0.0
CC			    KM(ID,K) = 0.0
	         Enddo
	      Enddo
	   Endif
C-----   velocity boundary condition
	   If(NUMVBC .NE. 0) Then
	      Do N = 1, NUMVBC
	         ID = IVBC(N)
	         Do K = 1, KBM
			    AAM(ID,K) = 0.0
CC			    KM(ID,K) = 0.0
	         Enddo
	      Enddo
	   Endif
c---------------------------------------------------------------------------c
	Return
	End
