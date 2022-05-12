C##########################################################################C
C                        SGS MODEL                                         C
C                   Version 1.0 (01/02/2011)                               C
C                 AUTHORIZED BY ZHANG JINGXIN                              C
C                               SHANGHAI JIAO TONG UNIVERSITY              C
C                               SHANGHAI, CHINA                            C
C--------------------------------------------------------------------------C         
C     SUBGRID STRESS MODEL FOR THE LARGE EDDY SIMULATION                   C
C                                                                          C
c##########################################################################C
	Subroutine SGSMODEL
      Include './Include/OCERM_INF'
	Dimension SGS(IJM,KB)
C---------------------------------------------------------------------------C
C         Different SGS models                                              c
c---------------------------------------------------------------------------c
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,DIVUX,DIVUY,DIVUZ,
!$OMP&  DIVVX,DIVVY,DIVVZ,DIVWX,DIVWY,DIVWZ)
	Goto (1) INDEX_SGS
C-----   Smagorinsky model  ------------------------------------------------c
1	Continue
	Do K = 1, KBM
C-----   wall unit  --------------------------------------------------------c
!$OMP DO
	   Do I = 1, IJM
	      SGS(I,K) = 0.0
	      FILTER_GRID = (AREA(I) * DZ(K) * DC(I)) ** 0.3333
	      DIVUX = 0.0                          ! du/dx
	      DIVUY = 0.0                          ! du/dy
	      DIVUZ = 0.0                          ! du/dz
	      DIVVX = 0.0                          ! dv/dx
	      DIVVY = 0.0                          ! dv/dy
	      DIVVZ = 0.0                          ! dv/dz
	      DIVWX = 0.0                          ! dw/dz
	      DIVWY = 0.0                          ! dw/dy
	      DIVWZ = 0.0                          ! dw/dz
	      If(CCM(I) .EQ. 1.0) Then
	         Do J = 1, CELL_POLYGEN(I)
	            If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.) Then
		             DIVUX = DIVUX + CELL_CUV(I,J,6) * 0.5 *
     &				   (UR(I,K) + UR(CELL_SIDE(I,J,2),K)) * CELL_CUV(I,J,7)
	               DIVUY = DIVUY + CELL_CUV(I,J,6) * 0.5 * 
     &				   (UR(I,K) + UR(CELL_SIDE(I,J,2),K)) * CELL_CUV(I,J,8)
	               DIVVX = DIVVX + CELL_CUV(I,J,6) * 0.5 *
     &				   (VR(I,K) + VR(CELL_SIDE(I,J,2),K)) * CELL_CUV(I,J,7)
		             DIVVY = DIVVY + CELL_CUV(I,J,6) * 0.5 *
     &				   (VR(I,K) + VR(CELL_SIDE(I,J,2),K)) * CELL_CUV(I,J,8)
	               DIVWX = DIVWX + CELL_CUV(I,J,6) * 
     &			  	   (WR(CELL_SIDE(I,J,2),K) + WR(I,K)) / 2. *
     &                   CELL_CUV(I,J,7)
	               DIVWY = DIVWY + CELL_CUV(I,J,6) * 
     &                   (WR(CELL_SIDE(I,J,2),K) + WR(I,K)) / 2. *
     &                   CELL_CUV(I,J,8)
C	               DIVUX = DIVUX + 
C     &			  	  WIX(I,J) * (UR(CELL_SIDE(I,J,2),K) - UR(I,K))
C	               DIVUY = DIVUY +
C     &                  WIY(I,J) * (UR(CELL_SIDE(I,J,2),K) - UR(I,K))
C	               DIVVX = DIVVX + 
C     &			  	  WIX(I,J) * (VR(CELL_SIDE(I,J,2),K) - VR(I,K))
C                     DIVVY = DIVVY + 
C     &			  	  WIY(I,J) * (VR(CELL_SIDE(I,J,2),K) - VR(I,K))
C	               DIVWX = DIVWX + 
C     &			  	  WIX(I,J) * (WR(CELL_SIDE(I,J,2),K) - WR(I,K))
C	               DIVWY = DIVWY +
C     &                  WIY(I,J) * (WR(CELL_SIDE(I,J,2),K) - WR(I,K))
                  Endif
                  If(CFM(CELL_SIDE(I,J,1)) .EQ. 0.0 .OR.     
     &			     CFM(CELL_SIDE(I,J,1)) .EQ. -1.0) Then
	               If(ISLIP .EQ. 1) Then
	                  DIVUY = DIVUY + CELL_CUV(I,J,6) *
     &					        UR(I,K) * CELL_CUV(I,J,8)
	                  DIVVX = DIVVX + CELL_CUV(I,J,6) *
     &                          VR(I,K) * CELL_CUV(I,J,7)
	                  DIVWX = DIVWX + CELL_CUV(I,J,6) *
     &                          WR(I,K) * CELL_CUV(I,J,7)
	                  DIVWY = DIVWY + CELL_CUV(I,J,6) *
     &                          WR(I,K) * CELL_CUV(I,J,8)
                     Endif
                  Endif
                  If(CFM(CELL_SIDE(I,J,1)) .EQ. -3.0) Then
	                  DIVUY = DIVUY + CELL_CUV(I,J,6) *
     &					        UR(I,K) * CELL_CUV(I,J,8)
	                  DIVVX = DIVVX + CELL_CUV(I,J,6) *
     &                          VR(I,K) * CELL_CUV(I,J,7)
	                  DIVWX = DIVWX + CELL_CUV(I,J,6) *
     &                          WR(I,K) * CELL_CUV(I,J,7)
	                  DIVWY = DIVWY + CELL_CUV(I,J,6) *
     &                          WR(I,K) * CELL_CUV(I,J,8)
                  Endif 
	         Enddo
	         DIVUX = DIVUX / AREA(I)
	         DIVUY = DIVUY / AREA(I)
	         DIVVX = DIVVX / AREA(I)
	         DIVVY = DIVVY / AREA(I)
	         DIVWX = DIVWX / AREA(I)
	         DIVWY = DIVWY / AREA(I)

c============================================================================c
c                    open boundary treatments                                c
c============================================================================c	   
C-----   elevation boundary condition
	         If(NUMEBC .NE. 0) Then
		        Do N = 1, NUMEBC
	               ID = IEBC(N)
	               If(ID .EQ. I) Then
	                  DIVUX = 0.0
	                  DIVUY = 0.0
	                  DIVVX = 0.0
	                  DIVVY = 0.0
	                Endif
	            Enddo
	         Endif
C-----   astrotidal boundary condition
  	         If(NUMAST .NE. 0) Then
	            Do N = 1, NUMAST
	               ID = IABC(N)
	               If(ID .EQ. I) Then
	                  DIVUX = 0.0
	                  DIVUY = 0.0
	                  DIVVX = 0.0
	                  DIVVY = 0.0
	                Endif
	            Enddo
	         Endif
c-----   discharge boundary condition
	         If(NUMQBC .NE. 0) Then
	            Do N = 1, NUMQBC
	               ID = IQBC(N)
	               If(ID .EQ. I) Then
	                  DIVUX = 0.0
	                  DIVUY = 0.0
	                  DIVVX = 0.0
	                  DIVVY = 0.0
	               Endif
	            Enddo
	         Endif
C-----   velocity boundary condition
	         If(NUMVBC .NE. 0) Then
	            Do N = 1, NUMVBC
	               ID = IVBC(N)
	               If(ID .EQ. I) Then
	                  DIVUX = 0.0
	                  DIVUY = 0.0
	                  DIVVX = 0.0
	                  DIVVY = 0.0
	               Endif
	            Enddo
	         Endif
			   If(KBM .GT. 1 .AND. K .EQ. 1) Then
	            DIVUZ = (UR(I,1) - .5 * (UR(I,1) + UR(I,2))) /
     &				    (.5 * DZ(1) * DC(I))
	            DIVVZ = (VR(I,1) - .5 * (VR(I,1) + VR(I,2))) / 
     &				    (.5 * DZ(1) * DC(I))

C			    DIVUZ = (UR(I,1) - 
C     &				     (UR(I,1) * DZ(2) + UR(I,2) * DZ(1)) / 2. /
C     &                     DZZ(1)) / (0.5 * DZ(1) * DC(I))
C			    DIVVZ = (VR(I,1) - 
C     &				     (VR(I,1) * DZ(2) + VR(I,2) * DZ(1)) / 2. /
C     &                     DZZ(1)) / (0.5 * DZ(1) * DC(I))
			    DIVWZ = (WR(I,1) - .5 * (WR(I,1) + WR(I,2))) /
     &				    (.5 * DZ(1) * DC(I))
C     &				     (WR(I,1) * DZ(2) + WR(I,2) * DZ(1)) / 2. /
C     &                     DZZ(1)) / (0.5 * DZ(1) * DC(I))
			 Endif
			 If(K .EQ. KBM) Then
		       If(KBM .GT. 1) Then
			       DIVUZ = (.5 * (UR(I,KBM-1) + UR(I,KBM)) - UR(I,KBM))/
     &		             (0.5 * DZ(KBM) * DC(I))
	             DIVVZ = (.5 * (VR(I,KBM-1) + VR(I,KBM)) - VR(I,KBM))/ 
     &		             (0.5 * DZ(KBM) * DC(I))
		           DIVWZ = (.5 * (WR(I,KBM-1) + WR(I,KBM)) - WR(I,KBM))/
     &		             (0.5 * DZ(KBM) * DC(I))
	           Else
                 DIVUZ = UR(I,K) / (0.5 * DZ(KBM) * DC(I))
			     DIVVZ = VR(I,K) / (0.5 * DZ(KBM) * DC(I))
			     DIVWZ = WR(I,K) / (0.5 * DZ(KBM) * DC(I))
	           Endif
C			    DIVUZ = UR(I,K) / 2. / DZ(KBM)
C			    DIVVZ = VR(I,K) / 2. / DZ(KBM)
C			    DIVWZ = WR(I,K) / 2. / DZ(KBM)
			 Endif
			 If(K .GT. 1 .AND. K .LT. KBM) Then
			    DIVUZ = ((UR(I,K-1) - UR(I,K)) / DZZ(K-1) / DC(I) + 
     &                     (UR(I,K) - UR(I,K+1)) / DZZ(K) / DC(I)) / 2.
			    DIVVZ = ((VR(I,K-1) - VR(I,K)) / DZZ(K-1) / DC(I) + 
     &                     (VR(I,K) - VR(I,K+1)) / DZZ(K) / DC(I)) / 2.
			    DIVWZ = ((WR(I,K-1) - WR(I,K)) / DZZ(K-1) / DC(I) + 
     &                     (WR(I,K) - WR(I,K+1)) / DZZ(K) / DC(I)) / 2.       
			 Endif

	         S11 = (DIVUX + DIVUX) / 2.
	         S12 = (DIVUY + DIVVX) / 2.
	         S13 = (DIVUZ + DIVWX) / 2.
			   S22 = (DIVVY + DIVVY) / 2.
			   S23 = (DIVVZ + DIVWY) / 2.
			   S33 = (DIVWZ + DIVWZ) / 2.  

	         SGS(I,K) = (UNICON * FILTER_GRID) ** 2. * 
     &			        Sqrt(2. * (S11**2. + S22**2. + S33**2. + 
     &                         2.*S12**2. + 2.*S13**2. + 2.*S23**2.))
C	PRINT*, SGS(I,K)
		       ZSTAR = DC(I) * (1. + ZZ(K)) / 1.E-6 * 
     &			     Sqrt(Sqrt(TBX(I)**2.+TBY(I)**2.)/RMEAN(I,KBM)) 
C	         SGS(I,K) = SGS(I,K) * Sqrt(1. - EXP(-(ZSTAR/25.)**3.))
		  Endif
	   Enddo
!$OMP END DO NOWAIT  
	Enddo
	Goto 100
C------------------------------------------------------------------------------C
c-------------------    end   -------------------------------------------------c
100	Continue
	Do K = 1, KBM
!$OMP DO	
	   Do I = 1, IJM
	      AAM(I,K) = SGS(I,K)
	      KM(I,K)  = SGS(I,K)
	   Enddo
!$OMP END DO NOWAIT
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
	Return
	End
