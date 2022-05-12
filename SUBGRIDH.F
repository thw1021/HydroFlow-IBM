C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (12/07/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                Horizontal eddy coefficient based on the Subgrid model     #
c                                                                           #
c############################################################################
      Subroutine SUBGRIDH
      Include './Include/OCERM_INF'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,N,ID,DIVUX,DIVUY,DIVVX,
!$OMP&                                 DIVVY) 
	Do K = 1, KBM
!$OMP DO
	   Do I = 1, IJM
	      DIVUX = 0.0                          ! du/dx
	      DIVUY = 0.0                          ! du/dy
	      DIVVX = 0.0                          ! dv/dx
	      DIVVY = 0.0                          ! dv/dy
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
CC	               DIVUX = DIVUX + 
CC     &			  	  WIX(I,J) * (UR(CELL_SIDE(I,J,2),K) - UR(I,K))
CC	               DIVUY = DIVUY +
CC     &                  WIY(I,J) * (UR(CELL_SIDE(I,J,2),K) - UR(I,K))
CC	               DIVVX = DIVVX + 
CC   &			  	  WIX(I,J) * (VR(CELL_SIDE(I,J,2),K) - VR(I,K))
CC                     DIVVY = DIVVY + 
CC     &			  	  WIY(I,J) * (VR(CELL_SIDE(I,J,2),K) - VR(I,K))
                  Endif
                  If(CFM(CELL_SIDE(I,J,1)) .EQ. 0.0 .OR.     
     &			     CFM(CELL_SIDE(I,J,1)) .EQ. -1.0) Then
                     If(ISLIP .EQ. 1) Then
	                  DIVUY = DIVUY + CELL_CUV(I,J,6) *
     &                          UR(I,K) * CELL_CUV(I,J,8)
	                  DIVVX = DIVVX + CELL_CUV(I,J,6) *
     &                          VR(I,K) * CELL_CUV(I,J,7)
CC	                  DIVUX = DIVUX - 2. * WIX(I,J) * UR(I,K)
CC	                  DIVUY = DIVUY - 2. * WIY(I,J) * UR(I,K)
CC	                  DIVVX = DIVVX - 2. * WIX(I,J) * VR(I,K)
CC	                  DIVVY = DIVVY - 2. * WIY(I,J) * VR(I,K)
	               Endif
                  Endif
                  If(CFM(CELL_SIDE(I,J,1)) .EQ. -3.0) Then
                     DIVUY = DIVUY + CELL_CUV(I,J,6) *
     &                       UR(I,K) * CELL_CUV(I,J,8)
	               DIVVX = DIVVX + CELL_CUV(I,J,6) *
     &                       VR(I,K) * CELL_CUV(I,J,7)
                  Endif 
	         Enddo	
               DIVUX = DIVUX / AREA(I)
	         DIVUY = DIVUY / AREA(I)
	         DIVVX = DIVVX / AREA(I)
	         DIVVY = DIVVY / AREA(I)
	         AAM(I,K) = HORCON * AREA(I) * Sqrt(
     &			 DIVUX **2. + DIVVY **2. + .5 * (DIVVX + DIVUY) **2.)
		  Endif
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
C----------------------------------------------------------------------------C
	Return
	End
