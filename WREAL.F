C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (12/07/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c           computes the velocity in the vertical direction                 #
c                                                                           #
c############################################################################
      Subroutine WREAL
      Include './Include/OCERM_INF'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,FXH,FYH,FXE,FYE)
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      FXH = 0.0
	      FYH = 0.0
	      FXE = 0.0
	      FYE = 0.0
		    Do J = 1, CELL_POLYGEN(I)
		       FXH = FXH + 
     &	  		     HS(CELL_SIDE(I,J,1)) *
     &               CELL_CUV(I,J,7) * CELL_CUV(I,J,6)
	         FYH = FYH + 
     &               HS(CELL_SIDE(I,J,1)) *
     &               CELL_CUV(I,J,8) * CELL_CUV(I,J,6)
		       If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
			      FXE = FXE + 
     &		 	      WIX(I,J) * (ELF(CELL_SIDE(I,J,2)) - ELF(I))
	            FYE = FYE + 
     &			      WIY(I,J) * (ELF(CELL_SIDE(I,J,2)) - ELF(I))
	         Endif
		    Enddo
	      Do K = 1, KBM
	         WR(I,K) = (W(I,K) + W(I,K+1)) / 2. +
     &                   ZZ(K) * (U(I,K) / (HC(I) + ELF(I)) * 
     &                            FXH / AREA(I) + 
     &                            V(I,K) / (HC(I) + ELF(I)) *
     &                            FYH / AREA(I)) +
     &                   (1. + ZZ(K)) * (U(I,K) / (HC(I)+ELF(I)) * FXE +
     &                                   V(I,K) / (HC(I)+ELF(I)) * FYE +
     &                   (ELF(I) - EL(I)) / DTI)
C     &                  + DHT / DTI
c                QZ(I,K) = WR(I,K) * (HC(I) + ELF(I))

	      Enddo
	   Endif
	Enddo
!$OMP END DO
!$OMP END PARALLEL
C---------------------------------------------------------------------------c
       Return
       End
