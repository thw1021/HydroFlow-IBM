C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (18/11/2013)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c      computes the variable gradation                                      #
c                                                                           #
c############################################################################
      Subroutine GRAD_XY(VAR_EDGE,GRADX,GRADY)
      Include './Include/OCERM_INF'
	Dimension VAR_EDGE(IJE,KB),GRADX(IJM,KB),GRADY(IJM,KB)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K) 
	Do K = 1, KBM
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         GRADX(I,K) = 0.0
	         GRADY(I,K) = 0.0
 	         Do J = 1, CELL_POLYGEN(I)
	            GRADX(I,K) = GRADX(I,K) + 
     &	                       VAR_EDGE(CELL_SIDE(I,J,1),K) *
     &                         CELL_CUV(I,J,7) * CELL_CUV(I,J,6)
	            GRADY(I,K) = GRADY(I,K) + 
     &	                       VAR_EDGE(CELL_SIDE(I,J,1),K) *
     &                         CELL_CUV(I,J,8) * CELL_CUV(I,J,6)
	         Enddo
	         GRADX(I,K) = GRADX(I,K) / AREA(I)
	         GRADY(I,K) = GRADY(I,K) / AREA(I)
	      Endif
	   Enddo
!$OMP END DO 
	Enddo
!$OMP END PARALLEL	
	Return
	End
	
	Subroutine GRAD_Z(VAR_UD,GRADZ)
      Include './Include/OCERM_INF'
	Dimension VAR_UD(IJM,KB),GRADZ(IJM,KB)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K) 
	Do K = 1, KBM
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         GRADZ(I,K) = (VAR_UD(I,K) - VAR_UD(I,K+1)) / DZ(K) / DC(I)
	      Endif
	   Enddo
!$OMP END DO 
	Enddo
!$OMP END PARALLEL	
	Return
	End
	
      Subroutine GRAD2D(VAR_VERTEX,GRADX,GRADY)
      Include './Include/OCERM_INF'
	Dimension VAR_VERTEX(IJP),GRADX(IJM),GRADY(IJM)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K) 
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         GRADX(I) = 0.0
	         GRADY(I) = 0.0
 	         Do J = 1, CELL_POLYGEN(I)
	            GRADX(I) = GRADX(I) + 
     &	                       (VAR_VERTEX(CELL_SIDE(I,J,4)) +
     &	                        VAR_VERTEX(CELL_SIDE(I,J,3))) / 2. *
     &                         CELL_CUV(I,J,7) * CELL_CUV(I,J,6)
	            GRADY(I) = GRADY(I) + 
     &	                       (VAR_VERTEX(CELL_SIDE(I,J,4)) +
     &	                        VAR_VERTEX(CELL_SIDE(I,J,3))) / 2. *
     &                         CELL_CUV(I,J,8) * CELL_CUV(I,J,6)
	         Enddo
	         GRADX(I) = GRADX(I) / AREA(I)
	         GRADY(I) = GRADY(I) / AREA(I)
	      Endif
	   Enddo
!$OMP END DO 
!$OMP END PARALLEL	
	Return
	End
	
      Subroutine GRADXYADD(GRADX,GRADY,GRADZ,ADDX,ADDY)
      Include './Include/OCERM_INF'
	Dimension GRADX(IJM,KB),GRADY(IJM,KB),GRADZ(IJM,KB)
	Dimension ADDX(IJM,KB), ADDY(IJM,KB)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,ELFX,ELFY,FXH,FYH) 
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         ELFX = 0.0
	         ELFY = 0.0
	         FXH = 0.0
	         FYH = 0.0
	         Do J = 1, CELL_POLYGEN(I)
                  If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
                     ELFX = ELFX + CELL_CUV(I,J,6) *
     &                   (ELF(CELL_SIDE(I,J,2)) + ELF(I)) / 2. *
     &                    CELL_CUV(I,J,7)
                     ELFY = ELFY + CELL_CUV(I,J,6) *
     &                   (ELF(CELL_SIDE(I,J,2)) + ELF(I)) / 2. *
     &                    CELL_CUV(I,J,8)                    
                  Else
                    ELFX = ELFX + CELL_CUV(I,J,6) *
     &                   ELF(I) * CELL_CUV(I,J,7)
                    ELFY = ELFY + CELL_CUV(I,J,6) *
     &                   ELF(I) * CELL_CUV(I,J,8)     
                  Endif            
			      FXH = FXH + HS(CELL_SIDE(I,J,1)) *
     &                        CELL_CUV(I,J,7) * CELL_CUV(I,J,6)
	            FYH = FYH + HS(CELL_SIDE(I,J,1)) *
     &                        CELL_CUV(I,J,8) * CELL_CUV(I,J,6)
	         Enddo
	         ELFX = ELFX / AREA(I)
	         ELFY = ELFY / AREA(I)
	         FXH = FXH / AREA(I)
	         FYH = FYH / AREA(I)
	         Do K = 1, KBM
	            ADDX(I,K) = GRADZ(I,K) / DC(I) * (
     &                       (1. + ZZ(K)) * ELFX + ZZ(K) * FXH)
                  ADDY(I,K) = GRADZ(I,K) / DC(I) * (
     &                       (1. + ZZ(K)) * ELFY + ZZ(K) * FYH)       
	         Enddo
	      Endif
	   Enddo
!$OMP END DO 
!$OMP END PARALLEL	
	Return
	End
	

	
	
	
