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
      Subroutine WSIGMA
      Include './Include/OCERM_INF'
      Dimension TEMP(KB)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,FXH,FYH,FXE,FYE,
!$OMP&                                 QZS,US,VS,TEMP)
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

C		       If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
C			      FXE = FXE + 
C     &		 	        WIX(I,J) * (ELF(CELL_SIDE(I,J,2)) - ELF(I))
C	            FYE = FYE + 
C     &			        WIY(I,J) * (ELF(CELL_SIDE(I,J,2)) - ELF(I))
C	         Endif
               If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
                  FXE = FXE + CELL_CUV(I,J,6) *
     &                 (ELF(CELL_SIDE(I,J,2)) + ELF(I)) / 2. *
     &                  CELL_CUV(I,J,7)
                  FYE = FYE + CELL_CUV(I,J,6) *
     &                 (ELF(CELL_SIDE(I,J,2)) + ELF(I)) / 2. *
     &                  CELL_CUV(I,J,8)                    
               Else
                  FXE = FXE + CELL_CUV(I,J,6) *
     &                  ELF(I) * CELL_CUV(I,J,7)
                  FYE = FYE + CELL_CUV(I,J,6) *
     &                  ELF(I) * CELL_CUV(I,J,8)     
               Endif            
		    Enddo
            FXE = FXE / AREA(I)
            FYE = FYE / AREA(I)
	      W(I,1) = 0.0
	      W(I,KB) = 0.0
	      TEMP(1) = 0.0
	      TEMP(KB) = 0.0
	      Do K = 1, KBM
	         PT(I,K) = QZ(I,K) / (HC(I) + ELF(I)) - 
     &                  ZZ(K) * (U(I,K) /(HC(I)+ELF(I)) * FXH / AREA(I)+
     &                          V(I,K) /(HC(I)+ELF(I)) * FYH / AREA(I))-
     &            (1.+ZZ(K))*(U(I,K) / (HC(I)+ELF(I)) * FXE +
     &                        V(I,K) / (HC(I)+ELF(I)) * FYE +
     &                      (ELF(I) - EL(I)) / DTI)
     
C	         W(I,K) =  .5*(QZ(I,K-1)+QZ(I,K)) / (HC(I) + ELF(I)) - 
C     &                   Z(K) * (.5*(U(I,K-1)+U(I,K))/(HC(I)+ELF(I)) * 
C     &                            FXH / AREA(I) + 
C     &                           .5*(V(I,K-1)+V(I,K))/(HC(I)+ELF(I)) *
C     &                            FYH / AREA(I)) -
C     &            (1.+Z(K))*(.5*(U(I,K-1)+U(I,K))/(HC(I)+ELF(I)) * FXE +
C     &                       .5*(V(I,K-1)+V(I,K))/(HC(I)+ELF(I)) * FYE +
C     &                    (ELF(I) - EL(I)) / DTI)

	      Enddo
CCC	      Do K = 2, KBM
CCC	         W(I,K) = (TEMP(K-1)*(Z(K)-ZZ(K)) + TEMP(K)*(ZZ(K-1)-Z(K)))/
CCC     &                  (ZZ(K-1)-ZZ(K))     
C	         W(I,K) = (TEMP(K-1) + TEMP(K)) / 2.
C	         W(I,K) = TEMP(K-1) * W_INTER(K,1) + 
C     &                  TEMP(K) * W_INTER(K,2) +
C     &                  TEMP(K+1) * W_INTER(K,3)        
            
CCC	      Enddo
CCC	      W(I,2) = 
CCC     &               TEMP(1)*W_INTER(2,2)+TEMP(2)*W_INTER(2,3)	      
CC	      W(I,KBM) = TEMP(KBM-1) * W_INTER(KBM,1) + 
CC     &                  TEMP(KBM) * W_INTER(KBM,2) +
CC     &                  TEMP(KB) * W_INTER(KBM,3) 
c            w(i,2) = TEMP(1)
C	      W(I,1) = QZ(I,1) / (HC(I) + ELF(I)) - 
C     &               ZZ(1) * (U(I,1) /(HC(I) + ELF(I)) * FXH / AREA(I) +
C     &                        V(I,1) /(HC(I) + ELF(I)) * FYH / AREA(I))-
C     &              (1.+ZZ(1))*(U(I,1) / (HC(I) + ELF(I)) * FXE +
C     &                          V(I,1) / (HC(I) + ELF(I)) * FYE +
C     &              (ELF(I) - EL(I)) / DTI)
CCC            W(I,1) = 2. * W(I,2) - W(I,3)
	   Endif
	Enddo
!$OMP END DO 
!$OMP END PARALLEL
C---------------------------------------------------------------------------c
       Return
       End
