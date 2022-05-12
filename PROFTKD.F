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
      Subroutine PROFTKD(F1)
      Include './Include/OCERM_INF'
      Common/VISW/VISSOURCE(IJM,KB),DISSMARK(IJM,KB)
      Parameter (SIGMAO1=0.5,SIGMAO2=0.856,BETA1=0.075,BETA2=0.0828)
     	Dimension F1(IJM,KB),VISCOE(IJM,KB)
     	Dimension AAAA(KBM,KBM), BBBB(KBM)
	Dimension AA(3*KBM-2), BB(KBM)
      
!!!##################################################################!!! 
!!!   Feature    :Type of BC at free surface.                        !!!
!!!               (1: zeore-gradient, 2: value)                      !!!
!!!   UpdatedBy  :WangJian                                           !!!
!!!   UpdatedDate:2020-3-9 22:47:29                                  !!!
!!!##################################################################!!!      
      Itype_BCtop = 1
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
         RKM_TOP(K1)= 0.0
	Enddo
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K)
	Do K = 1, KB
!$OMP DO	
	   Do I = 1, IJM
		    VISCOE(I,K) = F1(I,K) * SIGMAO1 + (1. - F1(I,K)) * SIGMAO2
	   Enddo
!$OMP END DO	   
	Enddo
!$OMP END PARALLEL	
c===========================================================================c
c                  arrays for the turbulent kinetic governing equation      c
c===========================================================================c
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,K,K1,K2,AAAA,BBBB,AA,BB,L)
	Do I = 1, IJM
	   IF(CCM(I) .EQ. 1.0) Then
	 	   Do K = 2, KBM 
	        AAAA(K,K-1) = -DTI * (UMOL + (KM(I,K) + KM(I,K-1)) / 2. *
     &    	               (VISCOE(I,K) + VISCOE(I,K-1)) / 2.) / 
     &			            DC(I) ** 2. / DZZ(K-1) * DISSMARK(I,K)
	        AAAA(K-1,K) = AAAA(K,K-1) 
	     Enddo
	     Do K = 2, KBM - 1
	        AAAA(K,K) = DZ(K) - AAAA(K,K-1) - AAAA(K,K+1)
     &	                  + DTI * DZ(K) * VISSOURCE(I,K)
C     &                    (F1(I,K)*BETA1+(1.-F1(I,K))*BETA2)*TDISS(I,K)	
	     Enddo
	     Do K = 1, KBM
	        BBBB(K) = VISF(I,K) / AREA(I)
		 Enddo
c----      Boundary condition   --------------------------------------------c      
          IF (Itype_BCtop .EQ. 1 )THEN           
	     If(KBM .GT. 1) Then
c              AAAA(1,2) = 0.0
              AAAA(1,1) = DZ(1) - AAAA(1,2) + DTI * DZ(1) *
     &                    VISSOURCE(I,1)
c              IF(IBC_CONSISTENT .EQ. 1) THEN
c               AAAA(1,1) = DZ(1) - AAAA(1,2) + DTI*( UMOL + KM_TOP(I,1)
c     &                  * VISCOE(I,1) )/ DC(I)**2. * DIVTDISSZ_TOP(I,1)
c     &                  + DTI * DZ(1)* VISSOURCE(I,1)
c              ENDIF
	     Else
	        AAAA(1,1) = DZ(1) +
     &		              DTI * (KM(I,1)+UMOL) * VISCOE(I,1) /     
     &			          (.5 * DZ(1)) / DC(I) ** 2. + VISSOURCE(I,1)
	     Endif
          Endif           
c--------------------------------value------------------------------
          IF (Itype_BCtop .EQ. 2 )THEN 
	     If(KBM .GT. 1) Then
c---------------value at top layer center--------------------------------                  
c                AAAA(1,1) = 1.0
c                AAAA(1,2) = 0.0
c                BBBB(1) = TENSOR(I,1) / 0.3 
C              BBBB(1) = sqrt(TKE(I,1)+1.0E-20)*(0.09**0.75)/
C     &                    ( 0.07*0.41* DC(I) )                
c---------------value at top layer---------------------------------------              
             RKM_TOP(I) = VIS(I,1)-0.5*DZ(1)/DZZ(1)*(VIS(I,2)-VIS(I,1))
             AAAA(1,1) = DZ(1) - AAAA(1,2) 
     &    + DTI*(UMOL + RKM_TOP(I)*VISCOE(I,1))/DC(I)** 2./(0.5*DZ(1))
     &    + DTI * DZ(1)* VISSOURCE(I,1)
             
             BBBB(1) = BBBB(1)  
     &    +DTI*(UMOL + RKM_TOP(I)*VISCOE(I,1))/DC(I)**2./(0.5*DZ(1))
     &    *sqrt(TKE(I,1))*(0.09**0.75)/(0.07*0.41* DC(I))/0.09
      
	     Else
	        AAAA(1,1) = DZ(1)
              BBBB(1) = TENSOR(I,1) / 0.3   
	     Endif
          ENDIF 
c-----     at the bottom  ----------------------------------------c
		   If(WFBC .EQ. '      FUN1') Then      ! wall function for RANS
	        AAAA(KBM,KBM) = 1.0
	        BBBB(KBM) = Sqrt(Sqrt(TBX(I)** 2.+TBY(I)** 2.)/RMEAN(I,KBM))/ 
     &	                  0.3 / 0.41 / (DC(I) * DZZ(KBM))
	     Endif  
	     If(WFBC .EQ. '      FUN2') Then      ! wall function for LES

	     Endif
	     If(WFBC .EQ. '      FUN3') Then      ! no-slip B.C.
	       If(KBM .GT. 1) Then
	        AAAA(KBM,KBM) = DZ(KBM) - AAAA(KBM,KBM-1) *DISSMARK(I,KBM) + 
     &                        DTI * DISSMARK(I,KBM) *
     &                       (UMOL+KM(I,KBM)*VISCOE(I,KBM)) / 
     &                        DC(I) ** 2. / DZZ(KBM)
     &	                  + DTI * DZ(KBM) * VISSOURCE(I,KBM)
C     &                    (F1(I,KBM)*BETA1+(1.-F1(I,KBM))*BETA2)*
C     &                    TDISS(I,KBM)
     	
              BBBB(KBM) = BBBB(KBM) + DISSMARK(I,KBM)*
     &                   DTI * (UMOL+KM(I,KBM)*VISCOE(I,KBM))/
     &	                 DC(I) ** 2. / DZZ(KBM) * 60 * UMOL / BETA1 /
     &                   (DC(I) * DZZ(KBM)) ** 2.
	       Else
	        AAAA(KBM,KBM) = 1.0
	        BBBB(KBM) = TENSOR(I,KBM) / 0.3
C	        AAAA(KBM,KBM) = DZ(KBM) + DTI * 
C     &                        (UMOL+KM(I,KBM)*VISCOE(I,KBM)) / 
C     &                        DC(I) ** 2. / DZZ(KBM)
C     &	                      + DTI * DZ(KBM) *
C     &                       (F1(I,KBM)*BETA1+(1.-F1(I,KBM))*BETA2)*
C     &                        TDISS(I,KBM)

C              BBBB(KBM) = BBBB(KBM) + 
C     &                        DTI * (UMOL+KM(I,KBM)*VISCOE(I,KBM))/
C     &	                 DC(I) ** 2. / DZZ(KBM) * 60 * UMOL / BETA1 /
C     &                   (DC(I) * DZZ(KBM)) ** 2.
	       Endif
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
	        TDISS(I,K) = Max(BB(K),1.0E-10)
C	        IF(TDISS(I,K) .GT. 1E7) THEN
C	        PRINT*, I,K
C	        PAUSE
C	        ENDIF
	     Enddo
	   Endif
	Enddo
!$OMP END PARALLEL DO
C---------------------------------------------------------------------------C
      Return
      End

