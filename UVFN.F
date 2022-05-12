C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (18/05/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c      1. compute the velocities on the sides by means of adding            #
c         anti-diffusive terms;                                             #
c      2. mark the upwind and downwind sides of every sides based           #
c         on the calculated normal velocities on the sides                  #
c############################################################################
      Subroutine UVFN
	Include './Include/OCERM_INF'
	Parameter(SCHEME = 1)
	Dimension PPPP(IJM,KB),DELPV(IJM,KB),AC(IJM,KB),AB(IJM,KB)
      Dimension PRO_NX(IJE), PRO_NY(IJE),PRO_XI(IJE),PRO_ETA(IJE)
C===========================================================================C
C               compute the velocities on the sides                         c
C===========================================================================C
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&         PRIVATE(I,J,K,IL,IR,UD,DISXI,EXIX,EXIY,ANTIDIFFU,N1,N2,
!$OMP&                 ANTIDIFFV,ELFX,ELFY,FXH,FYH,DIVUX,DIVUY,DIVUZ,
!$OMP&                 DIVVX,DIVVY,DIVVZ,DIVWX,DIVWY,DIVWZ,S11,S12,
!$OMP&                 S13,S22,S23,S33,SXX,SYY,SZZ,DISKXI,DISETA,
!$OMP&                 X1,Y1,X2,Y2,X3,Y3,X4,Y4,CROSS,UFN) 

C	Goto (1, 2) SCHEME
C---------------------------------------------------------------------------C
C                  gradient of variant P                                    c
c---------------------------------------------------------------------------c	 
C1	Continue
!$OMP DO
      Do I = 1, IJE
         If(CFM(I) .EQ. 1.0) Then
            IL = INDEX_EDGE(I,1,1)
	      IR = INDEX_EDGE(I,1,2) 
	      N1 = IEND_EDGE(I,1)
	      N2 = IEND_EDGE(I,2)
	      X1 = CXY(IL,1)
	      Y1 = CXY(IL,2)
	      X2 = CXY(IR,1)
	      Y2 = CXY(IR,2)
	      X3 = PXY(N1,1)
	      Y3 = PXY(N1,2)
	      X4 = PXY(N2,1)
	      Y4 = PXY(N2,2)
	      CROSS = (X2 - X1) * (Y4 - Y3) - (X4 - X3) * (Y2 - Y1)
	      If(CROSS .LT. 0.0) Then
	         X3 = PXY(N2,1)
	         Y3 = PXY(N2,2)
	         X4 = PXY(N1,1)
	         Y4 = PXY(N1,2)	         
	      Endif    
	      DISKXI = Sqrt((X1 - X2) ** 2. + (Y1 - Y2) ** 2.)
            DISETA = Sqrt((X3 - X4) ** 2. + (Y3 - Y4) ** 2.)
            PRO_XI(I) = (X2 - X1) / DISKXI
            PRO_ETA(I) = (Y2 - Y1) / DISKXI
            PRO_NX(I) = (Y4 - Y3) / DISETA
            PRO_NY(I) = (X3 - X4) / DISETA
	   Endif   
      Enddo
!$OMP END DO
!$OMP DO
	Do I = 1, IJP
	   ELFV(I) = 0.0
C	   Do K = 1, KB
C	      PNV(I,K) = 0.0
C         Enddo
      Enddo
!$OMP END DO         
      Do K = 1, 1
!$OMP DO
	   Do I = 1, IJP
	      Do J = 1, POINT_CELL(I)
	         If(CCM(INDEX_POINT(I,J)) .EQ. 1.0) Then
                  If(K .EQ. 1)
     &			  ELFV(I) = ELFV(I) + ELF(INDEX_POINT(I,J)) * RCOE(I,J)
C                  PNV(I,K) = PNV(I,K)+PN(INDEX_POINT(I,J),K) * RCOE(I,J)
	         Endif
	      Enddo
	   Enddo
!$OMP END DO
      Enddo	         
      Do K = 1, KB
!$OMP DO
         Do I = 1, IJM 
            AB(I,K)=0.0
            AC(I,K)=0.0
            DELP(I,K,1) = 0.0
            DELP(I,K,2) = 0.0
         Enddo
!$OMP END DO NOWAIT
      Enddo
      If(SCHEME .EQ.1) Then
	Do K = 1, 1
!$OMP DO	
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	        DELP(I,K,1) = 0.0
	        DELP(I,K,2) = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	             DELP(I,K,1) = DELP(I,K,1) + 
     &			   (ELF(I) + ELF(CELL_SIDE(I,J,2))) / 2. *
     &               CELL_CUV(I,J,6) * CELL_CUV(I,J,7) 
	             DELP(I,K,2) = DELP(I,K,2) + 
     &			   (ELF(I) + ELF(CELL_SIDE(I,J,2))) / 2. *
     &               CELL_CUV(I,J,6) * CELL_CUV(I,J,8) 
	           Else
	             DELP(I,K,1) = DELP(I,K,1) + 
     &			   ELF(I) * CELL_CUV(I,J,6) * CELL_CUV(I,J,7) 
	             DELP(I,K,2) = DELP(I,K,2) + 
     &			   ELF(I) * CELL_CUV(I,J,6) * CELL_CUV(I,J,8) 
	           Endif
     	        Enddo
	        DELP(I,K,1) = DELP(I,K,1) / AREA(I)
	        DELP(I,K,2) = DELP(I,K,2) / AREA(I)
              
	        GRADP(I,K,1) = RMEAN(I,K)*GRAV*DELP(I,K,1)   
	        GRADP(I,K,2) = RMEAN(I,K)*GRAV*DELP(I,K,2)               
              
	      Endif
	   Enddo
!$OMP END DO NOWAIT 
	Enddo
	Endif
C	Goto 100
C2	Continue
      If(SCHEME .EQ. 2) Then
	Do K = 1, 1
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         DELP(I,K,1) = 0.0
	         DELP(I,K,2) = 0.0
	         Do J = 1, CELL_POLYGEN(I)
	            If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.) Then
	               DELP(I,K,1) = DELP(I,K,1) + 
     &			  	  WIX(I,J) * (ELF(CELL_SIDE(I,J,2)) - ELF(I))
	               DELP(I,K,2) = DELP(I,K,2) +
     &                  WIY(I,J) * (ELF(CELL_SIDE(I,J,2)) - ELF(I)) 
                    Endif 
	         Enddo
	         GRADP(I,K,1) = RMEAN(I,K)*GRAV*DELP(I,K,1)   
	         GRADP(I,K,2) = RMEAN(I,K)*GRAV*DELP(I,K,2)          	        
	      Endif
	   Enddo
!$OMP END DO NOWAIT
	Enddo
	Endif
C	Goto 100
c----------------------------------------------------------------------------c
C100	Continue
	Do K = 1, KBM
!$OMP DO	
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         AC(I,K) = GRAV * DTI * DC(I)
	      Endif
	   Enddo
!$OMP END DO NOWAIT
	Enddo
!$OMP BARRIER
!$OMP DO
	Do I = 1, IJE
         If(CFM(I) .EQ. 1.0) Then
	      IL = INDEX_EDGE(I,1,1)
	      IR = INDEX_EDGE(I,1,2)
	      N1 = IEND_EDGE(I,1)
	      N2 = IEND_EDGE(I,2)
              DISETA = Sqrt((PXY(N2,1) - PXY(N1,1)) ** 2. + 
     &                    (PXY(N2,2) - PXY(N1,2)) ** 2.)
	      DISXI = Sqrt((CXY(IR,1) - CXY(IL,1)) ** 2. +
     &                   (CXY(IR,2) - CXY(IL,2)) ** 2.)
	      Do K = 1, KBM
	         UN(I,K) = 0.0
	         VN(I,K) = 0.0

	         ANTIDIFFU = 0.5 * (AC(IL,K) + AC(IR,K)) *  
     &                     ((DELP(IL,1,1) + DELP(IR,1,1)) / 2. -            
     &                (ELF(IR) - ELF(IL)) / DISXI * PRO_XI(I))
	         ANTIDIFFV = 0.5 * (AC(IL,K) + AC(IR,K)) *  
     &                     ((DELP(IL,1,2) + DELP(IR,1,2)) / 2. -
     &                (ELF(IR) - ELF(IL)) / DISXI * PRO_ETA(I))

C	         ANTIDIFFU = -0.5 * (AC(IL,K) + AC(IR,K)) *  
C     &                     (ELF(IR) - ELF(IL)) / DISXI * PRO_XI(I) +
C     &                     (AC(IL,K) * DELP(IL,1,1) + 
C     &                      AC(IR,K) * DELP(IR,1,1)) / 2.
C	         ANTIDIFFV = -0.5 * (AC(IL,K) + AC(IR,K)) *  
C     &                     (ELF(IR) - ELF(IL)) / DISXI * PRO_ETA(I) +
C     &                     (AC(IL,K) * DELP(IL,1,2) + 
C     &                      AC(IR,K) * DELP(IR,1,2)) / 2.
     
C	         ANTIDIFFU = -0.5 * (AC(IL,K) + AC(IR,K)) *  
C     &                     (ELF(IR) - ELF(IL)) / DISXI  -
C     &                     ((AC(IL,K) * DELP(IL,1,1) + 
C     &                       AC(IR,K) * DELP(IR,1,1)) / 2. * PRO_XI(I) +
C     &                      (AC(IL,K) * DELP(IL,1,2) +
C     &                       AC(IR,K) * DELP(IR,1,2)) / 2. * PRO_ETA(I))
c	         UN(I,K) = (U(IL,K) * Sqrt(AREA(IR)) + 
c     &			        U(IR,K) * Sqrt(AREA(IL))) /
c     &                   (Sqrt(AREA(IL)) + Sqrt(AREA(IR))) + ANTIDIFFU
c	         VN(I,K) = (V(IL,K) * Sqrt(AREA(IR)) + 
c     &			        V(IR,K) * Sqrt(AREA(IL))) /
c     &                   (Sqrt(AREA(IL)) + Sqrt(AREA(IR))) + ANTIDIFFV
C 	        UFN =  0.5 * (U(IL,K) + U(IR,K)) * PRO_NX(I) +
C     &               0.5 * (V(IL,K) + V(IR,K)) * PRO_NY(I) + ANTIDIFFU
C			  UN(I,K) = UFN * PRO_NX(I)
C			  VN(I,K) = UFN * PRO_NY(I)
			  UN(I,K) = 0.5 * (U(IL,K) + U(IR,K)) + ANTIDIFFU 
			  VN(I,K) = 0.5 * (V(IL,K) + V(IR,K)) + ANTIDIFFV 
	      Enddo
	   Endif
	Enddo
!$OMP END DO

C===========================================================================C
C            hydrostatic fluid flow                                         c
C===========================================================================C
	Do K = 1, KBM
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then                
	        GRADP(I,K,1) = GRADP(I,1,1)
	        GRADP(I,K,2) = GRADP(I,1,2)              
                GRADP(I,K,3) = -RMEAN(I,K)*GRAV 
            Endif
         ENDDO
!$OMP END DO NOWAIT         
      ENDDO                 

C===========================================================================C
C            non-hydrostatic fluid flow                                     c
C===========================================================================C
	If (HYDTYPE .EQ. 'NONSTATIC') Then
         
C	Goto (10, 20) SCHEME
C---------------------------------------------------------------------------C
C                  gradient of variant P                                    c
c---------------------------------------------------------------------------c	 
C10	Continue	
      If(SCHEME .EQ. 1) Then
	Do K = 1, KBM
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then                                      
	        DELP(I,K,1) = 0.0
	        DELP(I,K,2) = 0.0
                DELP(I,K,3) = 0.0
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	             DELP(I,K,1) = DELP(I,K,1) + 
     &			   (PN(I,K) + PN(CELL_SIDE(I,J,2),K)) / 2. *
     &               CELL_CUV(I,J,6) * CELL_CUV(I,J,7) 
	             DELP(I,K,2) = DELP(I,K,2) + 
     &		  	   (PN(I,K) + PN(CELL_SIDE(I,J,2),K)) / 2. *
     &               CELL_CUV(I,J,6) * CELL_CUV(I,J,8) 
	           Else
	             DELP(I,K,1) = DELP(I,K,1) + 
     &			   PN(I,K) * CELL_CUV(I,J,6) * CELL_CUV(I,J,7) 
	             DELP(I,K,2) = DELP(I,K,2) + 
     &		  	   PN(I,K) * CELL_CUV(I,J,6) * CELL_CUV(I,J,8) 
	           Endif
     	        Enddo
	        DELP(I,K,1) = DELP(I,K,1) / AREA(I)
	        DELP(I,K,2) = DELP(I,K,2) / AREA(I)
              
c-----z pressure gradient              
              If(K .EQ. 1)Then
                  DELP(I,K,3) = ( PN(I,K) - PN(I,K+1) ) / 
     &                          ( DZZ(K) * DC(I) ) 
              Elseif(K .EQ. KBM)Then
                  DELP(I,K,3) = ( PN(I,K-1) - PN(I,K) ) / 
     &                          ( DZZ(K-1) * DC(I) ) 
              Else                                    
                  DELP(I,K,3) = ( PN(I,K-1) - PN(I,K+1) ) / 
     &                          ( (Z(K-1) - Z(K+1)) * DC(I) )
              Endif  
c------total pressure gradient              
	      GRADP(I,K,1) = GRADP(I,K,1)+DELP(I,K,1)   
	      GRADP(I,K,2) = GRADP(I,K,2)+DELP(I,K,2)   
              GRADP(I,K,3) = GRADP(I,K,3)+DELP(I,K,3)                 
	      Endif
	   Enddo
!$OMP END DO NOWAIT
	Enddo
	Endif
C	Goto 200
C20	Continue
      If(SCHEME .EQ. 2) Then
	Do K = 1, KBM
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         DELP(I,K,1) = 0.0
	         DELP(I,K,2) = 0.0
                 DELP(I,K,3) = 0.0
	         Do J = 1, CELL_POLYGEN(I)
		          If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	               DELP(I,K,1) = DELP(I,K,1) + 
     &			  	     WIX(I,J) * (PN(CELL_SIDE(I,J,2),K) - PN(I,K))
	               DELP(I,K,2) = DELP(I,K,2) +
     &                   WIY(I,J) * (PN(CELL_SIDE(I,J,2),K) - PN(I,K))  
	            Endif
	         Enddo
	         DELP(I,K,1) = DELP(I,K,1) / AREA(I)
	         DELP(I,K,2) = DELP(I,K,2) / AREA(I)
c-----z pressure gradient              
                IF(K .EQ. 1)THEN
                  DELP(I,K,3) = ( PN(I,K) - PN(I,K+1) ) / 
     &                          ( DZZ(K) * DC(I) ) 
                ELSEIF(K .EQ. KBM)THEN
                  DELP(I,K,3) = ( PN(I,K-1) - PN(I,K) ) / 
     &                          ( DZZ(K-1) * DC(I) ) 
                ELSE                                    
                  DELP(I,K,3) = ( PN(I,K-1) - PN(I,K+1) ) / 
     &                          ( (Z(K-1) - Z(K+1)) * DC(I) )
                ENDIF 
c------total pressure gradient              
	        GRADP(I,K,1) = GRADP(I,K,1)+DELP(I,K,1)   
	        GRADP(I,K,2) = GRADP(I,K,2)+DELP(I,K,2)   
                GRADP(I,K,3) = GRADP(I,K,3)+DELP(I,K,3)  
	      Endif
	   Enddo
!$OMP END DO NOWAIT
	Enddo
	Endif
C	Goto 200
c----------------------------------------------------------------------------c
C200	Continue
	Do K = 1, KBM
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
C	         AC(I,K) = DTI * DZ(K) * DC(I) / RMEAN(I,K) / AB(I,K)
	         AC(I,K) = DTI * DC(I) / RMEAN(I,K)
	      Endif
	   Enddo
!$OMP END DO NOWAIT
	Enddo
!$OMP BARRIER
!$OMP DO
	Do I = 1, IJE
	   If(CFM(I) .EQ. 1.0) Then
	      IL = INDEX_EDGE(I,1,1)
	      IR = INDEX_EDGE(I,1,2)
	      N1 = IEND_EDGE(I,1)
	      N2 = IEND_EDGE(I,2)
            DISETA = Sqrt((PXY(N2,1) - PXY(N1,1)) ** 2. + 
     &                    (PXY(N2,2) - PXY(N1,2)) ** 2.)
	      DISXI = Sqrt((CXY(IR,1) - CXY(IL,1)) ** 2. +
     &                   (CXY(IR,2) - CXY(IL,2)) ** 2.)
            Do K = 1, KBM 
C-----      anti-diffusion
	         ANTIDIFFU = 0.5 * (AC(IL,K) + AC(IR,K)) *  
     &                     ((DELP(IL,K,1) + DELP(IR,K,1)) / 2. -
     &                     (PN(IR,K) - PN(IL,K)) / DISXI * PRO_XI(I))
     
	         ANTIDIFFV = 0.5 * (AC(IL,K) + AC(IR,K)) *  
     &                     ((DELP(IL,K,2) + DELP(IR,K,2)) / 2. -
     &                      (PN(IR,K) - PN(IL,K)) / DISXI * PRO_ETA(I))
     
			   UN(I,K) = UN(I,K) + ANTIDIFFU 
			   VN(I,K) = VN(I,K) + ANTIDIFFV 
	      Enddo
	   Endif
	Enddo
!$OMP END DO
c---------------------------------------------------------------------------c
c                adding diffusive terms to the vertical velocity            c
c---------------------------------------------------------------------------c
c	Do I = 1, IJM
c	   If(CCM(I) .EQ. 1.0) Then
c	      Do K = 2, KBM - 1
c	         DELPV(I,K) = 0.5 * ((PN(I,K-1) - PN(I,K)) / DZZ(K-1) + 
c     &                             (PN(I,K) - PN(I,K+1)) / DZZ(K))
c	      Enddo
c	      DELPV(I,1) = (PN(I,1) - PN(I,2)) / DZZ(1)
c	      DELPV(I,KBM) = 0.5 * (PN(I,KBM-1) - PN(I,KBM)) / DZZ(KBM-1)
c	   Endif
c	Enddo
c	Do I = 1, IJM
c	   If(CCM(I) .EQ. 1.0) Then
c	      COE = DTI / DC(I) 
c		  Do K = 2, KBM 
c	         W(I,K) = W(I,K) - 
c     &			      0.5 * (COE / RMEAN(I,K-1) + COE / RMEAN(I,K)) *
c     &                  (PN(I,K-1) - PN(I,K)) / DZZ(K-1) +
c     &                  0.5 * (COE / RMEAN(I,K-1) * DELPV(I,K-1) +
c     &                         COE / RMEAN(I,K) * DELPV(I,K)) 
c	      Enddo
c	   Endif
c	EnddO
	Endif
	Do K = 1, KBM
!$OMP	DO
	   Do I = 1, IJE
	      If(CFM(I) .EQ. 1.0) Then
	         UN(I,K) = UN(I,K) / DS(I)
	         VN(I,K) = VN(I,K) / DS(I)
	      Endif
	   Enddo
!$OMP END DO
	Enddo
C======================================================================================C
C-----   TKE  CALCULATION   VORTICITY CALCULATION   -----------------------------------c
C======================================================================================C
	If(ITK .EQ. 1) Then
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         ELFX = 0.0
	         ELFY = 0.0
	         FXH = 0.0
	         FYH = 0.0
	         Do J = 1, CELL_POLYGEN(I)
	            If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	               ELFX = ELFX + 
     &				      WIX(I,J) * (ELF(CELL_SIDE(I,J,2)) - ELF(I))
	               ELFY = ELFY + 
     &				      WIY(I,J) * (ELF(CELL_SIDE(I,J,2)) - ELF(I))
	            Endif
			      FXH = FXH + HS(CELL_SIDE(I,J,1)) *
     &                        CELL_CUV(I,J,7) * CELL_CUV(I,J,6)
	            FYH = FYH + HS(CELL_SIDE(I,J,1)) *
     &                        CELL_CUV(I,J,8) * CELL_CUV(I,J,6)
	         Enddo
	         FXH = FXH / AREA(I)
	         FYH = FYH / AREA(I)
	         Do K = 1, KBM
	            DIVUX = 0.0                          ! du/dx
	            DIVUY = 0.0                          ! du/dy
	            DIVUZ = 0.0                          ! du/dz
	            DIVVX = 0.0                          ! dv/dx
	            DIVVY = 0.0                          ! dv/dy
	            DIVVZ = 0.0                          ! dv/dz
	            DIVWX = 0.0                          ! dw/dx
	            DIVWY = 0.0                          ! dw/dy
	            DIVWZ = 0.0                          ! dw/dz
	            S11 = 0.0
	            S12 = 0.0
	            S13 = 0.0
	            S22 = 0.0
	            S23 = 0.0
	            S33 = 0.0
	            SXX = 0.0
	            SYY = 0.0
	            SZZ = 0.0
	            Do J = 1, CELL_POLYGEN(I)
	               If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.) Then
	                  DIVUX = DIVUX + CELL_CUV(I,J,6) * 0.5 *
     &				          (UR(I,K) + UR(CELL_SIDE(I,J,2),K)) * CELL_CUV(I,J,7)

	                  DIVUY = DIVUY + CELL_CUV(I,J,6) * 0.5 * 
     &				          (UR(I,K) + UR(CELL_SIDE(I,J,2),K)) * CELL_CUV(I,J,8)

                        DIVVX = DIVVX + CELL_CUV(I,J,6) * 0.5 *
     &				          (VR(I,K) + VR(CELL_SIDE(I,J,2),K)) * CELL_CUV(I,J,7)

                        DIVVY = DIVVY + CELL_CUV(I,J,6) * 0.5 *
     &				          (VR(I,K) + VR(CELL_SIDE(I,J,2),K)) * CELL_CUV(I,J,8)

	                  DIVWX = DIVWX + CELL_CUV(I,J,6) * 
     &			  	            (WR(CELL_SIDE(I,J,2),K) + WR(I,K)) / 2.*
     &                          CELL_CUV(I,J,7) 

	                  DIVWY = DIVWY + CELL_CUV(I,J,6) *
     &                          (WR(CELL_SIDE(I,J,2),K) + WR(I,K)) / 2.*
     &                          CELL_CUV(I,J,8)
                     Endif
                     If(CFM(CELL_SIDE(I,J,1)) .EQ. 0.0 .OR.     
     &			        CFM(CELL_SIDE(I,J,1)) .EQ. -1.0) Then
                        If(ISLIP .EQ. 1) Then
	                     DIVUY = DIVUY + CELL_CUV(I,J,6) *
     &                             UR(I,K) * CELL_CUV(I,J,8)

	                     DIVVX = DIVVX + CELL_CUV(I,J,6) *
     &                             VR(I,K) * CELL_CUV(I,J,7)

	                     DIVWX = DIVWX + CELL_CUV(I,J,6) *
     &  					       WR(I,K) * CELL_CUV(I,J,7)

	                     DIVWY = DIVWY + CELL_CUV(I,J,6) *
     &					           WR(I,K) * CELL_CUV(I,J,8)
	                  Endif
                     Endif
                     If(CFM(CELL_SIDE(I,J,1)) .EQ. -3.0) Then
                        DIVUY = DIVUY + CELL_CUV(I,J,6) *
     &                          UR(I,K) * CELL_CUV(I,J,8)

	                  DIVVX = DIVVX + CELL_CUV(I,J,6) *
     &                          VR(I,K) * CELL_CUV(I,J,7)

                        DIVWX = DIVWX + CELL_CUV(I,J,6) *
     &  				        WR(I,K) * CELL_CUV(I,J,7)

	                  DIVWY = DIVWY + CELL_CUV(I,J,6) *
     &					        WR(I,K) * CELL_CUV(I,J,8)
                     Endif
	            Enddo
	            DIVUX = DIVUX / AREA(I)
	            DIVUY = DIVUY / AREA(I)
	            DIVVX = DIVVX / AREA(I)
	            DIVVY = DIVVY / AREA(I)
	            DIVWX = DIVWX / AREA(I)
	            DIVWY = DIVWY / AREA(I)
		    If(KBM .GT. 1 .AND. K .EQ. 1) Then
	               DIVUZ = (UR(I,1) - UR(I,2)) / DZZ(1) / DC(I) 
	               DIVVZ = (VR(I,1) - VR(I,2)) / DZZ(1) / DC(I)
	               DIVWZ = (WR(I,1) - WR(I,2)) / DZZ(1) / DC(I)
			      Endif
			      If(K .EQ. KBM) Then
		           If(KBM .GT. 1) Then
C                    DIVUZ = (UR(I,KBM-1) + UR(I,KBM)) / 2.0 /
C     &                      (DZ(KBM)*DC(I))
C	              DIVVZ = (VR(I,KBM-1) + VR(I,KBM)) / 2.0 /
C     & 				      (DZ(KBM)*DC(I))
C	              DIVWZ = (WR(I,KBM-1) + WR(I,KBM)) / 2.0 /
C     & 				      (DZ(KBM)*DC(I))
 		          DIVUZ = ((UR(I,KBM-1)-UR(I,KBM)) / (DZZ(KBM)*DC(I))+
     & 				       UR(I,KBM) / (.5 * DZ(KBM) * DC(I))) / 2.
	              DIVVZ = ((VR(I,KBM-1)-VR(I,KBM)) / (DZZ(KBM)*DC(I))+
     & 				       VR(I,KBM) / (.5 * DZ(KBM) * DC(I))) / 2.
	              DIVWZ = ((WR(I,KBM-1)-WR(I,KBM)) / (DZZ(KBM)*DC(I))+
     & 				       WR(I,KBM) / (.5 * DZ(KBM) * DC(I))) / 2.

	              Else
                      DIVUZ = UR(I,K) / (0.5 * DZ(KBM) * DC(I))
			        DIVVZ = VR(I,K) / (0.5 * DZ(KBM) * DC(I))
			        DIVWZ = WR(I,K) / (0.5 * DZ(KBM) * DC(I))
	              Endif
			      Endif
			      If(K .GT. 1 .AND. K .LT. KBM) Then
			      DIVUZ = ((UR(I,K-1) - UR(I,K)) / DZZ(K-1) / DC(I) + 
     &                       (UR(I,K) - UR(I,K+1)) / DZZ(K) / DC(I)) /2.
			      DIVVZ = ((VR(I,K-1) - VR(I,K)) / DZZ(K-1) / DC(I) + 
     &                       (VR(I,K) - VR(I,K+1)) / DZZ(K) / DC(I)) /2.
			      DIVWZ = ((WR(I,K-1) - WR(I,K)) / DZZ(K-1) / DC(I) + 
     &                       (WR(I,K) - WR(I,K+1)) / DZZ(K) / DC(I)) /2.
			      Endif
c-----            !#Sigma coordinate correction#, WangJian, 2019-12-16 15:02:34 !    
c	            DIVUX = DIVUX - DIVUZ * (ZZ(K)*FXH + (1.+ZZ(K))*ELFX)/DC(I)
c	            DIVUY = DIVUY - DIVUZ * (ZZ(K)*FYH + (1.+ZZ(K))*ELFY)/DC(I)
c	            DIVVX = DIVVX - DIVVZ * (ZZ(K)*FXH + (1.+ZZ(K))*ELFX)/DC(I)
c	            DIVVY = DIVVY - DIVVZ * (ZZ(K)*FYH + (1.+ZZ(K))*ELFY)/DC(I)
c	            DIVWX = DIVWX - DIVWZ * (ZZ(K)*FXH + (1.+ZZ(K))*ELFX)/DC(I)
c	            DIVWY = DIVWY - DIVWZ * (ZZ(K)*FYH + (1.+ZZ(K))*ELFY)/DC(I)
c-----            Stress tension (DUi/DXj+DUj/DXi)/2.				 
	            S11 = (DIVUX + DIVUX) / 2.
	            S12 = (DIVUY + DIVVX) / 2.
	            S13 = (DIVUZ + DIVWX) / 2.
			      S22 = (DIVVY + DIVVY) / 2.
			      S23 = (DIVVZ + DIVWY) / 2.
			      S33 = (DIVWZ + DIVWZ) / 2. 
                    
			      STRESS(I,K,1) = DIVUX
			      STRESS(I,K,2) = DIVUY
			      STRESS(I,K,3) = DIVUZ
			      STRESS(I,K,4) = DIVVX
			      STRESS(I,K,5) = DIVVY
			      STRESS(I,K,6) = DIVVZ
			      STRESS(I,K,7) = DIVWX
			      STRESS(I,K,8) = DIVWY
			      STRESS(I,K,9) = DIVWZ
				 
	            TENSOR(I,K) = Sqrt(2. * (S11**2. + S22**2. + S33**2. + 
     &                            2.*S12**2. + 2.*S13**2. + 2.*S23**2.))
c-----            Vorticity (DUi/DXj-DUj/DXi)
                  SXX = (DIVWY - DIVVZ) / 2.
	            SYY = (DIVUZ - DIVWX) / 2.
	            SZZ = (DIVVX - DIVUY) / 2.

	            VORX(I,K) =  SXX * 2.
	            VORY(I,K) =  SYY * 2.
	            VORZ(I,K) =  SZZ * 2.
	            VORTICITY(I,K) = 
     &              Sqrt(VORX(I,K)**2. + VORY(I,K)**2. + VORZ(I,K)**2.)
	         Enddo
	      Endif
	   Enddo
!$OMP END DO
	Endif
!$OMP END PARALLEL
c-----   open boundary treatments                                
C-----   elevation boundary condition
	   If(NUMEBC .NE. 0) Then
		 Do N = 1, NUMEBC
	        ID = IEBC(N)
	        Do K = 1, KBM
C	           TENSOR(ID,K) = 0.0
C	           VORTICITY(ID,K) = 0.0
	        Enddo
	     Enddo
	   Endif
C-----   astrotidal boundary condition
  	   If(NUMAST .NE. 0) Then
	     Do N = 1, NUMAST
	        ID = IABC(N)
	        Do K = 1, KBM
C	           TENSOR(ID,K) = 0.0
C	           VORTICITY(ID,K) = 0.0
	        Enddo
	     Enddo
	   Endif
c-----   discharge boundary condition
	   If(NUMQBC .NE. 0) Then
	      Do N = 1, NUMQBC
	         ID = IQBC(N)
	         Do K = 1, KBM
C			    TENSOR(ID,K) = 0.0
C			    VORTICITY(ID,K) = 0.0
	         Enddo
	      Enddo
	   Endif
C-----   velocity boundary condition
	   If(NUMVBC .NE. 0) Then
	      Do N = 1, NUMVBC
	         ID = IVBC(N)
	         Do K = 1, KBM
C			    TENSOR(ID,K) = 0.0
C			    VORTICITY(ID,K) = 0.0
	         Enddo
	      Enddo
	   Endif
C===========================================================================C
C         discharge boundary conditions                                     c
c===========================================================================c
	If(NUMQBC .NE. 0) Call BCOND(3)	 
	If(NUMVBC .NE. 0) Call BCOND(2)  
C===========================================================================C
C         update the upwind and downwind direction for every side           c
C===========================================================================C
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&         PRIVATE(I,J,K,IL,IR,UD,N1,N2) 
	Do K = 1, KBM
!$OMP DO
	   Do I = 1, IJE
	      If(CFM(I) .EQ. 1.0) Then
	        IL = INDEX_EDGE(I,K,1)
	        IR = INDEX_EDGE(I,K,2)
	        UD = UN(I,K) * (CXY(IR,1) - CXY(IL,1)) + 
     &			   VN(I,K) * (CXY(IR,2) - CXY(IL,2))
		      If(UD .LT. 0.0) Then
	           INDEX_EDGE(I,K,1) = IR
	           INDEX_EDGE(I,K,2) = IL
	        Endif
	      Endif
	   Enddo
!$OMP END DO NOWAIT
	Enddo
C===========================================================================C
C                       variable on the cell vertex                         c
C===========================================================================C
!$OMP DO
	Do I = 1, IJP
C	   ELFV(I) = 0.0
C	   DV(I) = 0.0
	   Do K = 1, KB
	      UV(I,K) = 0.0
	      VV(I,K) = 0.0
	      WV(I,K) = 0.0
	      QZV(I,K) = 0.0
	      URV(I,K) = 0.0
	      VRV(I,K) = 0.0
	      WRV(I,K) = 0.0
	      PNV(I,K) = 0.0
C	      PTV(I,K) = 0.0
	      TENSORV(I,K) = 0.0
	      KMV(I,K) = 0.0
	      VISV(I,K) = 0.0
	      VORTICITYV(I,K) = 0.0
	      VORXV(I,K) = 0.0
	      VORYV(I,K) = 0.0
	      VORZV(I,K) = 0.0
	      OMIGAV(I,K) = 0.0
	      TKEV(I,K) = 0.0
	      TDISSV(I,K) = 0.0
            CV(I,K) = 0.0
            
            POREV(I,K) = 0.0     

            DEMFORCEXV(I,K) = 0.0   
            DEMFORCEYV(I,K) = 0.0   
            DEMFORCEZV(I,K) = 0.0            
            
	      TKERV(I,K) = 0.0
	      TDISSRV(I,K) = 0.0
            CRV(I,K) = 0.0 
                               
            RRF1V(I,K) = 0.0 
            RRF2V(I,K) = 0.0
            RRF3V(I,K) = 0.0
            Do K1 = 1, 9
               STRESSV(I,K,K1) = 0.0
            Enddo
            
	   Enddo
	Enddo
!$OMP END DO NOWAIT
	DO K = 1, KBM
!$OMP DO
	   Do I = 1, IJP
C	   RTOL = 0.0
C	   IIR = 0
	     Do J = 1, POINT_CELL(I)
	        If(CCM(INDEX_POINT(I,J)) .EQ. 1.0) Then
!---------------Variables of vertex at middle layer--------------------
	  	     UV(I,K) = UV(I,K) + U(INDEX_POINT(I,J),K) * RCOE(I,J) 
		         VV(I,K) = VV(I,K) + V(INDEX_POINT(I,J),K) * RCOE(I,J) 
		         WV(I,K) = WV(I,K) + W(INDEX_POINT(I,J),K) * RCOE(I,J)
		         QZV(I,K) = QZV(I,K)+ QZ(INDEX_POINT(I,J),K)* RCOE(I,J) 
                  
                 VISV(I,K) =VISV(I,K)+VIS(INDEX_POINT(I,J),K)*RCOE(I,J)
                 TKEV(I,K) = TKEV(I,K)+TKE(INDEX_POINT(I,J),K)*RCOE(I,J)
                 TDISSV(I,K) = TDISSV(I,K) +
     &                         TDISS(INDEX_POINT(I,J),K) * RCOE(I,J) 
                  CV(I,K) = CV(I,K)+C(INDEX_POINT(I,J),K)*RCOE(I,J)
                 
		     Endif	       
	     Enddo
	  Enddo
!$OMP END DO
	Enddo
C=============================================================================C
c-----   for saving the date at grid points  ---------------------------------c
C=============================================================================C
	Do K = 2, KBM
!$OMP DO
	   Do I = 1, IJP
	      Do J = 1, POINT_CELL(I)
	         If(CCM(INDEX_POINT(I,J)) .EQ. 1.0) Then
	         
                  URV(I,K) = URV(I,K) + 
     &                       (UR(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                        UR(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J) 
     
	            VRV(I,K) = VRV(I,K)+ 
     &                       (VR(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                        VR(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J)
      
	            WRV(I,K) = WRV(I,K) + 
     &                       (WR(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                        WR(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J) 
     
	            TENSORV(I,K)=TENSORV(I,K)+ 
     &                   (TENSOR(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                    TENSOR(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J)
      
	            KMV(I,K) = KMV(I,K)+ 
     &                       (KM(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                        KM(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J) 
                                  
	            VORTICITYV(I,K) = VORTICITYV(I,K) + 
     &                (VORTICITY(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                 VORTICITY(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J) 
	            
     	            VORXV(I,K) = VORXV(I,K) +  
     &                (VORX(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                 VORX(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J)
     	            
     	            VORYV(I,K) = VORYV(I,K) +
     & 	              (VORY(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                 VORY(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J)
     
     	            VORZV(I,K) = VORZV(I,K) +
     &	              (VORZ(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                 VORZ(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J)

                  PNV(I,K) = PNV(I,K) +
     &                    (PN(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                     PN(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J)
                  
                  TKERV(I,K) = TKERV(I,K) + 
     &                      (TKE(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                       TKE(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J) 
                   
                  TDISSRV(I,K) = TDISSRV(I,K) +
     &                    (TDISS(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                     TDISS(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J) 

                  CRV(I,K) = CRV(I,K) + 
     &                    (C(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                     C(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J) 

                  POREV(I,K) = POREV(I,K) + 
     &                    (PORE(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                     PORE(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J)                  

                  DEMFORCEXV(I,K) = DEMFORCEXV(I,K) + 
     &             (DEMFORCEX(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &              DEMFORCEX(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J) 
                  DEMFORCEYV(I,K) = DEMFORCEYV(I,K) + 
     &             (DEMFORCEY(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &              DEMFORCEY(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J)   
                  DEMFORCEZV(I,K) = DEMFORCEZV(I,K) + 
     &             (DEMFORCEZ(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &              DEMFORCEZ(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J)                    
                  
                  RRF1V(I,K) = RRF1V(I,K) + 
     &                    (RRF1(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                     RRF1(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J) 
                                    
                  RRF2V(I,K) = RRF2V(I,K) + 
     &                    (RRF2(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                     RRF2(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J) 
                  
                  RRF3V(I,K) = RRF3V(I,K) + 
     &                    (RRF3(INDEX_POINT(I,J),K-1) * (Z(K)-ZZ(K)) + 
     &                     RRF3(INDEX_POINT(I,J),K) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J) 
                                    
                  
                  Do K1 = 1, 9
                      STRESSV(I,K,K1) = STRESSV(I,K,K1) +
     &                (STRESS(INDEX_POINT(I,J),K-1,K1) * (Z(K)-ZZ(K)) + 
     &                 STRESS(INDEX_POINT(I,J),K,K1) * (ZZ(K-1)-Z(K))) /
     &                       (ZZ(K-1) - ZZ(K)) * RCOE(I,J) 
                  Enddo
               Endif
             Enddo 
             QCRIV(I,K) = 
     &                (VORTICITYV(I,K) ** 2.-TENSORV(I,K) ** 2.) / 4.
             OMIGAV(I,K) = VORTICITYV(I,K) ** 2. / (1.E-20 + 
     &              VORTICITYV(I,K) **2. + TENSORV(I,K) ** 2. +
     &              0.001*Abs(VORTICITYV(I,K) **2.-TENSORV(I,K) ** 2.))
          Enddo
!$OMP END DO
      Enddo    
!$OMP DO
      Do I = 1, IJP
         If(KBM .GT. 1) Then
            URV(I,1) = 2. * URV(I,2) - URV(I,3)
            VRV(I,1) = 2. * VRV(I,2) - VRV(I,3)
            WRV(I,1) = 2. * WRV(I,2) - WRV(I,3)
            KMV(I,1) = 2. * KMV(I,2) - KMV(I,3)
            
            TENSORV(I,1) = Dmax1(2. * TENSORV(I,2) - TENSORV(I,3),0.0)
            TENSORV(I,KB)=Dmax1(2.*TENSORV(I,KB-1)-TENSORV(I,KB-2),0.0)
            
           VORTICITYV(I,1)=Dmax1(2.*VORTICITYV(I,2)-VORTICITYV(I,3),0.0)
           VORTICITYV(I,KB) =
     &           Dmax1(2. * VORTICITYV(I,KB-1) - VORTICITYV(I,KB-2),0.0)
           
            VORXV(I,1) = 2. * VORXV(I,2) - VORXV(I,3)
            VORXV(I,KB) = 2. * VORXV(I,KB-1) - VORXV(I,KB-2)
            
            VORYV(I,1) = 2. * VORYV(I,2) - VORYV(I,3)
            VORYV(I,KB) = 2. * VORYV(I,KB-1) - VORYV(I,KB-2)
            
            VORZV(I,1) = 2. * VORZV(I,2) - VORZV(I,3)
C            VORZV(I,KB) = 2. * VORZV(I,KB-1) - VORZV(I,KB-2)
            
            PNV(I,1) = 2. * PNV(I,2) - PNV(I,3)
            PNV(I,KB) = 2. * PNV(I,KB-1) - PNV(I,KB-2)
            
            TKERV(I,1) = 2. * TKERV(I,2) - TKERV(I,3)

            QCRIV(I,1) = 2. * QCRIV(I,2) - QCRIV(I,3)
            QCRIV(I,KB) = 2. * QCRIV(I,KB-1) - QCRIV(I,KB-2)
		
		TDISSRV(I,1) = 2. * TDISSRV(I,2) - TDISSRV(I,3)
            TDISSRV(I,KB) = 2. * TDISSRV(I,KB-1) - TDISSRV(I,KB-2) 
            
            CRV(I,1) = 2. * CRV(I,2) - CRV(I,3)
            CRV(I,KB) = 2. * CRV(I,KB-1) - CRV(I,KB-2)
            
            OMIGAV(I,1) = Dmax1(2. * OMIGAV(I,2) - OMIGAV(I,3),0.0)
            OMIGAV(I,KB) = Dmax1(2.*OMIGAV(I,KB-1)-OMIGAV(I,KB-2),0.0)

            POREV(I,1) = Dmin1(Dmax1(2.*POREV(I,2)-POREV(I,3),0.),1.)
            POREV(I,KB) = Dmin1(Dmax1(2.*POREV(I,KB-1)-
     &                  POREV(I,KB-2),0.), 1.)            
 
            DEMFORCEXV(I,1)=2.*DEMFORCEXV(I,2)-DEMFORCEXV(I,3)  
            DEMFORCEXV(I,KB)=2.*DEMFORCEXV(I,KB-1)-DEMFORCEXV(I,KB-2) 
            DEMFORCEYV(I,1)=2.*DEMFORCEYV(I,2)-DEMFORCEYV(I,3)  
            DEMFORCEYV(I,KB)=2.*DEMFORCEYV(I,KB-1)-DEMFORCEYV(I,KB-2) 
            DEMFORCEZV(I,1)=2.*DEMFORCEZV(I,2)-DEMFORCEZV(I,3)  
            DEMFORCEZV(I,KB)=2.*DEMFORCEZV(I,KB-1)-DEMFORCEZV(I,KB-2)            
            
            RRF1V(I,1) = 2. * RRF1V(I,2) - RRF1V(I,3)
            RRF1V(I,KB) = 2. * RRF1V(I,KB-1) - RRF1V(I,KB-2)
            RRF2V(I,1) = 2. * RRF2V(I,2) - RRF2V(I,3)
            RRF2V(I,KB) = 2. * RRF2V(I,KB-1) - RRF2V(I,KB-2)
            RRF3V(I,1) = 2. * RRF3V(I,2) - RRF3V(I,3)
            RRF3V(I,KB) = 2. * RRF3V(I,KB-1) - RRF3V(I,KB-2)
            Do K1 = 1, 9
               STRESSV(I,1,K1) = 2. * STRESSV(I,2,K1) - STRESSV(I,3,K1)
               STRESSV(I,KB,K1) = 
     &                     2. * STRESSV(I,KB-1,K1) - STRESSV(I,KB-2,K1)
            Enddo
         Endif    
      Enddo
!$OMP ENDDO                

!$OMP DO
	Do I = 1, IJE
	   If(CFM(I) .EQ. 0.0 .OR. CFM(I) .EQ. -1.0) Then
	      If(ISLIP .EQ. 0.0) Then
	         N1 = IEND_EDGE(I,1)
	         N2 = IEND_EDGE(I,2)
	         Do K = 1, KBM
	            UV(N1,K) = 0.0
	            UV(N2,K) = 0.0
	            VV(N1,K) = 0.0
	            VV(N2,K) = 0.0
	            QZV(N1,K) = 0.0
	            QZV(N2,K) = 0.0
	            WV(N1,K) = 0.0
	            WV(N2,K) = 0.0
	            URV(N1,K) = 0.0
	            URV(N2,K) = 0.0
	            VRV(N1,K) = 0.0
	            VRV(N2,K) = 0.0
	            WRV(N1,K) = 0.0
	            WRV(N2,K) = 0.0
	            TENSORV(N1,K) = 0.0
	            TENSORV(N2,K) = 0.0
	            KMV(N1,K) = 0.0
	            KMV(N2,K) = 0.0
	            VISV(N1,K) = 0.0
	            VISV(N2,K) = 0.0
	            TKEV(N1,K) = 0.0
	            TKEV(N2,K) = 0.0
                  TKERV(N1,K) = 0.0
	            TKERV(N2,K) = 0.0
                  
                  RRF1V(N1,K) = 0.0
	            RRF1V(N2,K) = 0.0                  
                  RRF2V(N1,K) = 0.0
	            RRF2V(N2,K) = 0.0
                  RRF3V(N1,K) = 0.0
	            RRF3V(N2,K) = 0.0                  
C	            TDISSV(N1,K) = 0.0
C	            TDISSV(N2,K) = 0.0
	         Enddo
	      Endif
	   Endif
	   If(CFM(I) .EQ. -2.0) Then
            N1 = IEND_EDGE(I,1)
	      N2 = IEND_EDGE(I,2)
	      Do K = 1, KBM
	         UV(N1,K) = 0.0
	         UV(N2,K) = 0.0
	         VV(N1,K) = 0.0
	         VV(N2,K) = 0.0
	         QZV(N1,K) = 0.0
	         QZV(N2,K) = 0.0
	         WV(N1,K) = 0.0
	         WV(N2,K) = 0.0
	         URV(N1,K) = 0.0
	         URV(N2,K) = 0.0
	         VRV(N1,K) = 0.0
	         VRV(N2,K) = 0.0
	         WRV(N1,K) = 0.0
	         WRV(N2,K) = 0.0
	         TENSORV(N1,K) = 0.0
	         TENSORV(N2,K) = 0.0
	         KMV(N1,K) = 0.0
	         KMV(N2,K) = 0.0
	         VISV(N1,K) = 0.0
	         VISV(N2,K) = 0.0
	         TKEV(N1,K) = 0.0
	         TKEV(N2,K) = 0.0
               TKERV(N1,K) = 0.0
	         TKERV(N2,K) = 0.0
               
               RRF1V(N1,K) = 0.0
	         RRF1V(N2,K) = 0.0                  
               RRF2V(N1,K) = 0.0
               RRF2V(N2,K) = 0.0
               RRF3V(N1,K) = 0.0
	         RRF3V(N2,K) = 0.0                
C	            TDISSV(N1,K) = 0.0
C	            TDISSV(N2,K) = 0.0
	      Enddo
	  Endif
	Enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
	Return
	End
