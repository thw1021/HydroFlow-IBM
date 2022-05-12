C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (26/05/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                update the variables, including the velocities,            #
c                masks for wet and dry points, and so on.                   #
c############################################################################
      Subroutine REUV
      Include './Include/OCERM_INF'
	Parameter(IWD=3,TOL=0.5)
CC	Call BCOND(2)
c===========================================================================C
c                       reset velocities                                    C
c===========================================================================C
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,DTMAX,DT,WD,
!$OMP&    DP,UD,VD,DX,DY,VFT,IL,IR,INOR)
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     Do K = 1, KBM
	        UR(I,K) = U(I,K) / (HC(I) + ELF(I))
	        VR(I,K) = V(I,K) / (HC(I) + ELF(I))
	        W(I,K) = W(I,K) / (HC(I) + ELF(I))
	        If(HYDTYPE .EQ. 'NONSTATIC') Then
	           WR(I,K) = (QZ(I,K) + QZ(I,K+1))/ 2. / (HC(I) + ELF(I))
CCC	           qz(I,K) = wr(I,K) * (HC(I) + ELF(I))

	        Endif
	     Enddo
		   If(WFBC .EQ. '      FUN1') Then                      ! RANS module
	        TBX(I) = CBC(I) * RMEAN(I,KBM) * 
     &		   Sqrt(UR(I,KBM) ** 2. + VR(I,KBM) ** 2.) * UR(I,KBM)
	        TBY(I) = CBC(I) * RMEAN(I,KBM) * 
     &		   Sqrt(UR(I,KBM) ** 2. + VR(I,KBM) ** 2.) * VR(I,KBM)
	     Endif
	     If(WFBC .EQ. '      FUN2') Then                      ! LES module
	        TBX(I) = CBC(I) * RMEAN(I,KBM) * UR(I,KBM)
	        TBY(I) = CBC(I) * RMEAN(I,KBM) * VR(I,KBM)
	     Endif
	     If(WFBC .EQ. '      FUN3') Then                      ! no-slip B.C.
	        TBX(I) = RMEAN(I,KBM) * (UMOL) * 
     &			     UR(I,KBM) / (DC(I) * DZ(KBM) * .5)
C      - Z01(I))
	        TBY(I) = RMEAN(I,KBM) * (UMOL) *
     &                 VR(I,KBM) / (DC(I) * DZ(KBM) * .5)
C      - Z01(I))
	     Endif
	   Endif
	Enddo
!$OMP END DO 
!$OMP DO	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     UAVE(I) = 0.0
	     VAVE(I) = 0.0
	     Do K = 1, KBM
	        UAVE(I) = UAVE(I) + UR(I,K) * DZ(K)
	        VAVE(I) = VAVE(I) + VR(I,K) * DZ(K)
	     Enddo
	   Endif
	Enddo
!$OMP END DO 
c===========================================================================C
c                    update water depth at cell centers                     c
c===========================================================================C
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      DC(I) = HC(I) + ELF(I)
	   Else
	      DTMAX = -1.E6
	      WD = 0.0
	      DT = 0.0
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .NE. -1. .AND. 
     &			  CFM(CELL_SIDE(I,J,1)) .NE. -2. .AND.
     &            CFM(CELL_SIDE(I,J,1)) .NE. -3. ) Then
	            If(CCM(CELL_SIDE(I,J,2)) .EQ. 1.0) Then
	               DT = ELF(CELL_SIDE(I,J,2))
	               If(DT .GT. DTMAX) Then
				        DTMAX = DT                   ! for the minimal 
	                  INOR = CELL_SIDE(I,J,2)      ! mark for the cell
				     Endif	     
c	               DT2 = DT2 + ELF(CELL_SIDE(I,J,2))
C	               DT = DT + ELF(CELL_SIDE(I,J,2))
C	               WD = WD + 1.
	            Endif
	         Endif
	      Enddo
  	      DT = DTMAX
	      If(DT .GT. -1.E6) Then
	         DC(I) = HC(I) + DT       ! new water depth on dry cell
	         If(DC(I) .GE. TOL) Then      ! new wetting
c	            ELF(I) = DT2 / WD
c	            DC(I) = HC(I) + ELF(I)
                  ELF(I) = TOL - HC(I)
C                  ELF(I) = AREA(INOR) * (ELF(INOR) - ELF(I)) / 
C     &                     (AREA(INOR) + AREA(I)) + ELF(I)
C                  ELF(INOR) = ELF(I)
                  DC(I) = TOL
C                  DC(INOR) = HC(INOR) + ELF(INOR)
	         Endif
	      Endif
	   Endif
	Enddo
!$OMP END DO	
!$OMP DO	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      If(DC(I) .LE. 0.5 * TOL) CCM(I) = 0.0
	   Else
	      If(DC(I) .GE. TOL) CCM(I) = 1.0
	   Endif
	Enddo
!$OMP END DO
!$OMP BARRIER
c---------------------------------------------------------------------------c
c                elevation boundary cells                                   c
c---------------------------------------------------------------------------c
C	Call BCOND(1)
C	Do N = 1, NUMEBC
C	   ID = IEBC(N)
C	   CCM(ID) = -1.0
C	Enddo
c===========================================================================C
c                    update water depth at cell surfaces                    c
c===========================================================================C
C-----     Define the upwind cells                                          c
!$OMP DO
      Do I = 1, IJE
         If(CFM(I) .EQ. 1.0) Then
           IL = INDEX_EDGE(I,1,1)
           IR = INDEX_EDGE(I,1,2)
           UD = .5 * (UAVE(IL) + UAVE(IR)) * (CXY(IR,1) - CXY(IL,1)) + 
     &			.5 * (VAVE(IL) + VAVE(IR)) * (CXY(IR,2) - CXY(IL,2))
           If(UD .LT. 0.0) Then
              INDEX_EDGE(I,1,1) = IR
              INDEX_EDGE(I,1,2) = IL
           Endif
        Endif
      Enddo
!$OMP END DO
!$OMP DO
	Do I = 1, IJE
	   If(CFM(I) .NE. -1.0 .AND. CFM(I) .NE. -2.0 .AND.
     &	    CFM(I) .NE. -3.0 ) Then
	      ELFM(I) = .5 * 
     &		        (ELF(INDEX_EDGE(I,1,1)) + ELF(INDEX_EDGE(I,1,2)))
	      DS(I) = HS(I) + ELFM(I)
	      If(DS(I) .GT. TOL) Then
	         UD = .5 * (UAVE(INDEX_EDGE(I,1,1)) + 
     &			        UAVE(INDEX_EDGE(I,1,2)))
	         VD = .5 * (VAVE(INDEX_EDGE(I,1,1)) + 
     &			        VAVE(INDEX_EDGE(I,1,2)))
               FR_LOCAL = Sqrt(UD**2.+VD**2.) / Sqrt(9.8*DS(I))
	         If(FR_LOCAL .GE. 1.0) Then
	            DP = 0.5
	         Else
	            DP = 0.5	         
	         Endif
c               DP=0.5                                   ! Averaging bew. cells
C               DP = 1.0                                 ! Upwind scheme
C	         DP = Sqrt(AREA(INDEX_EDGE(I,1,2))) /
C     &			  (Sqrt(AREA(INDEX_EDGE(I,1,1))) +
C     &               Sqrt(AREA(INDEX_EDGE(I,1,2))))
	      Else
	         UD = .5 * (UAVE(INDEX_EDGE(I,1,1)) + 
     &			        UAVE(INDEX_EDGE(I,1,2)))
	         VD = .5 * (VAVE(INDEX_EDGE(I,1,1)) + 
     &			        VAVE(INDEX_EDGE(I,1,2)))
	         DX = (CXY(INDEX_EDGE(I,1,2),1)-CXY(INDEX_EDGE(I,1,1),1))
	         DY = (CXY(INDEX_EDGE(I,1,2),2)-CXY(INDEX_EDGE(I,1,1),2))
	         VFT = UD * DX + VD * DY 
               If(VFT .GT. 0.0) DP=1.0
               If(VFT .LT. 0.0) DP=0.0
               If(VFT .EQ. 0.0 .AND. ELF(INDEX_EDGE(I,1,1)) .GE. 
     &			                   ELF(INDEX_EDGE(I,1,2))) DP=1.0
               If(VFT .EQ. 0.0 .AND. ELF(INDEX_EDGE(I,1,1)) .LT. 
     &			                   ELF(INDEX_EDGE(I,1,2))) DP=0.0
	      Endif
	      ELFM(I) = DP * ELF(INDEX_EDGE(I,1,1)) + 
     &		        (1.0 - DP) * ELF(INDEX_EDGE(I,1,2))

	      If(CCM(INDEX_EDGE(I,1,1)) .EQ. 0.0 .AND.    ! dry cells neighbouring
     &         CCM(INDEX_EDGE(I,1,2)) .EQ. 1.0 ) 
     &         ELFM(I) = ELF(INDEX_EDGE(I,1,2))
	      If(CCM(INDEX_EDGE(I,1,1)) .EQ. 1.0 .AND.    ! dry cells neighbouring
     &         CCM(INDEX_EDGE(I,1,2)) .EQ. 0.0 ) 
     &         ELFM(I) = ELF(INDEX_EDGE(I,1,1))

	      DS(I) = HS(I) + ELFM(I)                     ! new depth ar cell face
	   Else
	      If(INDEX_EDGE(I,1,1) .EQ. -999) 
     &         ELFM(I) = ELF(INDEX_EDGE(I,1,2))
	      If(INDEX_EDGE(I,1,2) .EQ. -999)
     &	       ELFM(I) = ELF(INDEX_EDGE(I,1,1)) 
      	       DS(I) = HS(I) + ELFM(I)
	   Endif
	Enddo
!$OMP END DO NOWAIT
c===========================================================================C
C                redefine the drying and flooding points                    c
c===========================================================================C
!$OMP DO
	Do I = 1, IJE
	   If(CFM(I) .EQ. 1.0) Then
	      If(DS(I) .LE. .5 * TOL) CFM(I) = 0.0
	   Else
	      If(CFM(I) .EQ. 0.0) Then
	         If(DS(I) .GE. TOL) CFM(I) = 1.0
	      Endif
	   Endif
	Enddo
!$OMP END DO 
!$OMP BARRIER	
C===========================================================================C
C                  redefine the cell status                                 c
c===========================================================================c
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 0.0) Then
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .NE. -1.0 .AND. 
     &			CFM(CELL_SIDE(I,J,1)) .NE. -2.0 .AND.
     &          CFM(CELL_SIDE(I,J,1)) .NE. -3.0) 
     &	        CFM(CELL_SIDE(I,J,1)) = 0.0
	      Enddo
	      Do K = 1, KB
	         UR(I,K) = 0.0
	         VR(I,K) = 0.0
	         WR(I,K) = 0.0
	         U(I,K) = 0.0
	         V(I,K) = 0.0
	         W(I,K) = 0.0
	         QZ(I,K) = 0.0
	         PN(I,K) = 0.0
	         VIS(I,K) = 0.0
	      Enddo
	      UAVE(I) = 0.0
	      VAVE(I) = 0.0
	   Endif	
	Enddo	
!$OMP END DO NOWAIT
!$OMP END PARALLEL
c----------------------------------------------------------------------------C
      Return
      End
