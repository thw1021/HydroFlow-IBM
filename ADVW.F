C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (28/07/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c      computes the advective,Coriolis,horizontal dispersive terms in the   #
c      momentum equation of w                                               #
c                                                                           #
c############################################################################
      Subroutine ADVW
      Include './Include/OCERM_INF'
      Include './Include/VORGEN_INF'

	Parameter(Z0WALL=1.E-5)
	Dimension HQ(IJE,KB),VQ(IJM,KB)
	Dimension VAR_EDGE(IJE,KB),VAR_UD(IJM,KB),GRADX(IJM,KB),
     &	         GRADY(IJM,KB),GRADZ(IJM,KB),WW(IJM,KB)
        Dimension REY_STRESS_WU(NUM_CELL,NUM_VER),
     &           REY_STRESS_WV(NUM_CELL,NUM_VER),
     &           REY_STRESS_WW(NUM_CELL,NUM_VER) 
     	Dimension FLU_SOURCE(NUM_CELL,NUM_VER)

      IIII = 0
      IJM_B = IJM_DYN_B + IIII
      IJM_E = IJM_DYN_E + IIII
c===========================================================================c
c                  initialiing arrays                                       c
c===========================================================================c
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(I,K) SHARED(WF,HQ,VQ)
	Do K = 1, KB
!$OMP DO	
	   Do I = 1, IJM
            WF(I,K) = 0.0
	   Enddo
!$OMP END DO 
!$OMP DO 	
	   Do I = 1, IJE
	      HQ(I,K) = 0.0
	   Enddo
!$OMP END DO  
!$OMP DO	
	   Do I = 1, IJM
	      VQ(I,K) = 0.0
	   Enddo
!$OMP END DO
	Enddo
!$OMP END PARALLEL
C============================================================================c
c            TVD schemes for the calculation of convective fluxes            c
c============================================================================c
C---- Variables on the midpoint of the cell surfaces
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,IL,IR,WIS)
c!$OMP DO
c      Do I = 1, IJM
c	   Do K = 1, KB
c	      VAR_UD(I,K) = 0.0
c	   Enddo
c	Enddo
c!$OMP END DO
C      Do K = 1, KBM
C!$OMP DO
C         Do I = 1, IJE
C            VAR_EDGE(I,K) = 0.0
C            IL = INDEX_EDGE(I,K,1)
C            IR = INDEX_EDGE(I,K,2)
C            If(CFM(I) .EQ. 1.0) Then
C               VAR_EDGE(I,K) = (QZ(IL,K) + QZ(IR,K)) / 2.
C            Endif
C            If(CFM(I) .EQ. 0.0 .OR. CFM(I) .EQ. -1.0) Then
C               If(ISLIP .EQ. 1) Then
C                  If(IL .GT. 0) Then
C                     If(CCM(IL) .EQ. 1.0) VAR_EDGE(I,K) = QZ(IL,K)
C                  Endif   
C                  If(IR .GT. 0) Then
C                     If(CCM(IR) .EQ. 1.0) VAR_EDGE(I,K) = QZ(IR,K)
C                  Endif
C               Endif
C            Endif
C            If(CFM(I) .EQ. -3.0) Then
C               If(IL .GT. 0) Then
C                  If(CCM(IL) .EQ. 1.0) VAR_EDGE(I,K) = QZ(IL,K)
C               Endif   
C               If(IR .GT. 0) Then
C                  If(CCM(IR) .EQ. 1.0) VAR_EDGE(I,K) = QZ(IR,K)
C               Endif
C            Endif
C         Enddo
C!$OMP END DO         
C      Enddo
!$OMP DO
      Do I = IJM_B, IJM_E 
         If(CCM(I) .EQ. 1.0) Then
            Do K = 1, KBM
               WW(I,K) = (QZ(I,K) + QZ(I,K+1)) / 2.
            Enddo
         Endif
      Enddo
!$OMP END DO
!$OMP DO 
      Do I = IJM_B, IJM_E
         If(CCM(I) .EQ. 1.0) Then
            Do K = 1, KBM
               GRADX(I,K) = 0.0
               GRADY(I,K) = 0.0
               Do J = 1, CELL_POLYGEN(I)
		  If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	               WIS = (WW(I,K) + WW(CELL_SIDE(I,J,2),K)) / 2.
	               GRADX(I,K) = GRADX(I,K) + WIS *
     &                         CELL_CUV(I,J,7) * CELL_CUV(I,J,6)
	               GRADY(I,K) = GRADY(I,K) + WIS *
     &                         CELL_CUV(I,J,8) * CELL_CUV(I,J,6)
                  Endif
                  If(CFM(CELL_SIDE(I,J,1)) .EQ. 0.0 .OR.
     &               CFM(CELL_SIDE(I,J,1)) .EQ. -1.0) Then
                     If(ISLIP .EQ. 1) Then
                        WIS = WW(I,K)
                        GRADX(I,K) = GRADX(I,K) + WIS *
     &                               CELL_CUV(I,J,7) * CELL_CUV(I,J,6)
	                  GRADY(I,K) = GRADY(I,K) + WIS *
     &                               CELL_CUV(I,J,8) * CELL_CUV(I,J,6)
                     Endif
                  Endif
                  If(CFM(CELL_SIDE(I,J,1)) .EQ. -3.0) Then
                     WIS = WW(I,K)
                     GRADX(I,K) = GRADX(I,K) + WIS *
     &                            CELL_CUV(I,J,7) * CELL_CUV(I,J,6)
	               GRADY(I,K) = GRADY(I,K) + WIS *
     &                            CELL_CUV(I,J,8) * CELL_CUV(I,J,6)
                  Endif
               Enddo
               GRADX(I,K) = GRADX(I,K) / AREA(I)
               GRADY(I,K) = GRADY(I,K) / AREA(I)  
            Enddo
         Endif
      Enddo
!$OMP END DO
c!$OMP DO 
c      Do I = 1, IJM
c         If(CCM(I) .EQ. 1.0) Then
c            Do K = 2, KBM
c               GRADZ(I,K) = QZ(I,K) - QZ(I,K+1)
c            Enddo
c            GRADZ(I,1) = QZ(I,1) - QZ(I,2)
c         Endif
c      Enddo
c!$OMP END DO      

C!$OMP DO	
C	Do I = 1, IJM
C	   If(CCM(I) .EQ. 1.0) Then
C	      Do K = 2, KBM
C	         VAR_UD(I,K) = (QZ(I,K-1) + QZ(I,K)) / 2.0
C	      Enddo
C	      If(KB .GE. 3) Then
C	         VAR_UD(I,1) = 2. * QZ(I,1) - QZ(I,2)
C	         VAR_UD(I,KB) = 0.0
C	      Else
C	         VAR_UD(I,1) = QZ(I,1)
C	      Endif   
C	   Endif
C	Enddo 
C!$OMP END DO
!$OMP END PARALLEL	
C	Call GRAD_XY(VAR_EDGE,GRADX,GRADY)
C	Call GRAD_Z(VAR_UD,GRADZ)
      If(IWENOSCHEME .NE. 0) Then
         Do K = 1, KBM
            Do I = -1, NUM_GHOST, -1
               QGHOST(I,K) = WW(INDEX_GHOST(I),K)
            Enddo
         Enddo   
	Endif
      Call TVDSCHEMEH(HQ,WW,GRADX,GRADY,IH_TVD)
CC      Call TVDSCHEMEV(VQ,QZ,GRADZ,8)
c============================================================================c
C                Advection Descrization by 2nd Order TVD                     c
c============================================================================c
!$OMP PARALLEL DEFAULT(SHARED) SHARED(HQ,VQ) 
!$OMP&         PRIVATE(I,J,K,AAMF,FLUX1,FLUX2,FLUX3,D2,CDWALL,
!$OMP&         ID,CROSS_LENGTH,UMEAN,VMEAN,WMEAN,IS)
      If (ADVECT.EQ.'NON-LINEAR') Then
c----------------------------------------------------------------------------c
c                    horizontal advective terms                              c
c----------------------------------------------------------------------------c
	   Do K = 2, KBM
!$OMP DO	   
	      Do I = IJM_B, IJM_E
	         If(CCM(I) .EQ. 1.0) Then
	            Do J = 1, CELL_POLYGEN(I)
		           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	                  
	                WF(I,K) = WF(I,K) + ( PORE_HF(CELL_SIDE(I,J,1),K-1) +
     &                          PORE_HF(CELL_SIDE(I,J,1),K) ) / 2. *
     &			      (HQ(CELL_SIDE(I,J,1),K) + HQ(CELL_SIDE(I,J,1),K-1)) / 2.*
     &          CELL_CUV(I,J,6) *
     &	     ((UN(CELL_SIDE(I,J,1),K) + UN(CELL_SIDE(I,J,1),K-1)) / 2. *
     &         CELL_CUV(I,J,7) +
     &        (VN(CELL_SIDE(I,J,1),K) + VN(CELL_SIDE(I,J,1),K-1)) / 2. *
     &         CELL_CUV(I,J,8))
	               Endif
	            Enddo
	         Endif
	         WF(I,K) = WF(I,K) * DZZ(K-1)
	      Enddo
!$OMP END DO NOWAIT	      
	   Enddo   
!$OMP BARRIER
c============================================================================c
c                    open boundary treatments                                c
c============================================================================c	   
C-----   elevation boundary condition
!$OMP MASTER
	   If(NUMEBC .NE. 0) Then
		 Do N = 1, NUMEBC
	        ID = IEBC(N)
	        IS = IEBCINX(N)
	        Do K = 2, KBM
	           UNEBC = (UR(ID,K) + UR(ID,K-1)) / 2. * CELL_CUV(ID,IS,7) + 
     & 			         (VR(ID,K) + VR(ID,K-1)) / 2. * CELL_CUV(ID,IS,8)
	           If(UNEBC .GT. 0.0) Then
C	              WF(ID,K) = 0.0
C	           Else
	              WF(ID,K) = WF(ID,K) + DZZ(K-1) * CELL_CUV(ID,IS,6) *
     &			             QZ(ID,K) * UNEBC
	           Endif
c	           UF(ID,K) = 0.0
	        Enddo
	     Enddo
	   Endif
C-----   astrotidal boundary condition
  	   If(NUMAST .NE. 0) Then
	     Do N = 1, NUMAST
	        ID = IABC(N)
	        IS = IABCINX(N)
	        Do K = 2, KBM
	           UNEBC = (UR(ID,K) + UR(ID,K-1)) / 2. * CELL_CUV(ID,IS,7) + 
     & 			         (VR(ID,K) + VR(ID,K-1)) / 2. * CELL_CUV(ID,IS,8)
	           If(UNEBC .GT. 0.0) Then
C	              WF(ID,K) = 0.0
C	           Else
	              WF(ID,K) = WF(ID,K) + DZZ(K-1) * CELL_CUV(ID,IS,6) *
     &			             QZ(ID,K) * UNEBC
	           Endif
c	           UF(ID,K) = 0.0
	        Enddo
	     Enddo
	   Endif
c-----   discharge boundary condition
c	   If(NUMQBC .NE. 0) Then
c	      Do N = 1, NUMQBC
c	         ID = IQBC(N)
c	         IS = IQBCINX(N)
c	         Do K = 1, KBM
c	            WF(ID,K) = 0.0
C			    WF(ID,K) = WF(ID,K) + DZ(K) * DS(CELL_SIDE(ID,IS,1)) *
C     &				  UN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,6) *
C     &				  (UN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,7) +
C     &                   VN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,8))
c	         Enddo
c	      Enddo
c	   Endif
C-----   velocity boundary condition
C	   If(NUMVBC .NE. 0) Then
C	      Do N = 1, NUMVBC
C	         ID = IVBC(N)
C	         IS = IVBCINX(N)
C	         Do K = 1, KBM
c	            WF(ID,K) = 0.0
C			    WF(ID,K) = WF(ID,K) + DZ(K) * DS(CELL_SIDE(ID,IS,1)) *
C     &				  UN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,6) *
C     &				  (UN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,7) +
C     &                   VN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,8))
C	         Enddo
C	      Enddo
C	   Endif

!$OMP END MASTER
!$OMP BARRIER
c----------------------------------------------------------------------------c
c                    vertical advective terms                                c
c----------------------------------------------------------------------------c
	   Do K = 2, KBM
!$OMP DO	   
	      Do I = IJM_B, IJM_E
	         If(CCM(I) .EQ. 1.0) Then
	            WF(I,K) = -WF(I,K) - AREA(I) * ( PORE(I,K-1) *
     &		         (QZ(I,K)+QZ(I,K-1))/2. * (W(I,K)+W(I,K-1))/2. -
     &                  PORE(I,K) * (QZ(I,K)+QZ(I,K+1))/2. * 
     &                  (W(I,K)+W(I,K+1))/2.)
	         Endif
	      Enddo
!$OMP END DO NOWAIT	      
	   Enddo
!$OMP BARRIER
	Endif
c============================================================================c
C                   Horizontal Diffusion Descrization by CS                  c
c============================================================================c
	Do K = 2, KBM
!$OMP DO	
	   Do I = IJM_B, IJM_E
	      If(CCM(I) .EQ. 1.0) Then
	        Do J = 1, CELL_POLYGEN(I)
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	             
	              FLUX1 = 0.0
	              FLUX2 = 0.0

	              AAMF = UMOL + (AAM(I,K) + AAM(CELL_SIDE(I,J,2),K) +
     &	                    AAM(I,K-1)+AAM(CELL_SIDE(I,J,2),K-1)) / 4. 

	              FLUX1 = (DISCOE(I,J,1) - DISCOE(I,J,8)) * AAMF *
     &                      (QZ(CELL_SIDE(I,J,2),K) - QZ(I,K)) *
     &                      ( PORE_HF(CELL_SIDE(I,J,1),K-1) +
     &                       PORE_HF(CELL_SIDE(I,J,1),K) ) / 2.

	              FLUX2 = (DISCOE(I,J,7) - DISCOE(I,J,2)) * AAMF *
     &                      (QZV(CELL_SIDE(I,J,4),K) - 
     &				       QZV(CELL_SIDE(I,J,3),K)) *
     &                      ( PORE_HF(CELL_SIDE(I,J,1),K-1) +
     &                       PORE_HF(CELL_SIDE(I,J,1),K) ) / 2.

				    WF(I,K) = WF(I,K) + (FLUX1 + FLUX2) * DZZ(K-1)
	           Endif
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 0.0 .OR.     
     &			    CFM(CELL_SIDE(I,J,1)) .EQ. -1.0) Then
C-----   refering to the rough side wall

                    If(ISLIP .EQ. 1) Then			   
	                 FLUX1 = 0.0
		            Else
c-----   distance to the wall
C                       D2 = Abs((PXY(CELL_SIDE(I,J,3),1) - CXY(I,1)) * 
C     &                       (PXY(CELL_SIDE(I,J,4),2) - CXY(I,2)) - 
C     &                       (PXY(CELL_SIDE(I,J,4),1) - CXY(I,1)) *
C     &                       (PXY(CELL_SIDE(I,J,3),2) - CXY(I,2))) / 
C     &                      CELL_CUV(I,J,6)	
                        D2 = D2D(I)
c-----   drag coefficient
                        If(VERTMIX .EQ. 'SSTMODEL  ') Then
                             ZSTAR = Dmax1(Sqrt(0.3)*
     *                               Sqrt(TKE(I,K)) * D2 / 1.E-6,15.0)
                             CDWALL = 0.41 * Sqrt(0.3)*
     *                         Sqrt(TKE(I,KBM)+1.E-20)/Log(9.81*ZSTAR)
                         Endif
                         FLUX1 = -CDWALL * QZ(I,K) * D2 *
     &                         (DISCOE(I,J,1) - DISCOE(I,J,8)) *
     &                          PORE(I,K)

C                       CDWALL = 0.4 ** 2. / (Dlog(D2 / Z0WALL)) ** 2.
C	                 FLUX1 = -CDWALL * Abs(QZ(I,K)) * QZ(I,K) * D2 *
C     &                         (DISCOE(I,J,1) - DISCOE(I,J,8)) / DC(I)
	              Endif
				    WF(I,K) = WF(I,K) + FLUX1 * DZZ(K-1)
				 Endif  
c-----           non-slip wall				  
			     If(CFM(CELL_SIDE(I,J,1)) .EQ. -2.0) Then  
			        FLUX1 = -(UMOL + (AAM(I,K)+AAM(I,K-1))/2.) * 
     &   			        (DISCOE(I,J,1)-DISCOE(I,J,8)) * 
     &				         QZ(I,K) * PORE(J,K)
                    WF(I,K) = WF(I,K) + FLUX1 * DZZ(K-1)
			     Endif
c-----           slip wall 
			     If(CFM(CELL_SIDE(I,J,1)) .EQ. -3.0) Then  
			        FLUX1 = 0.0
                    WF(I,K) = WF(I,K) + FLUX1 * DZZ(K-1) * PORE(J,K)
			     Endif
	        Enddo
	      Endif
	   Enddo
!$OMP END DO  
	Enddo
!$OMP BARRIER
c============================================================================c
c                    open boundary treatments                                c
c============================================================================c	   
c	Do N = 1, NUMEBC
c	   ID = IEBC(N)
c	   Do K = 1, KBM
C	      WF(ID,K) = 0.0
c	   Enddo
c	Enddo
c	Do N = 1, NUMAST
c	   ID = IABC(N)
c	   Do K = 1, KBM
C	      WF(ID,K) = 0.0
c	   Enddo
c	Enddo
c============================================================================c
c                         source terms treatments                            c
c============================================================================c     
C                   Coriolis Terms And  POREMODULE                           c
c----------------------------------------------------------------------------c
      If(POREMODULE .EQ. 'INCLUDE' .AND. DEM .EQ. 'NEGLECT') Then     
!$OMP DO
        Do I =2,IJM
           If(CCM(I).EQ.1.0)Then
              Do K=1,KBM
                 If(PORE(I,K) .NE. 1.0)Then
         !               WF(I,K) = WF(I,K) - QZ(I,K) * AREA(I) * DZ(K) *
         !&                     SQRT(UR(I,K)**2.+VR(I,K)**2.+WR(I,K)**2.)*
         !&                     1. / 2. * APU(I,K) * CDC 
                    WF(I,K) = WF(I,K) - QZ(I,K) * AREA(I) * DZZ(K-1) *
     &                        SQRT( ( (UR(I,K) + UR(I,K-1))/2. )**2. +
     &                              ( (VR(I,K) + VR(I,K-1))/2. )**2. +
     &                              ( (WR(I,K) + WR(I,K-1))/2. )**2. ) *
     &                        1. / 2. * APU(I,K) * CDC                   
                 Endif
              Enddo
            Endif
        Enddo
!$OMP END DO 
	Endif

c----------------------------------------------------------------------------c
C                    DEMFORCE                                                c
c----------------------------------------------------------------------------c
      If(DEM .NE. 'NEGLECT' .AND. I_PTF .NE. 0) Then
!$OMP DO
        Do I = IJM_B, IJM_E
         If( CCM(I) .EQ. 1.0) Then
            Do K = 2 ,KBM
                WF(I,K) = WF(I,K) + ( DEMFORCEZ(I,K) + 
     &                    DEMFORCEZ(I,K-1) )/2. * AREA(I) * DZZ(K-1)  
            Enddo
         Endif
        Enddo
!$OPM END DO 
      Endif
c----------------------------------------------------------------------------c
C                    Atmaspher Pressure and wind stress                      c
c----------------------------------------------------------------------------c

c----------------------------------------------------------------------------c
C    Explicited part of Dynamic pressure and Hydrostatic pressrure           C
C----------------------------------------------------------------------------C
      If (HYDTYPE .EQ. 'NONSTATIC1') Then
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         Do K = 2, KBM
	             GRADZ(I,K) = (PN(I,K-1) - PN(I,K)) / DZZ(K-1)
	         Enddo
	         GRADZ(I,1) = (0.0 - PN(I,1)) / (.5 * DZ(1))
               GRADZ(I,KB) = 0.0
	      Endif
	   Enddo
!$OMP END DO
C	   Do I = 1, IJM
C	      If(CCM(I) . EQ. 1.0) Then
C	         Do K = 2, KBM - 1
C	            WF(I,K) = WF(I,K) - (.5 * (PN(I,K-1) + PN(I,K)) -
C    &			                  .5 * (PN(I,K) + PN(I,K+1))) /
C     &                      RMEAN(I,K) * AREA(I)
C	         Enddo
C	         WF(I,1) = WF(I,1) + .5 * (PN(I,1) + PN(I,2)) /
C     &                   RMEAN(I,1) * AREA(I)
C	         WF(I,KBM) = WF(I,KBM) - (PN(I,KBM-1) - PN(I,KBM)) / 
C     &		             RMEAN(I,KBM) * AREA(I)
c		     WF(I,1) = WF(I,1) - (PN(I,1) - PN(I,2)) /
c     &                   RMEAN(I,1) * AREA(I)

c	         WF(I,KBM) = WF(I,KBM) - (PN(I,KBM-1) - PN(I,KBM)) / 2. / 
c     &			         RMEAN(I,KBM) * AREA(I)
C	      Endif
C	   Enddo
!$OMP DO
         Do I = 1, IJM
            If(CCM(I) .EQ. 1.0) Then
	         Do K = 1, KBM - 1
	            WF(I,K) = WF(I,K) - AREA(I) * DZ(K) * 
     &                      PORE_VF(I,K) / RMEAN(I,K) *
     &	                    (QZ_INTER(K,1) * GRADZ(I,K) + 
     &                       QZ_INTER(K,2) * GRADZ(I,K+1) +
     &                       QZ_INTER(K,3) * GRADZ(I,K+2))
	         Enddo
	         WF(I,KBM) = WF(I,KBM) - AREA(I) * DZ(KBM) *
     &                      PORE_VF(I,KB) / RMEAN(I,KBM) *
     &      	           (QZ_INTER(KBM,1) * GRADZ(I,KBM-1) + 
     &                      QZ_INTER(KBM,2) * GRADZ(I,KBM) +
     &                      QZ_INTER(KBM,3) * GRADZ(I,KB))
            Endif
         Enddo   
!$OMP END DO
      Endif      
C----------------------------------------------------------------------------C
c                    wave radiation stress                                   C
C----------------------------------------------------------------------------C
	If(WAVEDYN.NE.'NEGLECT') Then

	Endif
C----------------------------------------------------------------------------C
C                   Baroclinic pressure term                                 C
C----------------------------------------------------------------------------C
      If(TOR.EQ.'BAROCLINIC') Then

	Endif
C----------------------------------------------------------------------------C
C                   Synthetize eddy in ZDES model                            c
C----------------------------------------------------------------------------C
!$OMP BARRIER
C!$OMP MASTER
      If(DES .EQ. 'SAZDES ') Then 
!$OMP DO 
         Do I = 1, NUM_CELL
            ID = ID_CELL(I)
            IS = ID_CELL_EDGE(I) 
            Do K = 1, NUM_VER

               UMEAN = 0.0
               VMEAN = 0.0
               WMEAN = 0.0
               Do J = 1, CELL_POLYGEN(ID)
                  UMEAN = UMEAN + UR(CELL_SIDE(ID,J,2),K)
                  VMEAN = VMEAN + VR(CELL_SIDE(ID,J,2),K)
                  WMEAN = WMEAN + WR(CELL_SIDE(ID,J,2),K)
               Enddo
               UMEAN = UMEAN + UR(ID,K)
               VMEAN = VMEAN + VR(ID,K)
               WMEAN = WMEAN + WR(ID,K)
               UMEAN = UMEAN / (1. + CELL_POLYGEN(ID))
               VMEAN = VMEAN / (1. + CELL_POLYGEN(ID))
               WMEAN = WMEAN / (1. + CELL_POLYGEN(ID))

               REY_STRESS_WU(I,K) = -(UMEAN * CELL_CUV(ID,IS,7) +
     &                                VMEAN * CELL_CUV(ID,IS,8)) * 
     &                               WDIS(I,K)   
   
               REY_STRESS_WV(I,K) = -(UDIS(I,K) * CELL_CUV(ID,IS,7) +
     &                            VDIS(I,K) * CELL_CUV(ID,IS,8)) * WMEAN
    
               REY_STRESS_WW(I,K) = -(UDIS(I,K) * CELL_CUV(ID,IS,7) +
     &                                VDIS(I,K) * CELL_CUV(ID,IS,8)) * 
     &                               WDIS(I,K)

C               REY_STRESS_WU(I,K) = -WDIS(I,K) * UDIS(I,K)
c     &                              -WMEAN * UDIS(I,K) 
c     &                              -UMEAN * WDIS(I,K)    
C               REY_STRESS_WV(I,K) = -WDIS(I,K) * VDIS(I,K)
C     &                              -WMEAN * VDIS(I,K)
C     &                              -VMEAN * WDIS(I,K) 
C               REY_STRESS_WW(I,K) = -WDIS(I,K) * WDIS(I,K)
C     &                              -WMEAN * WDIS(I,K)
C     &                              -WMEAN * WDIS(I,K)
            Enddo
            CROSS_LENGTH = CELL_CUV(ID,IS,6)
            Do K = 1, NUM_VER - 1
CC               FLU_SOURCE(I,K) = DZ(K) * CROSS_LENGTH * 2. *  
CC     &             (UDIS(I,K) + VDIS(I,K)) +
CC     &              AREA(ID) / DC(ID) * WDIS(I,K) 
               FLU_SOURCE(I,K) = DZ(K) * CROSS_LENGTH * DC(ID) *
     &                    (REY_STRESS_WU(I,K) + REY_STRESS_WV(I,K) +
     &                     REY_STRESS_WW(I,K))
            Enddo
            Do K = 1, NUM_VER - 1
C                WF(ID,K) = WF(ID,K) + QZ(ID,K) * FLU_SOURCE(I,K)
               WF(ID,K) = WF(ID,K) + FLU_SOURCE(I,K) 
             Enddo
          Enddo  

C         Do J = 1, NUM_CELL
C            ID = ID_CELL(J)
C            Do K = 1, NUM_VER
C            ZWALL = DC(ID_CELL(J)) * (1. + ZZ(K))
C            If(ZWALL .GE. DWINTERFACE) Then
C               REYU = 0.0
C               REYV = 0.0
C               REYW = 0.0
C               Do I = 1, CELL_POLYGEN(ID)
C                  If(CFM(CELL_SIDE(ID,I,1)) .EQ. 1.0) Then
C	               REYU = REYU -
C     &                      WIX(ID,I) * REY_STRESS_UU(J,K)
C                     REYV = REYV - 
C     &                      WIY(ID,I) * REY_STRESS_UV(J,K)
C                     REYW = (REY_STRESS_UW(J,K) - REY_STRESS_UW(J,K+1))/
C     &                       DZZ(K) / DC(ID)               

C	            Endif
C	         Enddo   
C               WF(ID,K) = WF(ID,K) - DC(ID) * 
C     &                   (REYU + REYV + REYW) * AREA(ID) * DZ(K)
C            Endif
C            Enddo
C         Enddo
!$OMP END DO 
      Endif 
C!$OMP END MASTER
!$OMP BARRIER
C============================================================================C
C                     Step forward in time                                   C
C============================================================================C
	Do K = 2, KBM
!$OMP DO	
	   Do I = IJM_B, IJM_E
	      If(CCM(I) .EQ. 1.0) Then
	        WF(I,K) = QZ(I,K) * AREA(I) * DZZ(K-1) *  
     &                    PORE_VF(I,K) + DTI * WF(I,K)
	      Endif
         Enddo
!$OMP END DO        
	Enddo
!$OMP END PARALLEL	
C============================================================================C
C            Sponge layer setting for numerical wave flume                   C
C============================================================================C
	If(IFLUME .EQ. 1) Then
	   XL = -1.0
	   XR = 26.0
         XLL = -5.
         XRR = 30.0 	  
	   Do I = 1, IJM
	      If(CXY(I,1) .LE. XL) Then
	         Do K = 1, KBM
	            DIFFCOE = 10.0 * (XL - CXY(I,1)) / (XL - XLL)
				  WF(I,K) = WF(I,K) - 
     &				        DTI * DIFFCOE * QZ(I,K) * DZ(K) * AREA(I)
	         Enddo
	      Endif
	      If(CXY(I,1) .GE. XR) Then
	         Do K = 1, KBM
	            DIFFCOE = 10.0 * (CXY(I,1) - XR) / (XRR - XR)
				  WF(I,K) = WF(I,K) - 
     &				        DTI * DIFFCOE * QZ(I,K) * DZ(K) * AREA(I)
	         Enddo
	      Endif
	   Enddo
	Endif
c======================  end subroutine program  ============================C
      Return
      End
