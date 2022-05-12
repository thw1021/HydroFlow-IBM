C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (30/04/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c      computes the advective,Coriolis,horizontal dispersive terms in the   #
c      momentum equation of v                                               #
c                                                                           #
c############################################################################
      Subroutine ADVV
      Include './Include/OCERM_INF'
      Include './Include/VORGEN_INF'

	Parameter(Z0WALL=1.E-5,I_METHOD = 1)
	Dimension HQ(IJE,KB),VQ(IJM,KB),VFHYD(IJM)
	Dimension VAR_EDGE(IJE,KB),VAR_T(IJM,-1:KB+1),GRADX(IJM,KB),
     &	         GRADY(IJM,KB),GRADZ(IJM,KB)
      Dimension FYE(IJM),FYH(IJM),TEMP(KB+1)
      Dimension REY_STRESS_VU(NUM_CELL,NUM_VER),
     &           REY_STRESS_VV(NUM_CELL,NUM_VER),
     &           REY_STRESS_VW(NUM_CELL,NUM_VER) 
     	Dimension FLU_SOURCE(NUM_CELL,NUM_VER)

C===========================================================================c
c                  initialiing arrays                                       c
c===========================================================================c
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(I,K) SHARED(VF,HQ,VQ)
	Do K = 1, KB
!$OMP DO	
	   Do I = 1, IJM
            VF(I,K) = 0.0
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
      If(IWENOSCHEME .NE. 0) Then
         Do K = 1, KBM
            Do I = -1, NUM_GHOST, -1
               QGHOST(I,K) = V(INDEX_GHOST(I),K)
            Enddo
         Enddo   
	Endif
C============================================================================c
c            TVD schemes for the calculation of convective fluxes            c
c============================================================================c
C---- Variables on the midpoint of the cell surfaces
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,IL,IR,VW)
!$OMP DO
      Do I = 1, IJM
	   Do K = 1, KBM
	      VAR_T(I,K) = V(I,K)
	   Enddo
	   VAR_T(I,0) = 2. * V(I,1) - V(I,2)
	   VAR_T(I,-1) = VAR_T(I,0)
	   VAR_T(I,KB) = 0.0
	   VAR_T(I,KB+1) = 0.0
	Enddo
!$OMP END DO
!$OMP DO 
      Do I = 1, IJM
         If(CCM(I) .EQ. 1.0) Then
            Do K = 1, KBM
               GRADX(I,K) = 0.0
               GRADY(I,K) = 0.0
               Do J = 1, CELL_POLYGEN(I)
		          If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	               VW = (V(I,K) + V(CELL_SIDE(I,J,2),K)) / 2.
	               GRADX(I,K) = GRADX(I,K) + VW *
     &                         CELL_CUV(I,J,7) * CELL_CUV(I,J,6)
	               GRADY(I,K) = GRADY(I,K) + VW *
     &                         CELL_CUV(I,J,8) * CELL_CUV(I,J,6)
                  Endif
                  If(CFM(CELL_SIDE(I,J,1)) .EQ. 0.0 .OR.
     &               CFM(CELL_SIDE(I,J,1)) .EQ. -1.0) Then
                     If(ISLIP .EQ. 1) Then
                        VW = -(U(I,K)*CELL_CUV(I,J,8) - 
     &                         V(I,K)*CELL_CUV(I,J,7)) * CELL_CUV(I,J,7)
                        GRADX(I,K) = GRADX(I,K) + VW *
     &                               CELL_CUV(I,J,7) * CELL_CUV(I,J,6)
	                  GRADY(I,K) = GRADY(I,K) + VW *
     &                               CELL_CUV(I,J,8) * CELL_CUV(I,J,6)
                     Endif
                  Endif
                  If(CFM(CELL_SIDE(I,J,1)) .EQ. -3.0) Then
                     VW = -(U(I,K)*CELL_CUV(I,J,8) - 
     &                      V(I,K)*CELL_CUV(I,J,7)) * CELL_CUV(I,J,7)
                     GRADX(I,K) = GRADX(I,K) + VW *
     &                            CELL_CUV(I,J,7) * CELL_CUV(I,J,6)
	               GRADY(I,K) = GRADY(I,K) + VW *
     &                            CELL_CUV(I,J,8) * CELL_CUV(I,J,6)
                  Endif
               Enddo
               GRADX(I,K) = GRADX(I,K) / AREA(I)
               GRADY(I,K) = GRADY(I,K) / AREA(I)   
            Enddo
         Endif
      Enddo
!$OMP END DO
!$OMP DO 
      Do I = 1, IJM
         If(CCM(I) .EQ. 1.0) Then
            Do K = 2, KBM
               GRADZ(I,K) = (V(I,K-1) - V(I,K)) / DC(I) / DZZ(K-1)
            Enddo
            GRADZ(I,1) = 0.0
            GRADZ(I,KB) = V(I,KBM) / (0.5 * DZ(KBM) * DC(I))
         Endif
      Enddo
!$OMP END DO
!$OMP END PARALLEL    !#add#, WJ, 2020-3-5 23:10:42 !    
CCC	Call GRAD_XY(VAR_EDGE,GRADX,GRADY)
C	Call GRAD_Z(VAR_UD,GRADZ)
      Call TVDSCHEMEH(HQ,V,GRADX,GRADY,IH_TVD)                     
      Call TVDSCHEMEV(VQ,V,VAR_T,GRADZ,IV_TVD)
c============================================================================c
C                Advection Descrization by 2nd Order TVD                     c
c============================================================================c
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP&         PRIVATE(I,J,K,AAMF,FLUX1,FLUX2,FLUX3,D2,
!$OMP&                 CDWALL,TEMP,ID,IS,CROSS_LENGTH,
!$OMP&                 UMEAN,VMEAN,WMEAN)
      If (ADVECT.EQ.'NON-LINEAR') Then
c----------------------------------------------------------------------------c
c                    horizontal advective terms                              c
c----------------------------------------------------------------------------c
	   Do K = 1, KBM
!$OMP DO	   
	      Do I = 1, IJM
	         If(CCM(I) .EQ. 1.0) Then
	            Do J = 1, CELL_POLYGEN(I)
		           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	                  VF(I,K) = VF(I,K) + 
     &				     HQ(CELL_SIDE(I,J,1),K) * CELL_CUV(I,J,6) *
     &                     PORE_HF(CELL_SIDE(I,J,1),K) *             
     &				     (UN(CELL_SIDE(I,J,1),K) * CELL_CUV(I,J,7) +
     &                    VN(CELL_SIDE(I,J,1),K) * CELL_CUV(I,J,8)) 
	               Endif
	            Enddo
	         Endif
	         VF(I,K) = VF(I,K) * DZ(K)
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
	        Do K = 1, KBM
	           UNEBC = UR(ID,K) * CELL_CUV(ID,IS,7) + 
     & 			       VR(ID,K) * CELL_CUV(ID,IS,8)
	           If(UNEBC .LE. 0.0) Then
	              VF(ID,K) = 0.0
	           Else
	              VF(ID,K) = VF(ID,K) + DZ(K) * CELL_CUV(ID,IS,6) *
     &			             V(ID,K) * UNEBC
	           Endif
	        Enddo
	     Enddo
	   Endif
C-----   astrotidal boundary condition
  	   If(NUMAST .NE. 0) Then
	     Do N = 1, NUMAST
	        ID = IABC(N)
	        IS = IABCINX(N)
	        Do K = 1, KBM
	           UNEBC = UR(ID,K) * CELL_CUV(ID,IS,7) + 
     & 			       VR(ID,K) * CELL_CUV(ID,IS,8)
	           If(UNEBC .LE. 0.0) Then
	              VF(ID,K) = 0.0
	           Else
	              VF(ID,K) = VF(ID,K) + DZ(K) * CELL_CUV(ID,IS,6) *
     &			             V(ID,K) * UNEBC
	           Endif
	        Enddo
	     Enddo
	   Endif
c-----   discharge boundary condition
	   If(NUMQBC .NE. 0) Then
	      Do N = 1, NUMQBC
	         ID = IQBC(N)
	         IS = IQBCINX(N)
	         Do K = 1, KBM
			    VF(ID,K) = VF(ID,K) + DZ(K) * DS(CELL_SIDE(ID,IS,1)) *
     &				  VN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,6) *
     &				  (UN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,7) +
     &                   VN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,8))
	         Enddo
	      Enddo
	   Endif
c-----   velocity boundary condition
	   If(NUMVBC .NE. 0) Then
	      Do N = 1, NUMVBC
	         ID = IVBC(N)
	         IS = IVBCINX(N)
	         Do K = 1, KBM
			    VF(ID,K) = VF(ID,K) + DZ(K) * DS(CELL_SIDE(ID,IS,1)) *
     &				  VN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,6) *
     &				  (UN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,7) +
     &                   VN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,8))
	         Enddo
	      Enddo
	   Endif
!$OMP END MASTER
!$OMP BARRIER
c----------------------------------------------------------------------------c
c                    vertical advective terms                                c
c----------------------------------------------------------------------------c
	   Do K = 1, KBM
!$OMP DO	   
	      Do I = 1, IJM
	         If(CCM(I) .EQ. 1.0) Then
	            VF(I,K) = -VF(I,K) - AREA(I) *
     &			          (VQ(I,K) * W(I,K) * PORE_VF(I,K) - 
     &			           VQ(I,K+1) * W(I,K+1) * PORE_VF(I,K+1))
	         Endif
	      Enddo
!$OMP ENDDO NOWAIT	      
	   Enddo
!$OMP BARRIER
	Endif
c============================================================================c
C                   Horizontal Diffusion Descrization by CS                  c
c============================================================================c
	Do K = 1, KBM
!$OMP DO	
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	        Do J = 1, CELL_POLYGEN(I)
	           FLUX1 = 0.0
	           FLUX2 = 0.0
c                wet cell side
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then

	            AAMF = (AAM(I,K) + AAM(CELL_SIDE(I,J,2),K)) / 2.+ UMOL

		      FLUX1 = (DISCOE(I,J,7) - DISCOE(I,J,2)) * AAMF * 
     &				      (VV(CELL_SIDE(I,J,4),K) - 
     &				       VV(CELL_SIDE(I,J,3),K)) *
     &                      PORE_HF(CELL_SIDE(I,J,1),K)

		          FLUX2 = (DISCOE(I,J,1) - DISCOE(I,J,8)) * AAMF*
     &                      (V(CELL_SIDE(I,J,2),K) - V(I,K)) * 
     &                      PORE_HF(CELL_SIDE(I,J,1),K)

c		          FLUX3 = AAMF * DISCOE(I,J,3) * 
c     &                      (UV(CELL_SIDE(I,J,4),K) - 
c     &				       UV(CELL_SIDE(I,J,3),K))

c		          FLUX4 = -AAMF * DISCOE(I,J,4) * 
c     &                      (U(CELL_SIDE(I,J,2),K) - U(I,K))

				  VF(I,K) = VF(I,K) + (FLUX1 + FLUX2) * DZ(K)
	           Endif
c	           dry side or solid boundary
		       If(CFM(CELL_SIDE(I,J,1)) .EQ. 0.0 .OR.
     &			  CFM(CELL_SIDE(I,J,1)) .EQ. -1.0) Then
c-----   refering to the rough side wall
                    AAMF = AAM(I,K) + UMOL

	              FLUX1 = V(I,K) * DISCOE(I,J,8) * AAMF * PORE(I,K)

c	              FLUX2 = U(I,K) * DISCOE(I,J,4) * AAMF

	              If(Abs(CELL_CUV(I,J,7)) .LE. 1.E-6) Then
	                 FLUX3 = 0.0
	              Else
 	                 If(ISLIP .EQ. 1) Then
	                    FLUX3 = 0.0
	                 Else
C					    D2 = Abs((PXY(CELL_SIDE(I,J,3),1) - CXY(I,1))* 
C     &                             (PXY(CELL_SIDE(I,J,4),2) - CXY(I,2))- 
C     &                             (PXY(CELL_SIDE(I,J,4),1) - CXY(I,1))*
C     &                             (PXY(CELL_SIDE(I,J,3),2)-CXY(I,2)))/ 
C     &                         CELL_CUV(I,J,6)
                        D2 = D2D(I)
c-----   drag coefficient
                         If(VERTMIX .EQ. 'SSTMODEL  ') Then
                             ZSTAR = Dmax1(Sqrt(0.3)*
     *                          Sqrt(TKE(I,K)) * D2 / 1.E-6,15.0)
                             CDWALL = 0.41 * Sqrt(0.3)*
     *                         Sqrt(TKE(I,KBM)+1.E-20)/Log(9.81*ZSTAR)
                         Endif
                         FLUX3 = - DISCOE(I,J,1) * D2 * CDWALL * V(I,K) 
     *                           * PORE(I,K)

C                          CDWALL = 0.4 ** 2. / (Dlog(D2 / Z0WALL)) ** 2. 
C	                    FLUX3 = -DISCOE(I,J,1) * D2 / CELL_CUV(I,J,7)*
C     &                            CDWALL * Abs(V(I,K)*CELL_CUV(I,J,7))*
C     &                            V(I,K) * CELL_CUV(I,J,7) / DC(I)
	                 Endif
	              Endif
				  VF(I,K) = VF(I,K) + (FLUX1 + FLUX3) * DZ(K)
	           Endif
c	           non-slip boundary
      	         If(CFM(CELL_SIDE(I,J,1)) .EQ. -2.0) Then
                    FLUX2 = -(UMOL + AAM(I,K)) * PORE(I,K) *
     &                      (DISCOE(I,J,1)-DISCOE(I,J,8)) * V(I,K)
	              VF(I,K) = VF(I,K) + FLUX2 * DZ(K)
	           Endif
c-----           slip wall 
      	         If(CFM(CELL_SIDE(I,J,1)) .EQ. -3.0) Then
                    AAMF = AAM(I,K) + UMOL
	              FLUX1 = V(I,K) * DISCOE(I,J,8) * AAMF * PORE(I,K)
	              VF(I,K) = VF(I,K) + FLUX1 * DZ(K)
	           Endif
	        Enddo
	      Endif
	   Enddo
!$OMP END DO NOWAIT	   
	Enddo
!$OMP BARRIER
c============================================================================c
c                         source terms treatments                            c
c============================================================================c     
C                   Coriolis Terms and POREMODULE                            c
c----------------------------------------------------------------------------c
	If(CORLIS .EQ. 'INCLUDE') Then
!$OMP DO	
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	        Do K = 1, KBM
	           VF(I,K) = VF(I,K) - COR(I) * PORE(I,K) * U(I,K) * 
     &                     AREA(I) * DZ(K) 
	        Enddo
	      Endif
	   Enddo
!$OMP END DO
      Endif

      If(POREMODULE .EQ. 'INCLUDE' .AND. DEM .EQ. 'NEGLECT') Then       
!$OMP DO	
        Do I =1,IJM
           If(CCM(I).EQ.1.0)Then
              Do K=1,KBM
                 If(PORE(I,K) .NE. 1.0)Then
                    VF(I,K) =  VF(I,K) - V(I,K) * AREA(I) * DZ(K) *
     &                     SQRT(UR(I,K)**2.+VR(I,K)**2.+WR(I,K)**2.)*
     &                     1. / 2. * APU(I,K) * CDC 
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
        Do I = 1 , IJM
         If( CCM(I) .EQ. 1.0) Then
            Do K = 1 ,KBM
                VF(I,K) = VF(I,K) + DEMFORCEY(I,K) * AREA(I) * DZ(K)
            Enddo
         Endif
       Enddo
!$OPM END DO
      Endif
c----------------------------------------------------------------------------c
c                        immersed boundary force                             c
c----------------------------------------------------------------------------c
      If ((MODULEIBM .NE. 'NEGLECT') .AND. (IBFSCHEME.EQ.2)) Then
!$OMP DO
      Do K = 1, KBM
          Do I = 1, IJM
              VF(I,K) = VF(I,K) + BFY(I,K) * DZ(K) * AREA(I)
          Enddo
      Enddo
!$OMP END DO 
      Endif
c----------------------------------------------------------------------------c
C                    Atmaspher Pressure                                      c
c----------------------------------------------------------------------------c
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
            VF(I,1) = VF(I,1) +WVSURF(I) / RMEAN(I,1) 
	   Endif
	Enddo
!$OMP END DO
c----------------------------------------------------------------------------c
C    Explicited part of Dynamic pressure and Hydrostatic pressrure           C
C----------------------------------------------------------------------------C
!$OMP DO
	Do I = 1, IJM
	   If(CCM(I) . EQ. 1.0) Then
	      VFHYD(I) = 0.0
	      If(I_METHOD .EQ. 1) Then
                Do J = 1, CELL_POLYGEN(I)
                   If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
                     VFHYD(I) = VFHYD(I) + CELL_CUV(I,J,6) *
     &                   (ELF(CELL_SIDE(I,J,2)) + ELF(I)) / 2. *
     &                    CELL_CUV(I,J,8)
                   Else
                     VFHYD(I) = VFHYD(I) + CELL_CUV(I,J,6) *
     &                   ELF(I) * CELL_CUV(I,J,8)
                   Endif
	        Enddo
              Endif
	      If(I_METHOD .EQ. 2) Then
                Do J = 1, CELL_POLYGEN(I)
                   If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
                     VFHYD(I) = VFHYD(I) +
     &                    WIY(I,J) * (ELF(CELL_SIDE(I,J,2)) - ELF(I))
                   Endif
	        Enddo
                VFHYD(I) = VFHYD(I) * AREA(I)
              Endif
c             VFHYD(I) = VFHYD(I) * AREA(I)
C	      Do K = 1, KBM
C	         VF(I,K) = VF(I,K) - GRAV * DC(I) * (1.0 - THITA) * DZ(K)*
C     &                   VFHYD 
C	      Enddo
	   Endif
	Enddo
!$OMP END DO 
!$OMP MASTER
	If(NUMEBC .NE. 0) Then
		 Do N = 1, NUMEBC
	      ID = IEBC(N)
	      VFHYD(ID) = 0.0
	      Do J = 1, CELL_POLYGEN(ID)
	         If(CFM(CELL_SIDE(ID,J,1)) .EQ. 1.0) Then
	            VFHYD(ID) = VFHYD(ID) +
     &                    WIY(ID,J) * (ELF(CELL_SIDE(ID,J,2)) - ELF(ID))
	         Endif
	      Enddo
	      VFHYD(ID) = VFHYD(ID) * AREA(ID)
	    Enddo
	 Endif
C-----   astrotidal boundary condition
  	 If(NUMAST .NE. 0) Then
	    Do N = 1, NUMAST
	       ID = IABC(N)
	       VFHYD(ID) = 0.0
	       Do J = 1, CELL_POLYGEN(ID)
	          If(CFM(CELL_SIDE(ID,J,1)) .EQ. 1.0) Then
	             VFHYD(ID) = VFHYD(ID) +
     &                    WIY(ID,J) * (ELF(CELL_SIDE(ID,J,2)) - ELF(ID))
	          Endif
	       Enddo
	       VFHYD(ID) = VFHYD(ID) * AREA(ID)
	    Enddo
	 Endif
c-----   discharge boundary condition
	 If(NUMQBC .NE. 0) Then
	    Do N = 1, NUMQBC
	       ID = IQBC(N)
	       VFHYD(ID) = 0.0
	       Do J = 1, CELL_POLYGEN(ID)
	          If(CFM(CELL_SIDE(ID,J,1)) .EQ. 1.0) Then
	             VFHYD(ID) = VFHYD(ID) +
     &                  WIY(ID,J) * (ELF(CELL_SIDE(ID,J,2)) - ELF(ID))
	          Endif
	       Enddo
	       VFHYD(ID) = VFHYD(ID) * AREA(ID)
	    Enddo
	 Endif
C-----   velocity boundary condition
	 If(NUMVBC .NE. 0) Then
	    Do N = 1, NUMVBC
	       ID = IVBC(N)
	       VFHYD(ID) = 0.0
	       Do J = 1, CELL_POLYGEN(ID)
	         If(CFM(CELL_SIDE(ID,J,1)) .EQ. 1.0) Then
	            VFHYD(ID) = VFHYD(ID) +
     &                    WIY(ID,J) * (ELF(CELL_SIDE(ID,J,2)) - ELF(ID))
	         Endif
	       Enddo
	       VFHYD(ID) = VFHYD(ID) * AREA(ID)
	    Enddo
	 Endif
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
      Do I = 1, IJM
 	   If(CCM(I) . EQ. 1.0) Then
	      Do K = 1, KBM
	         VF(I,K) = VF(I,K) - GRAV * DC(I) * (1.0 - THITA) * DZ(K)*
     &                   VFHYD(I) * PORE(I,K) 
	      Enddo
	   Endif
      Enddo
 !$OMP END DO 
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


               REY_STRESS_VU(I,K) = -(UMEAN * CELL_CUV(ID,IS,7) +
     &                        VMEAN * CELL_CUV(ID,IS,8)) * VDIS(I,K)
C    &                               (UDIS(I,K) * CELL_CUV(ID,IS,8) +
C    &                                VDIS(I,K) * CELL_CUV(ID,IS,7))  
   
               REY_STRESS_VV(I,K) = -(UDIS(I,K) * CELL_CUV(ID,IS,7) +
     &                        VDIS(I,K) * CELL_CUV(ID,IS,8)) * VMEAN
    
               REY_STRESS_VW(I,K) = -(UDIS(I,K) * CELL_CUV(ID,IS,7) +
     &                        VDIS(I,K) * CELL_CUV(ID,IS,8)) * VDIS(I,K)
C     &                               (UDIS(I,K) * CELL_CUV(ID,IS,8) +
C     &                                VDIS(I,K) * CELL_CUV(ID,IS,7))

C               REY_STRESS_VU(I,K) = -VDIS(I,K) * UDIS(I,K)
c    &                              -VMEAN * UDIS(I,K) 
c    &                              -UMEAN * VDIS(I,K) 
C               REY_STRESS_VV(I,K) = -VDIS(I,K) * VDIS(I,K)
C     &                              -VMEAN * VDIS(I,K)
C     &                              -VMEAN * VDIS(I,K)   
C               REY_STRESS_VW(I,K) = -VDIS(I,K) * WDIS(I,K)
C     &                              -VMEAN * WDIS(I,K)
C     &                              -WMEAN * VDIS(I,K)
            Enddo
            CROSS_LENGTH = CELL_CUV(ID,IS,6)
            Do K = 1, NUM_VER - 1
CC              FLU_SOURCE(I,K) = DZ(K) * CROSS_LENGTH * 2. *  
CC     &             (UDIS(I,K) + VDIS(I,K)) + AREA(ID) / DC(ID) * WDIS(I,K)
 
               FLU_SOURCE(I,K) = DZ(K) * CROSS_LENGTH * DC(ID) *
     &                        (REY_STRESS_VU(I,K) + REY_STRESS_VV(I,K) +
     &                         REY_STRESS_VW(I,K))
            Enddo
            Do K = 1, NUM_VER - 1
C               VF(ID,K) = VF(ID,K) + V(ID,K) * FLU_SOURCE(I,K)
               VF(ID,K) = VF(ID,K) + FLU_SOURCE(I,K) 
            Enddo
         Enddo
!$OMP END DO
      Endif 
C!$OMP END MASTER
!$OMP BARRIER
C============================================================================C
C                     Step forward in time                                   C
C============================================================================C
	Do K = 1, KBM
!$OMP DO	
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	        VF(I,K) = V(I,K) * AREA(I) * DZ(K) * PORE(I,K) + DTI * VF(I,K)
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
				VF(I,K) = VF(I,K) - 
     &				      DTI * DIFFCOE * V(I,K) * DZ(K) * AREA(I)
	         Enddo
	      Endif
	      If(CXY(I,1) .GE. XR) Then
	         Do K = 1, KBM
	            DIFFCOE = 8.0 * (CXY(I,1) - XR) / (XRR - XR)
				VF(I,K) = VF(I,K) - 
     &				      DTI * DIFFCOE * V(I,K) * DZ(K) * AREA(I)
	         Enddo
	      Endif
	   Enddo
	Endif
c======================  end subroutine program  ============================C
      Return
      End

