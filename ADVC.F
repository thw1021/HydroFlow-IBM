C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (30/04/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c      computes the advective,horizontal dispersive terms in the mass       #
c      transportation                                                       #
c                                                                           #
c############################################################################
      Subroutine ADVC
      Include './Include/OCERM_INF'
	Dimension HQ(IJE,KB),VQ(IJM,KB)
	Dimension VAR_EDGE(IJE,KB),VAR_UD(IJM,KB),GRADX(IJM,KB),
     &	         GRADY(IJM,KB),GRADZ(IJM,KB),VAR_T(IJM,-1:KB+1)
c===========================================================================c
c                  initialiing arrays                                       c
c===========================================================================c
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,K)
	Do K = 1, KB
!$OMP DO	
	   Do I = 1, IJM
            CF(I,K) = 0.0
	   Enddo
!$OMP END DO NOWAIT	
!$OMP DO	
	   Do I = 1, IJE
	      HQ(I,K) = 0.0
	   Enddo
!$OMP END DO NOWAIT	
!$OMP DO	
	   Do I = 1, IJM
	      VQ(I,K) = 0.0
	   Enddo
!$OMP END DO NOWAIT	
	Enddo
!$OMP END PARALLEL	
      If(IWENOSCHEME .NE. 0) Then
         Do K = 1, KBM
            Do I = -1, NUM_GHOST, -1
               QGHOST(I,K) = C(INDEX_GHOST(I),K)
            Enddo
         Enddo   
	Endif
C============================================================================c
c            TVD schemes for the calculation of convective fluxes            c
c============================================================================c
C---- Variables on the midpoint of the cell surfaces
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,IL,IR)
!$OMP DO
      Do I = 1, IJM
	   Do K = 1, KB
	      VAR_UD(I,K) = 0.0
	   Enddo
	Enddo
!$OMP END DO
!$OMP DO
      Do I = 1, IJM
	   Do K = 1, KBM
	      VAR_T(I,K) = C(I,K)
	   Enddo
	   VAR_T(I,0) = 2. * C(I,1) - C(I,2)
	   VAR_T(I,-1) = VAR_T(I,0)
	   VAR_T(I,KB) = C(I,KBM)
	   VAR_T(I,KB+1) = C(I,KBM)
	Enddo
!$OMP END DO
      Do K = 1, KBM
!$OMP DO
         Do I = 1, IJE
            VAR_EDGE(I,K) = 0.0
            IL = INDEX_EDGE(I,K,1)
            IR = INDEX_EDGE(I,K,2)
            If(CFM(I) .EQ. 1.0) Then
               VAR_EDGE(I,K) = (C(IL,K) + C(IR,K)) / 2.
            Else
               If(IL .GT. 0) Then
                  If(CCM(IL) .EQ. 1.0) VAR_EDGE(I,K) = C(IL,K)
               Endif   
               If(IR .GT. 0) Then
                  If(CCM(IR) .EQ. 1.0) VAR_EDGE(I,K) = C(IR,K)
               Endif
            Endif
         Enddo
!$OMP END DO         
      Enddo
C!$OMP DO	
C	Do I = 1, IJE
C	   If(CFM(I) .EQ. 0.0 .OR. CFM(I) .EQ. -1.0) Then
C	      If(ISLIP .EQ. 0.0) Then
C	         N1 = IEND_EDGE(I,1)
C	         N2 = IEND_EDGE(I,2)
C	         Do K = 1, KBM
C	            VAR_VERTEX(N1,K) = 0.0
C	            VAR_VERTEX(N2,K) = 0.0
C	         Enddo
C	      Endif
C	   Endif
C	   If(CFM(I) .EQ. -2.0) Then
C            N1 = IEND_EDGE(I,1)
C	      N2 = IEND_EDGE(I,2)
C	      Do K = 1, KBM
C	         VAR_VERTEX(N1,K) = 0.0
C	         VAR_VERTEX(N2,K) = 0.0
C	      Enddo
C        Endif
C	Enddo
C!$OMP END DO
!$OMP DO	
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	      Do K = 2, KBM
	         VAR_UD(I,K) = (C(I,K-1) + C(I,K)) / 2.0
	      Enddo
	      If(KB .GE. 3) Then
               VAR_UD(I,1) = 2. * C(I,1) - VAR_UD(I,2)    !#surface variable#, WangJian, 2020-3-8 14:23:40 !    
	         VAR_UD(I,KB) = 2. * C(I,KBM) - VAR_UD(I,KBM)
	      Else
	         VAR_UD(I,1) = C(I,1)
	         VAR_UD(I,2) = C(I,1)
	      Endif   
C--------------Begin--------2019-12-13 11:19:53, added by WangJian------------------c            
            Do K = 2, KBM
               GRADZ(I,K) = (C(I,K-1) - C(I,K))/DC(I)/DZZ(K-1)
            Enddo
               GRADZ(I,1) = 0.0
               GRADZ(I,KB) = 0.0  
C--------------End--------2019-12-13 11:19:58, added by WangJian------------------c                
	   Endif
	Enddo 
!$OMP END DO
!$OMP END PARALLEL	
	Call GRAD_XY(VAR_EDGE,GRADX,GRADY)
!	Call GRAD_Z(VAR_UD,GRADZ)
      Call TVDSCHEMEH(HQ,C,GRADX,GRADY,8)
      Call TVDSCHEMEV(VQ,C,VAR_T,GRADZ,8)
c============================================================================c
C                Advection Descrization by 2nd Order TVD                     c
c============================================================================c
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP&         PRIVATE(I,J,K,AAMF,FLUX1,FLUX2,ID,IS,ISS)
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
					  CF(I,K) = CF(I,K) + DZ(K) * 
     &				     HQ(CELL_SIDE(I,J,1),K) * CELL_CUV(I,J,6) *
     &				     (UN(CELL_SIDE(I,J,1),K) * CELL_CUV(I,J,7) +
     &                      VN(CELL_SIDE(I,J,1),K) * CELL_CUV(I,J,8))
	               Endif
	            Enddo
	         Endif
	      Enddo
!$OMP END DO NOWAIT	      
	   Enddo
!$OMP BARRIER	   
c----------------------------------------------------------------------------c
c                    vertical advective terms                                c
c----------------------------------------------------------------------------c
	   Do K = 1, KBM
!$OMP DO	   
	      Do I = 1, IJM
	         If(CCM(I) .EQ. 1.0) Then
	            CF(I,K) = -CF(I,K) - AREA(I) *
     &			          (VQ(I,K) * W(I,K) - VQ(I,K+1) * W(I,K+1))
	         Endif
	      Enddo
!$OMP END DO NOWAIT	      
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
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then

	            AAMF = (AAM(I,K) + AAM(CELL_SIDE(I,J,2),K)) / 2. + UMOL
	 
		          FLUX1 = (DISCOE(I,J,1) - DISCOE(I,J,8))* AAMF * 
     &				      (C(CELL_SIDE(I,J,2),K) - C(I,K))

		          FLUX2 = (DISCOE(I,J,7) - DISCOE(I,J,2)) * AAMF*
     &                      (CV(CELL_SIDE(I,J,4),K) - 
     &                       CV(CELL_SIDE(I,J,3),K))

				  CF(I,K) = CF(I,K) + (FLUX1 + FLUX2) * DZ(K)
	           Endif
	        Enddo
	      Endif
	   Enddo
!$OMP END DO NOWAIT	   
	Enddo
!$OMP BARRIER	
c============================================================================c
c                    open boundary treatments                                c
c============================================================================c	   
!$OMP MASTER
C-----   elevation boundary condition
	   If(NUMEBC .NE. 0) Then
		   Do N = 1, NUMEBC
	        ID = IEBC(N)
	        IS = IEBCINX(N)
	        Do K = 1, KBM
	           UNEBC = UR(ID,K) * CELL_CUV(ID,IS,7) + 
     & 			         VR(ID,K) * CELL_CUV(ID,IS,8)
			     If(UNEBC .GT. 0.0) Then
	              CF(ID,K) = CF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &			               C(ID,K) * UNEBC
	           Else
                    CF(ID,K) = CF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &			               UNEBC * CEBC(N)
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
	           UNAST = UR(ID,K) * CELL_CUV(ID,IS,7) + 
     & 			       VR(ID,K) * CELL_CUV(ID,IS,8)
			     If(UNAST .GT. 0.0) Then
	              VISF(ID,K) = VISF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &			             C(ID,K) * UNAST
	           Else
                    VISF(ID,K) = VISF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &			             UNAST * CAST(N)
	           Endif
	         Enddo
	     Enddo
	   Endif
c-----   discharge boundary condition
	   If(NUMQBC .NE. 0) Then
	      Call BCOND(3)
	      Do N = 1, NUMQBC
	         ID = IQBC(N)
	         IS = IQBCINX(N)
	         ISS = CELL_SIDE(ID,IS,1)
	         Do K = 1, KBM
 	            UNQBC = UN(ISS,K) * CELL_CUV(ID,IS,7) + 
     & 			          VN(ISS,K) * CELL_CUV(ID,IS,8)
	            If(UNQBC .GT. 0.0) Then
	              CF(ID,K) = CF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &			                 C(ID,K) * UNQBC
	            Else
	              CF(ID,K) = CF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &                        UNQBC * CQBC(N,K) 
	            Endif
	         Enddo
	      Enddo
	   Endif
C-----   velocity boundary condition
	   If(NUMVBC .NE. 0) Then
	      Do N = 1, NUMVBC
	         ID = IVBC(N)
	         IS = IVBCINX(N)
	         Do K = 1, KBM
 	            UNVBC = UN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,7) + 
     & 			        VN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,8)
	            If(UNVBC .GT. 0.0) Then
	              CF(ID,K) = CF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6) *
     &			               C(ID,K) * UNVBC
	            Else
	              CF(ID,K) = CF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6) *
     &                         CVBC(N,K) * UNVBC
	            Endif
	         Enddo
	      Enddo
	   Endif
C-----   velocity boundary condition
	   If(NUMDBC .NE. 0) Then
	      Do N = 1, NUMDBC
	         ID = IDBC(N)
	         Do K = 1, KBM
 	            FRESH = QDIFF(N) * VDDIST(N,K) / 100.
	            CF(ID,K) = CF(ID,K) + DTI * CDBC(N) / DC(ID) *
     &		          FRESH * (1. + Sign(1.,FRESH)) / 2. +
     &                    DTI * AREA(ID) * DZ(K) * C(ID,K) *
     &		          FRESH * (1. - Sign(1.,FRESH)) / 2.

	         Enddo
	      Enddo
	   Endif
!$OMP END MASTER
!$OMP BARRIER
C============================================================================C
C                     Step forward in time                                   C
C============================================================================C
	Do K = 1, KBM
!$OMP DO	
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         CF(I,K) = C(I,K) * AREA(I) * DZ(K) + DTI * CF(I,K)
	      Endif
         Enddo
!$OMP END DO         
	Enddo
!$OMP END PARALLEL
c======================  end subroutine program  ============================C
      Return
      End
