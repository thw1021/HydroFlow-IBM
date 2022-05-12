C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (25/05/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c              TVD scheme for surface variants                              #
c                                                                           #
c############################################################################
      Subroutine TVDSCHEMEH(QSUR,Q,GRADX,GRADY,LIMTER)
      Include './Include/OCERM_INF'
	Integer LIMTER
C	Parameter(LIMTER = 8)
	Dimension QSUR(IJE,KB),Q(IJM,KB),GRADX(IJM,KB),GRADY(IJM,KB)
C	Dimension TVDCOE(IJM,KBM,2)
c----------------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,RF,I_UP,I_EDGE) 
	If(LIMTER .GT. 0) Then
	  Do K = 1, KBM
!$OMP DO
	     Do I = 1, IJE

	      If(CFM(I) .EQ. 1.0) Then

	        RF = 2.0 * GRADX(INDEX_EDGE(I,K,1),K) * 
     &		    (CXY(INDEX_EDGE(I,K,2),1) - CXY(INDEX_EDGE(I,K,1),1))+
     &             2.0 * GRADY(INDEX_EDGE(I,K,1),K) * 
     &		     (CXY(INDEX_EDGE(I,K,2),2) - CXY(INDEX_EDGE(I,K,1),2))

		      RF = RF / (Q(INDEX_EDGE(I,K,2),K) -
     &			   Q(INDEX_EDGE(I,K,1),K) + Sign(1.E-15,Q(INDEX_EDGE(I,K,2),K)-
     &			   Q(INDEX_EDGE(I,K,1),K))) - 1.0
       
		      QSUR(I,K) = Q(INDEX_EDGE(I,K,1),K) + 
     &			   0.5 * FUNLIMTER(LIMTER,RF) * 
     &             (Q(INDEX_EDGE(I,K,2),K) - Q(INDEX_EDGE(I,K,1),K))

	      Endif
	    Enddo
!$OMP END DO
	  Enddo
	Endif
	If(LIMTER .EQ. 0) Then                                ! Higher WENO order
		  Do K = 1, KBM
!$OMP DO
	      Do I = 1, IJE
	        If(CFM(I) .EQ. 1.0) Then
	          I_UP = INDEX_EDGE(I,K,1)
	          I_EDGE = 0
	          Do J = 1, CELL_POLYGEN(I_UP)
	             If(CELL_SIDE(I_UP,J,1) .EQ. I) I_EDGE = J 
	          Enddo
	        
                QSUR(I,K) = 0.0
                Do J = 1, NUM_STENCIL(I_UP)
                   If(I_STENCIL(I_UP,J) .GT. 0) Then
                      QSUR(I,K) = QSUR(I,K) + A_ENO(I_UP,I_EDGE,J) *
     &                            Q(I_STENCIL(I_UP,J),K) 
                   Else
                      QSUR(I,K) = QSUR(I,K) + A_ENO(I_UP,I_EDGE,J) *
     &                            QGHOST(I_STENCIL(I_UP,J),K)
                   Endif
                Enddo
                QSUR(I,K) = QSUR(I,K) + A_ENO(I_UP,I_EDGE,0) * Q(I_UP,K)
                
	        Endif
	      Enddo
!$OMP END DO
	   Enddo
	Endif  
!$OMP END PARALLEL
	Return
	End

	Function FUNLIMTER(LIMTER,RF)
	Double precision FUNLIMTER, RF
	Integer LIMTER
	  Goto (1,2,3,4,5,6,7,8,9,10,11) LIMTER
1	Continue                    ! SUPERBEE
	  FUNLIMTER = Dmax1(0.0,Dmin1(1.0,2.*RF), Dmin1(2.0,RF)) 
	Goto 100
2	Continue                    ! Van Leer
	  FUNLIMTER = (RF + DAbs(RF)) / (1.0 + RF + 1.E-10)
	Goto 100
3	Continue                    ! Van Albada
	  FUNLIMTER = (RF + RF ** 2.) / (1. + RF ** 2.)
	Goto 100
4	Continue                    ! Min-Mod
	  If(RF .GT. 0.0) Then
	    FUNLIMTER = Dmin1(RF,1.0)
	  Else
	    FUNLIMTER = 0.0
	  Endif
	Goto 100
5	Continue                    ! Sweby
        FUNLIMTER = Dmax1(0.0, Dmin1(1.0,1.5*RF), Dmin1(1.5,RF))
	Goto 100
6     Continue                    ! QUICK
	  FUNLIMTER = Dmax1(0.0, Dmin1(2.0*RF,(3.+RF)/4.,2.))   
	Goto 100
7	Continue                    ! UMIST
	  FUNLIMTER = Dmax1(0.0, Dmin1(2.0*RF,(1.+3.*RF)/4.,
     &	                (3.+RF)/4.,2.))
	Goto 100
8	Continue                    ! OSHER
	  FUNLIMTER = Dmax1(0.0, Dmin1(2.,RF))
	Goto 100
9     Continue                    ! MUSCL
	  FUNLIMTER = (RF + DAbs(RF)) / (1. + DAbs(RF)) 
	Goto 100
10	Continue                    ! 1ST UNWIND
        FUNLIMTER = 0.0
	Goto 100
11	Continue                    ! MC
	  FUNLIMTER = Dmax1(0.0, Dmin1((1.+RF)/2.,2.,2.*RF))
	Goto 100
100	Continue
	Return
	End	 
