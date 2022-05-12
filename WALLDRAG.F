C###########################################################################C
c                 SUBROUNTINE PROGRAM                                       #
C                 SOLID BOUNDARY CONDITION BASED ON WALL FUNCTION           #
c                               VERSION 1.0 (16/02/2011)                    #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C############################################################################
	Subroutine WALLDRAG
	Include './Include/OCERM_INF'
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(VELMEAN,NUMAVE)
	Do I = 1, IJM
	   VELMEAN = 0.0
	   NUMAVE = 0
	   If(DC(I) .GT. 0.0) Then
C-----   averaging the velocities on the nearest grids to the wall over the 
c-----   around cells instead of all cells on the grid level  	      
	      Do J = 1, CELL_POLYGEN(I)
		     If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
				VELMEAN = VELMEAN + 
     &				      Sqrt(UR(CELL_SIDE(I,J,2),KBM) ** 2. + 
     &                           VR(CELL_SIDE(I,J,2),KBM) ** 2.)
	            NUMAVE = NUMAVE + 1
		     Endif
	      Enddo
	      VELMEAN = VELMEAN + Sqrt(UR(I,KBM) ** 2. + VR(I,KBM) ** 2.)
	      VELMEAN = VELMEAN / (NUMAVE + 1)
	      CBC(I)=0.16 * VELMEAN /
     &              Log((ZZ(KBM)-Z(KB))*DC(I)/Z01(I))**2
	   Endif
	Enddo
!$OMP END PARALLEL DO
	Return
	End