C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (02/07/2015)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c      smoothing color function by Midpoint method                          #
c                                                                           #
c############################################################################
      Subroutine SMOOTHINGNOR(BB)
      Include './Include/OCERM_INF'
      Dimension BB(IJM,KB)
      Dimension TEMP(10)
C-----  SMOOTHING BY MIDPOINT METHOD  --------------------------------------C
C!$OMP PARALLEL DEFAULT(SHARED) 
C!$OMP&         PRIVATE(I,J,K,K1,K2,I_COUNT,TEMP)
C!$OMP DO
      Do I = 1, IJM 
          If(CCM(I) .EQ. 1.0) Then
             Do K = 2, KBM
                I_COUNT = 0
                Do J = 1, CELL_POLYGEN(I)                  ! Neibour points
	             If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	                I_COUNT = I_COUNT + 1
	                TEMP(I_COUNT) = BB(CELL_SIDE(I,J,2),K)
	             Else
	                I_COUNT = I_COUNT + 1
	                TEMP(I_COUNT) = BB(I,K)
	             Endif   
	          Enddo
	          I_COUNT = I_COUNT + 1                      ! Up and down points
	          TEMP(I_COUNT) = BB(I,K-1)
	          I_COUNT = I_COUNT + 1
	          TEMP(I_COUNT) = BB(I,K)
                If(K .EQ. KBM) Then
                   I_COUNT = I_COUNT + 1
                   TEMP(I_COUNT) = BB(I,K)
                Else
                   I_COUNT = I_COUNT + 1
                   TEMP(I_COUNT) = BB(I,K+1)
                Endif	                
c-------        arranging the neighbor points  in a ascending line   ---------c
                Do K1 = 1, I_COUNT - 1
                   Do K2 = 1, I_COUNT - K1
                      If(TEMP(K2) .GE. TEMP(K2+1)) Then
                         TEMPDATA = TEMP(K2)
                         TEMP(K2) = TEMP(K2+1)
                         TEMP(K2+1) = TEMPDATA
                      Endif
                   Enddo
                Enddo 
C-------        pickup the middle value as the smoothed value  ---------------c
                If(Mod(I_COUNT,2) .EQ. 0.0) Then
                   BB(I,K) = (TEMP(I_COUNT/2) + TEMP(I_COUNT/2+1)) / 2.
                Else
                   BB(I,K) = TEMP((I_COUNT-1) / 2 +1)
                Endif
             Enddo
          Endif   
      Enddo   
C!$OMP END DO     
C-----  UPDATING THE SOURCE TERMS  -----------------------------------------C
C!$OMP DO
C      Do I = 1, IJM
C         If(CCM(I) .EQ. 1.0) Then
C            Do K = 2, KBM
c              BB(I,K) = S_NEW(I,K)
C            Enddo    
C         Endif
C      Enddo
C!$OMP END DO
C!$OMP END PARALLEL     
      Return
      End