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
      Subroutine SMOOTHINGVER(BB)
      Include './Include/OCERM_INF'
      Dimension BB(KB+1)
      Dimension TEMP(5),BB_TEMP(0:KB+1)
C-----  SMOOTHING BY MIDPOINT METHOD  --------------------------------------C
      Do K = 1, KB + 1
         BB_TEMP(K) = BB(K)
      Enddo
      BB_TEMP(0) = BB(1)
      Do K = 2, KBM
c---    recording the neighbor data  ---------------------------------------c
         TEMP(1) = BB_TEMP(K-2)
         TEMP(2) = BB_TEMP(K-1)
         TEMP(3) = BB_TEMP(K)
         TEMP(4) = BB_TEMP(K+1)
         TEMP(5) = BB_TEMP(K+2)
C---    making a line  c----------------------------------------------------c               
         Do K1 = 1, 4
            Do K2 = 1, 5 - K1
               If(TEMP(K2) .GE. TEMP(K2+1)) Then
                  TEMPDATA = TEMP(K2)
                  TEMP(K2) = TEMP(K2+1)
                  TEMP(K2+1) = TEMPDATA
               Endif
            Enddo
         Enddo 
C---   picking up the middle value   ---------------------------------------c                
         BB(K) = TEMP(3)  
      Enddo
      Return
      End