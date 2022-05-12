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
      Subroutine SMOOTHING(BB)
      Include './Include/OCERM_INF'
      Dimension BB(IJM,KB)
      Dimension S_NEW(IJM,KBM),TEMP(5),BB_TEMP(0:KBM+2)
C-----  SMOOTHING BY MIDPOINT METHOD  --------------------------------------C
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP&         PRIVATE(I,J,K,K1,K2,I_COUNT,TEMP,TEMPDATA,FLAG,BB_TEMP)
!$OMP DO
      Do I = 1, IJM 
          If(CCM(I) .EQ. 1.0) Then
             Do K = 1, KBM
                BB_TEMP(K) = BB(I,K)
             Enddo
             BB_TEMP(0) = BB(I,1)
             BB_TEMP(KBM+1) = BB(I,KBM)
             BB_TEMP(KBM+2) = BB(I,KBM)
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
                BB(I,K) = TEMP(3)  
             Enddo
          Endif   
      Enddo   
!$OMP END DO     
C-----  UPDATING THE SOURCE TERMS  -----------------------------------------C
!$OMP DO
      Do I = 1, IJM
         If(CCM(I) .EQ. 1.0) Then
            Do K = 2, KBM
c              BB(I,K) = S_NEW(I,K)
            Enddo    
         Endif
      Enddo
!$OMP END DO
!$OMP END PARALLEL     
      Return
      End