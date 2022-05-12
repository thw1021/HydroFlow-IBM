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
      Dimension S_NEW(IJM,KBM),TEMP(20)
C-----  SMOOTHING BY MIDPOINT METHOD  --------------------------------------C
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP&         PRIVATE(I,J,K,K1,K2,I_COUNT,TEMP,TEMPDATA,FLAG)
!$OMP DO
      Do I = 1, IJM 
          If(CCM(I) .EQ. 1.0) Then
             S_NEW(I,1) = 0.0
             S_NEW(I,2) = (2. * BB(I,1) + 27. * BB(I,2) + 12. * BB(I,3)
     &                     - 8. * BB(I,4) + 2. * BB(I,5)) / 35.
             Do K = 3, KBM - 2
                S_NEW(I,K) =(-3.*BB(I,K-2) + 12.*BB(I,K-1) + 17.*BB(I,K)
     &                       + 12.*BB(I,K+1) - 3.*BB(I,K+2)) / 35.
             Enddo
             S_NEW(I,KBM-1) = (2.*BB(I,KBM-4) - 8.*BB(I,KBM-3) + 
     &           12.*BB(I,KBM-2) + 27.*BB(I,KBM-1) + 2.*BB(I,KBM)) / 35.
             S_NEW(I,KBM) = (-BB(I,KBM-4) + 4.*BB(I,KBM-3) - 
     &           6.*BB(I,KBM-2) + 4.*BB(I,KBM-1) + 69.*BB(I,KBM)) / 70.
          Endif   
      Enddo   
!$OMP END DO        
C-----  UPDATING THE SOURCE TERMS  -----------------------------------------C
!$OMP DO
      Do I = 1, IJM
         If(CCM(I) .EQ. 1.0) Then
            Do K = 2, KBM
              BB(I,K) = S_NEW(I,K)
            Enddo    
         Endif
      Enddo
!$OMP END DO
!$OMP END PARALLEL     
      Return
      End