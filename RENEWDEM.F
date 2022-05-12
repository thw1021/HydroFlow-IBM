c-----------------------------------------------------------------------------------c
c                   SUBROUTINE PROGRAM                                              #
C                 VERSION 1.0 (13/04/2015)                                          #
C                 AUTHORIZED BY CHE HAIOU                                           #
C                               SHANGHAI JIAO TONG UNIVERSITY                       #
C                               SHANGHAI, CHINA                                     #
c-----------------------------------------------------------------------------------#
c                       DEM  Renew the new coordinates and IP                       #
c-----------------------------------------------------------------------------------c
      Subroutine RENEWDEM
      Include './Include/OCERM_INF'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)
!$OMP DO      
      Do I = 1 , I_PARTICLE
        If ((IP(I) .LT. 1) .OR. (IP(I) .GT. IJM)) Then
            Print*, 'Error of IP in RENEWDEM!'
            Print*, 'IP=',IP(I)
            Pause
        Endif 
        If ((KP(I) .LT. 1) .OR. (KP(I) .GT. KBM+1)) Then
            Print*, 'Error of KP in RENEWDEM!'
            Print*, 'KP=',KP(I)
            Pause
        Endif 
           
        IPP(I) = IP(I)  
        KPP(I) = KP(I)
       
        XQDEM(I) = XPDEM(I)
        YQDEM(I) = YPDEM(I)
        ZQDEM(I) = ZPDEM(I)
      Enddo
!$OMP END DO
!$OMP END PARALLEL      
      Return
      End