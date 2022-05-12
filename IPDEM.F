c-----------------------------------------------------------------------------------c
c                   SUBROUTINE PROGRAM                                              #
C                 VERSION 1.0 (13/04/2015)                                          #
C                 AUTHORIZED BY CHE HAIOU                                           #
C                               SHANGHAI JIAO TONG UNIVERSITY                       #
C                               SHANGHAI, CHINA                                     #
c-----------------------------------------------------------------------------------#
c                       DEM position  of  particles  in  fluid grid                 #
c-----------------------------------------------------------------------------------c
      Subroutine IPDEM
      Include './Include/OCERM_INF'
      Double Precision X1DEM,X2DEM,Y1DEM,Y2DEM,DIJDEM,XPGDEM

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,
!$OMP&         X1DEM,X2DEM,Y1DEM,Y2DEM,DIJDEM,XPGDEM)
!$OMP DO     
      Do I = 1 , I_PARTICLE 
        IP(I) = 0         
        Do J = 1 , IJM
c-----distance between particle and center of a grid        
            DIJDEM = Sqrt( ( XPDEM(I) - CXY(J,1) )**2 + 
     &                     ( YPDEM(I) - CXY(J,2) )**2 )
            If( DIJDEM .LE. ( Ktimes * Sqrt(AREA(J)) ) ) Then
c----The grid node coordinates 
                Do K = 1 , Cell_POLYGEN(J)       
                    X1DEM = PXY(CELL_SIDE(J,K,3),1)
                    Y1DEM = PXY(CELL_SIDE(J,K,3),2)
                    X2DEM = PXY(CELL_SIDE(J,K,4),1)
                    Y2DEM = PXY(CELL_SIDE(J,K,4),2)           
c----calculate and judge whether a particle belongs the grid       
                    XPGDEM = ( XPDEM(I) - X1DEM ) * ( Y2DEM - Y1DEM ) - 
     &                       ( YPDEM(I) - Y1DEM ) * ( X2DEM - X1DEM )
                    If(XPGDEM .GT. 0.) Then
                        Exit
                    Endif    
                Enddo
                
                If(K .EQ. ( Cell_POLYGEN(J) + 1 ) ) Then             
                    IP(I) = J
                    Exit
                Endif
            Endif 
        Enddo
c------if a particle doesn't belong any grid, let the IP equals IJM+1         
        If (IP(I) .EQ. 0 ) Then
c            print*, 'A particle can not find his grid!'
c            Print*, I,IP(I),XPDEM(I),YPDEM(I)
            IP(I) = IJM + 1                      
        Endif
c------This is only a test, which is not necessary.        
        If((IP(I) .LT. 1) .OR. (IP(I) .GT. (IJM +1)) ) Then
            Print*, 'Error of IP in IPDEM.F!'
            Print*, I,IP(I),XPDEM(I),YPDEM(I)           
            Pause
        Endif 
      Enddo  
!$OMP END DO
!$OMP END PARALLEL       
C-----                
      Return 
      End