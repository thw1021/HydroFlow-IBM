c-----------------------------------------------------------------------------------c
c                   SUBROUTINE PROGRAM                                              #
C                 VERSION 1.0 (13/04/2015)                                          #
C                 AUTHORIZED BY CHE HAIOU                                           #
C                               SHANGHAI JIAO TONG UNIVERSITY                       #
C                               SHANGHAI, CHINA                                     #
c-----------------------------------------------------------------------------------#
c                      DEM   HardModule                                             #
c-----------------------------------------------------------------------------------c
      Subroutine HARDMODULE
      Include './Include/OCERM_INF'
      Double Precision UNDEM,DIJDEM,PNXDEM,PNYDEM,PNZDEM,PNLDEM
      Dimension DXDEM(I_PARTICLE),DYDEM(I_PARTICLE),DZDEM(I_PARTICLE)
      Dimension DUDEM(I_PARTICLE),DVDEM(I_PARTICLE),DWDEM(I_PARTICLE)
      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,UNDEM,DIJDEM,
!$OMP&         PNXDEM,PNYDEM,PNZDEM,PNLDEM,DXDEM,DYDEM,
!$OMP&         DZDEM,DUDEM,DVDEM,DWDEM)
!$OMP DO      
      Do I = 1 , I_PARTICLE
        DXDEM(I) = 0.
        DYDEM(I) = 0.
        DZDEM(I) = 0.
        DUDEM(I) = 0.
        DVDEM(I) = 0.
        DWDEM(I) = 0.
		Do J = 1 , I_PARTICLE
		    If ( I .NE. J ) Then
c-----distance between two particles		    
                DIJDEM = Sqrt( (XPDEM(I) - XPDEM(J) )**2 + ( YPDEM(I) - 
     &                   YPDEM(J) )**2 + ( ZPDEM(I) - ZPDEM(J) )**2 )		
			    If ( DIJDEM .LT. ( RPDEM(I) + RPDEM(J) ) ) Then
c------unit vector between particles' center			    
			        If ( Abs(DIJDEM - 0.) .LE. acdem ) Then
                        PNXDEM = XPDEM(I) - XQDEM(J)
                        PNYDEM = YPDEM(I) - YQDEM(J)
                        PNZDEM = ZPDEM(I) - ZQDEM(J)
                    Else
                        PNXDEM = XPDEM(I) - XPDEM(J)
                        PNYDEM = YPDEM(I) - YPDEM(J)
                        PNZDEM = ZPDEM(I) - ZPDEM(J) 
                    Endif
                    PNLDEM = Sqrt( PNXDEM**2 + PNYDEM**2 + PNZDEM**2 )
                    If ( Abs(PNLDEM - 0.) .LE. acdem ) Then
                        PNXDEM = XQDEM(I) - XPDEM(J)
                        PNYDEM = YQDEM(I) - YPDEM(J)
                        PNZDEM = ZQDEM(I) - ZPDEM(J)
                        PNLDEM = Sqrt( PNXDEM**2 + PNYDEM**2 +PNZDEM**2)
                    Endif
                    If(Abs(PNLDEM - 0.) .LE. acdem) Then
                        Print*, 'Error of PNLDEM in hardmodule.f!'
                        Pause
                    Endif
                    PNXDEM = PNXDEM / PNLDEM 
                    PNYDEM = PNYDEM / PNLDEM 
                    PNZDEM = PNZDEM / PNLDEM
c------ intermediate variable                 
                    UNDEM =( UPDEM(I) - UPDEM(J) ) * PNXDEM + (VPDEM(I)-
     &                      VPDEM(J))*PNYDEM+(WPDEM(I)-WPDEM(J)) *PNZDEM
c-------changes of particle' velocity			
			        DUDEM(I) = DUDEM(I) - PMDEM(J) / ( PMDEM(I) + PMDEM(J) )
     &	                    * (1. + ePPDEM) * UNDEM * PNXDEM  
				
			        DVDEM(I) = DVDEM(I) - PMDEM(J) / ( PMDEM(I) + PMDEM(J) )
     &	                    * (1. + ePPDEM) * UNDEM * PNYDEM 

			        DWDEM(I) = DWDEM(I) - PMDEM(J) / ( PMDEM(I) + PMDEM(J) ) 
     &	                    * (1. + ePPDEM) * UNDEM * PNZDEM 
     
C----- Impenetrability  Of PARTICLES: distribute the overlap to two particles based on  density
			        DXDEM(I) = DXDEM(I) + (RPDEM(I) + RPDEM(J) - DIJDEM) * 
     &	                    (RPMEAN(J) / (RPMEAN(I) + RPMEAN(J))) *PNXDEM
			        DYDEM(I) = DYDEM(I) + (RPDEM(I) + RPDEM(J) - DIJDEM) * 
     &	                    (RPMEAN(J) / (RPMEAN(I) + RPMEAN(J))) *PNYDEM
			        DZDEM(I) = DZDEM(I) + (RPDEM(I) + RPDEM(J) - DIJDEM) * 
     &	                    (RPMEAN(J) / (RPMEAN(I) + RPMEAN(J))) *PNZDEM
C---------                
			    Endif
		    Endif	    
		Enddo
	Enddo
!$OMP END DO
C-----Calculate the change of velocity and angular velocity of particles
!$OMP DO      
      Do I = 1 , I_PARTICLE
        XPDEM(I) = XPDEM(I) + DXDEM(I)
        YPDEM(I) = YPDEM(I) + DYDEM(I)
        ZPDEM(I) = ZPDEM(I) + DZDEM(I)
        UPDEM(I) = UPDEM(I) + DUDEM(I)
        VPDEM(I) = VPDEM(I) + DVDEM(I)
        WPDEM(I) = WPDEM(I) + DWDEM(I)        
	Enddo
!$OMP END DO
!$OMP END PARALLEL	
      
      Return
      End