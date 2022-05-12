C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (10/07/2009)                                  #
C                 AUTHORIZED BY Jian Wang                                   #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                 Boundary condition iteration module                       #
C                                                                           #
C############################################################################
      Subroutine BCITERATION
      
      Include './Include/OCERM_INF'
      
      Real(kind=8) total_v,Q_grid   
      Integer num
      Logical Path_alive
       
      Inquire(directory="./Include/Iteration_BC/",exist=Path_alive)
      If (.not. Path_alive)                                            
     &    Call system(trim(adjustl("mkdir .\Include\Iteration_BC"))) 
     
      Open(2001,file='./Include/Iteration_BC/infl.QBC')            
      Open(2002,file='./Include/Iteration_BC/VIS.QBC')             
      Open(2003,file='./Include/Iteration_BC/Temp_QKW.dat')  
      Open(2004,file='./Results/ZDES_SA.VEL',status='UNKNOWN')        
      
      total_v = 0.0                                                  
      num = 160                                                       
      Q_grid = 0.95*0.04*0.1*0.04    !percentage*U*H*Gridsize           
                                                    
      Do K = 1,KBM 
          total_v = total_v + DC(num) * DZ(K) * UR(num,K)  
      Enddo                                                            
      Do I = 1,2                                                      
          Write (2001,1000) I,1,(100.0*DC(num) * DZ(K) * UR(num,K) 
     &                                           /total_v,K=1,KBM) 
      Enddo 
      
      Write (2001,1100) 0.0                                          
      Write (2001,1100) Q_grid,Q_grid                                 
      Write (2001,1100) 9999.0                                         
      Write (2001,1100) Q_grid,Q_grid                                  

      Write (2002,1100) 0.0                                           
      Do I = 1,2                                                        
          Write (2002,1200) (TKE(num,K), K=1,KBM)                     
          Write (2002,1200) (TDISS(num,K), K=1,KBM)                    
      Enddo                                                            
      Write (2002,1100) 9999.0                                        
      Do I = 1,2                                                     
          Write (2002,1200) (TKE(num,K), K=1,KBM)                     
          Write (2002,1200) (TDISS(num,K), K=1,KBM)                   
      Enddo                                                          

      Write(2003,*) KBM                                                 
      Do K = 1,KBM                                                    
      Write(2003,'(I10,F14.6,E14.6,E14.6)')K,2001.0*DC(num)*DZ(K)*UR  
     &      (num,K)/total_v,TKE(num,K),TDISS(num,K)               
      Enddo                                                            
      
      num = 600                        
      Do K = KBM, 1, -1                                                 
          Write(2004,'(100F12.6)')                                     
     &            DC(num)*(1+ZZ(K))*sqrt(tbx(num)/1000.)/0.000001,  
     &            UR(num,K)/sqrt(tbx(num)/1000.)               
      Enddo   
      
      Close (2001)                                                     
      Close (2002)                                                  
      Close (2003)
      Close (2004)
                                              
1000  Format(2I14,1000F14.6)                                           
1100  Format(8F14.6)                                         
1200  Format(8E14.6)                                             
      
      End Subroutine
