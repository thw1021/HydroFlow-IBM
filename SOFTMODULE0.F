c-----------------------------------------------------------------------------------c
c                             SUBROUTINE PROGRAM                                    #
C                             DESINGED BY ChenJun                                   #
C                                                 SHANGHAI JIAO TONG UNIVERSITY     #
C                                                 SHANGHAI, CHINA                   #
c-----------------------------------------------------------------------------------#
c 子程序功能：
c 1------颗粒接触力和力矩计算，基于遍历碰撞接触搜索算法     
c-----------------------------------------------------------------------------------#
      
      Subroutine SOFTMODULE0
      Include './Include/OCERM_INF'
c      Parameter (ISCHEME = 2, N_SCHEME1 = 0)      
      Double Precision DIJDEM,PNXDEM,PNYDEM,PNZDEM,PNLDEM
      Double Precision UP_IJ,VP_IJ,WP_IJ,UNDEM,PKNDEM,CNDEM,ALFADEM
      Double Precision FNXDEM_IJ, FNYDEM_IJ, FNZDEM_IJ, FNDEM_IJ
      Double Precision VCTXDEM,VCTYDEM,VCTZDEM,VCTNLDEM
      Double Precision DIATAX,DIATAY,DIATAZ,PKTDEM,CTDEM
      Double Precision FTXDEM_IJ,FTYDEM_IJ,FTZDEM_IJ,FTDEM_IJ
      Double Precision TPXDEM_IJ,TPYDEM_IJ,TPZDEM_IJ,OMDEM        
      ALFADEM_MAX = 0.0      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,DIJDEM,PNXDEM,PNYDEM,PNZDEM,
!$OMP&         PNLDEM,UP_IJ,VP_IJ,WP_IJ,UNDEM,PKNDEM,CNDEM,ALFADEM,
!$OMP&         FNXDEM_IJ, FNYDEM_IJ, FNZDEM_IJ, FNDEM_IJ,
!$OMP&         VCTXDEM,VCTYDEM,VCTZDEM,VCTNLDEM,DIATAX,DIATAY,DIATAZ,
!$OMP&         PKTDEM,CTDEM,FTXDEM_IJ,FTYDEM_IJ,FTZDEM_IJ,FTDEM_IJ,
!$OMP&         TPXDEM_IJ,TPYDEM_IJ,TPZDEM_IJ,OMDEM,DIATATN,N_SYM1,TEMP,
!$OMP&         DOMX,DOMY,DOMZ,DOM,TRXDEM_IJ,TRYDEM_IJ,TRZDEM_IJ,
!$OMP&         TRX0,TRY0,TRZ0,VOMX,VOMY,VOMZ,VOM)            
!$OMP DO       
      Do I = 1 , I_PARTICLE
       IF(NSYM_CONT(I) .EQ. 1)THEN  
        FNXDEM(I) = 0.
        FNYDEM(I) = 0.
        FNZDEM(I) = 0.
        FTXDEM(I) = 0.
        FTYDEM(I) = 0.
        FTZDEM(I) = 0.
        TPXDEM(I) = 0.
        TPYDEM(I) = 0.
        TPZDEM(I) = 0.
        TPXDEM_IJ = 0.
        TPYDEM_IJ = 0.
        TPZDEM_IJ = 0.
        TRXDEM_IJ = 0.
        TRYDEM_IJ = 0.
        TRZDEM_IJ = 0.        
        N_P1(I) = 0        
         
        V_PARTICLE(I) = ( 4. / 3. ) *  ACOS(-1.0) * RPDEM(I)**3 
        OMDEM =  DMAX1( Sqrt(OMXDEM(I) * OMXDEM(I) + 
     &                       OMYDEM(I) * OMYDEM(I) +
     &                       OMZDEM(I) * OMZDEM(I)), acdem)        
C---------------------------------------------------------------------------------C
C                  contact among particles, chenjun, 20191124                     c
c---------------------------------------------------------------------------------c    
        Do J = 1 , I_PARTICLE
          If( I .NE. J ) Then
C----- normal contact displacement 
            ALFADEM = RPDEM(I) + RPDEM(J) - Sqrt(  
     &                ( XPDEM(I) - XPDEM(J) )**2 + 
     &                ( YPDEM(I) - YPDEM(J) )**2 + 
     &                ( ZPDEM(I) - ZPDEM(J) )**2 ) 
            
c----- collision judgement                
            If ( ALFADEM .GT. 0.0 ) Then 

C             If(ALFADEM .GT. 0.1*2*Min(RPDEM(I),RPDEM(J))) Then
C               Print*, 'Error1 in softmodule.f: contact error!'  
C               print*, 'ALFADEM=',ALFADEM 
C               print*, '(X,Y,Z)=',XPDEM(I),YPDEM(I),ZPDEM(I)
C               !Pause
C             Endif            

             IF(ALFADEM / RPDEM(I) .GT. ALFADEM_MAX)THEN
                 ALFADEM_MAX = ALFADEM / RPDEM(I)
             ENDIF

             DIATAX = 0.0
             DIATAY = 0.0
             DIATAZ = 0.0
             N_SYM1 = 0
             
             DO K = 1, N_P2(I)   
               IF(J .EQ. N_PP2(I,K))THEN
                  DIATAX = DELTAT_PP1X(I,K)
                  DIATAY = DELTAT_PP1Y(I,K)
                  DIATAZ = DELTAT_PP1Z(I,K)
                  N_SYM1 = 1
               ENDIF
             ENDDO 
                  
             N_P1(I) = N_P1(I) + 1
             If(N_P1(I) .GT. I_CPN) Then
               Print*, 'Error2 in softmodule.f: N_P1(I)'   
               Pause
             Endif
             N_PP1(I,N_P1(I)) = J
                          
c----- unit vector from particle i to j                   
             PNXDEM = XPDEM(J) - XPDEM(I)
             PNYDEM = YPDEM(J) - YPDEM(I)
             PNZDEM = ZPDEM(J) - ZPDEM(I) 
                    
             PNLDEM = Sqrt( PNXDEM**2 + PNYDEM**2 + PNZDEM**2 )
                                                            
             PNXDEM = PNXDEM / PNLDEM  
             PNYDEM = PNYDEM / PNLDEM 
             PNZDEM = PNZDEM / PNLDEM 
                    
c------ ralative velocity between particle i and j
             UP_IJ = UPDEM(I) - UPDEM(J)
             VP_IJ = VPDEM(I) - VPDEM(J)
             WP_IJ = WPDEM(I) - WPDEM(J)
                    
c------ UNDEM                    
             UNDEM = UP_IJ * PNXDEM + VP_IJ * PNYDEM + 
     &               WP_IJ * PNZDEM   
                    
C----- kn and CNDEM  
             IF(CONMOD .EQ. 'M1')THEN
                        
               PKNDEM = PKN_M1PP
               CNDEM = 2.0 * SQRT( PMDEM(I) * PMDEM(J) / 
     &		         ( PMDEM(I) + PMDEM(J) ) * PKNDEM ) * 
     &                 ABS( LOG(EN_M1PP) ) / 
     &                 SQRT( PI_CON**2 + ( LOG(EN_M1PP) )**2 )
                    
             ELSEIF(CONMOD .EQ. 'M2')THEN    
                        
               PKNDEM = 4./3. / ( (1. - PNIUPDEM(I)**2)/EPDEM(I) +
     &				  (1. - PNIUPDEM(J)**2) / EPDEM(J) ) /  
     &				  ( ( RPDEM(I) + RPDEM(J) ) / 
     &                    ( RPDEM(I) * RPDEM(J) ) )**0.5
                    
               CNDEM = 2. * Sqrt( PKNDEM * ( PMDEM(I) * PMDEM(J) / 
     &		         ( PMDEM(I) + PMDEM(J) ) ) ) 
                    
             ELSEIF(CONMOD .EQ. 'M3')THEN
                        
               PKNDEM = 2. * PI_CON**2 * PMDEM(I) / (20.0*TDEM)**2
               CNDEM = 2.0 * 0.3 * Sqrt(PMDEM(I)*PKNDEM)
                    
             ENDIF                        
                                       
C----- fn(i,j)
             FNDEM_IJ =  - PKNDEM * ALFADEM - CNDEM * UNDEM 
                    
             FNXDEM_IJ = FNDEM_IJ * PNXDEM				            
             FNYDEM_IJ = FNDEM_IJ * PNYDEM			           
             FNZDEM_IJ = FNDEM_IJ * PNZDEM    
             
             FNDEM_IJ = ABS( FNDEM_IJ )  		
                    
C----- VT(i,j)
             VCTXDEM = UP_IJ - UNDEM * PNXDEM +
     &			     RPDEM(I) * (OMYDEM(I)*PNZDEM - OMZDEM(I)*PNYDEM) +
     &			     RPDEM(J) * (OMYDEM(J)*PNZDEM - OMZDEM(J)*PNYDEM) 
     	
		   VCTYDEM = VP_IJ - UNDEM * PNYDEM -		
     &			     RPDEM(I) * (OMXDEM(I)*PNZDEM - OMZDEM(I)*PNXDEM) -
     &			     RPDEM(J) * (OMXDEM(J)*PNZDEM - OMZDEM(J)*PNXDEM) 
     			
		   VCTZDEM = WP_IJ - UNDEM * PNZDEM +
     &			     RPDEM(I) * (OMXDEM(I)*PNYDEM - OMYDEM(I)*PNXDEM) +
     &			     RPDEM(J) * (OMXDEM(J)*PNYDEM - OMYDEM(J)*PNXDEM)
     
             VCTNLDEM = Sqrt(VCTXDEM**2+VCTYDEM**2+VCTZDEM**2) 
                    
C----- DIATAX
             IF(N_SYM1 .EQ. 0)THEN
               DIATAX = VCTXDEM * DMIN1(
     &                  ALFADEM / Dmax1(UNDEM,acdem), TDEM )
               DIATAY = VCTYDEM * DMIN1(
     &                  ALFADEM / Dmax1(UNDEM,acdem), TDEM )
               DIATAZ = VCTZDEM * DMIN1(
     &                  ALFADEM / Dmax1(UNDEM,acdem), TDEM )
             ELSE
               DIATATN = DIATAX * PNXDEM + DIATAY * PNYDEM 
     &                                   + DIATAZ * PNZDEM
                        
               DIATAX = DIATAX - DIATATN * PNXDEM
     &                         + VCTXDEM * TDEM 
               DIATAY = DIATAY - DIATATN * PNYDEM
     &                         + VCTYDEM * TDEM
               DIATAZ = DIATAZ - DIATATN * PNZDEM
     &                         + VCTZDEM * TDEM   
             ENDIF                    
                    
C----- kt and CTDEM
             IF(CONMOD .EQ. 'M1')THEN
                        
               PKTDEM = RK_M1PP * PKNDEM
               CTDEM = RC_M1PP * CNDEM           
                    
             ELSEIF(CONMOD .EQ. 'M2')THEN
                        
      		 PKTDEM = 8. * ALFADEM**0.5 / ( (1. - PNIUPDEM(I)**2) /
     &				  GPDEM(I) + (1. - PNIUPDEM(J)**2) /GPDEM(J) ) /
     &				  ( (RPDEM(I) + RPDEM(J)) / (RPDEM(I) * RPDEM(J)) )**0.5
                    
			 CTDEM = 2. * Sqrt( ( PKTDEM * PMDEM(I) * PMDEM(J) / (PMDEM(I)
     &				 + PMDEM(J) ) ) )  
                    
             ELSEIF(CONMOD .EQ. 'M3')THEN
                        
               PKTDEM = PKNDEM / ( 2.0 * (1.0 + PNIUPDEM(I)) )
               CTDEM = CNDEM / ( 2.0 * (1.0 + PNIUPDEM(I)) )
                    
             ENDIF                    
                      
C----- ft(i,j) 
             FTXDEM_IJ = - PKTDEM * DIATAX - CTDEM * VCTXDEM 
             FTYDEM_IJ = - PKTDEM * DIATAY - CTDEM * VCTYDEM
             FTZDEM_IJ = - PKTDEM * DIATAZ - CTDEM * VCTZDEM
                    
             FTDEM_IJ =  Sqrt( FTXDEM_IJ * FTXDEM_IJ + 
     &                         FTYDEM_IJ * FTYDEM_IJ +
     &                         FTZDEM_IJ * FTZDEM_IJ )           
                    
C----- Whether slipe, if slipe then ft will equal to 'PMIUS * fn	'
		   If ( FTDEM_IJ .GT. ( PMIUS * FNDEM_IJ ) ) Then
                
               FTXDEM_IJ = PMIUS * FNDEM_IJ * FTXDEM_IJ / 
     &                     Dmax1(FTDEM_IJ,acdem)
			 FTYDEM_IJ = PMIUS * FNDEM_IJ * FTYDEM_IJ / 
     &                     Dmax1(FTDEM_IJ,acdem)
			 FTZDEM_IJ = PMIUS * FNDEM_IJ * FTZDEM_IJ / 
     &                     Dmax1(FTDEM_IJ,acdem)                   
                       
               IF(N_SYM1 .EQ. 1)THEN                       
                 DELTAT_PP1X(I,N_P1(I)) = DIATAX - VCTXDEM * TDEM
     &                                         + DIATATN * PNXDEM
                 DELTAT_PP1Y(I,N_P1(I)) = DIATAY - VCTYDEM * TDEM
     &                                         + DIATATN * PNYDEM
                 DELTAT_PP1Z(I,N_P1(I)) = DIATAZ - VCTZDEM * TDEM
     &                                         + DIATATN * PNZDEM
               ENDIF

             ELSE
                 
               DELTAT_PP1X(I,N_P1(I)) = DIATAX
               DELTAT_PP1Y(I,N_P1(I)) = DIATAY
               DELTAT_PP1Z(I,N_P1(I)) = DIATAZ  
               
             Endif                     

             TPXDEM_IJ = TPXDEM_IJ + RPDEM(I) *  
     &                   ( PNYDEM * FTZDEM_IJ - PNZDEM * FTYDEM_IJ )   
                    
             TPYDEM_IJ = TPYDEM_IJ + RPDEM(I) *  
     &                   ( PNZDEM * FTXDEM_IJ - PNXDEM * FTZDEM_IJ )   
                   
		   TPZDEM_IJ = TPZDEM_IJ + RPDEM(I) *  
     &                   ( PNXDEM * FTYDEM_IJ - PNYDEM * FTXDEM_IJ )   

C----- rolling friction torque               
             IF(ROLLMOD .EQ. 'CON1')THEN
                 
               TRXDEM_IJ = TRXDEM_IJ - RMIUS * FNDEM_IJ *
     &                                 OMXDEM(I) / OMDEM
               
               TRYDEM_IJ = TRYDEM_IJ - RMIUS * FNDEM_IJ * 
     &                                 OMYDEM(I) / OMDEM    
               
               TRZDEM_IJ = TRZDEM_IJ - RMIUS * FNDEM_IJ * 
     &                                 OMZDEM(I) / OMDEM 
               
             ELSEIF(ROLLMOD .EQ. 'CON2')THEN
                 
               DOMX = OMXDEM(I) - OMXDEM(J)
               DOMY = OMYDEM(I) - OMYDEM(J)
               DOMZ = OMZDEM(I) - OMZDEM(J)
               
               DOM =  DMAX1( Sqrt( DOMX**2 + DOMY**2 + DOMZ**2 ), acdem)
               
               TRXDEM_IJ = TRXDEM_IJ -  RMIUS * FNDEM_IJ * DOMX / DOM
               TRYDEM_IJ = TRYDEM_IJ -  RMIUS * FNDEM_IJ * DOMY / DOM
               TRZDEM_IJ = TRZDEM_IJ -  RMIUS * FNDEM_IJ * DOMZ / DOM
               
             ELSEIF(ROLLMOD .EQ. 'VIS1')THEN
                 
               VOMX = RPDEM(I) * ( OMYDEM(I)*PNZDEM-OMZDEM(I)*PNYDEM ) +
     &                RPDEM(J) * ( OMYDEM(J)*PNZDEM-OMZDEM(J)*PNYDEM )
               
               VOMY = RPDEM(I) * ( OMZDEM(I)*PNXDEM-OMXDEM(I)*PNZDEM ) +
     &                RPDEM(J) * ( OMZDEM(J)*PNXDEM-OMXDEM(J)*PNZDEM )
               
               VOMZ = RPDEM(I) * ( OMXDEM(I)*PNYDEM-OMYDEM(I)*PNXDEM ) +
     &                RPDEM(J) * ( OMXDEM(J)*PNYDEM-OMYDEM(J)*PNXDEM )
               
               VOM = SQRT( VOMX**2 + VOMY**2 + VOMZ**2 )
               
               TRXDEM_IJ = TRXDEM_IJ - RMIUS * VOM *
     &                     FNDEM_IJ * OMXDEM(I) / OMDEM
               
               TRYDEM_IJ = TRYDEM_IJ - RMIUS * VOM *
     &                     FNDEM_IJ * OMYDEM(I) / OMDEM    
               
               TRZDEM_IJ = TRZDEM_IJ - RMIUS * VOM *
     &                     FNDEM_IJ * OMZDEM(I) / OMDEM                 
               
             ELSEIF(ROLLMOD .EQ. 'VIS2')THEN      
                 
               VOMX = RPDEM(I) * ( OMYDEM(I)*PNZDEM-OMZDEM(I)*PNYDEM ) +
     &                RPDEM(J) * ( OMYDEM(J)*PNZDEM-OMZDEM(J)*PNYDEM )
               
               VOMY = RPDEM(I) * ( OMZDEM(I)*PNXDEM-OMXDEM(I)*PNZDEM ) +
     &                RPDEM(J) * ( OMZDEM(J)*PNXDEM-OMXDEM(J)*PNZDEM )
               
               VOMZ = RPDEM(I) * ( OMXDEM(I)*PNYDEM-OMYDEM(I)*PNXDEM ) +
     &                RPDEM(J) * ( OMXDEM(J)*PNYDEM-OMYDEM(J)*PNXDEM )  

               TRXDEM_IJ = TRXDEM_IJ -  RMIUS * VOMX * FNDEM_IJ 
               TRYDEM_IJ = TRYDEM_IJ -  RMIUS * VOMY * FNDEM_IJ
               TRZDEM_IJ = TRZDEM_IJ -  RMIUS * VOMZ * FNDEM_IJ      
               
             ENDIF             
                    
             FNXDEM(I) = FNXDEM(I) + FNXDEM_IJ                 
             FNYDEM(I) = FNYDEM(I) + FNYDEM_IJ  
             FNZDEM(I) = FNZDEM(I) + FNZDEM_IJ
                    
             FTXDEM(I) = FTXDEM(I) + FTXDEM_IJ
             FTYDEM(I) = FTYDEM(I) + FTYDEM_IJ
             FTZDEM(I) = FTZDEM(I) + FTZDEM_IJ
                                   
            Endif
          Endif 
        Enddo         
C---------------------------------------------------------------------------------C
C             contact between particle and boundary, chenjun, 20191124            c
c---------------------------------------------------------------------------------c                  
        IF(NP_BCD(I) .NE. 0)THEN               

          Do J = 1, NP_BCD(I) 
              
            IF(NCONT_INFO(I,J) .NE. 7)THEN  

              PNXDEM = BCDEM_INFO(I,3*J-2) - XPDEM(I)
              PNYDEM = BCDEM_INFO(I,3*J-1) - YPDEM(I)
              PNZDEM = BCDEM_INFO(I,3*J) - ZPDEM(I)
                    
              PNLDEM = Sqrt( PNXDEM**2 + PNYDEM**2 + PNZDEM**2 )    
                
              ALFADEM = 2 * RPDEM(I)- PNLDEM 
                
C              If(ALFADEM .GT. 0.1*2*RPDEM(I)) Then
C                Print*, 'Error1 in softmodule.f: contact error!' 
C                print*, 'ALFADEM=',ALFADEM 
C                print*, '(X,Y,Z)=',XPDEM(I),YPDEM(I),ZPDEM(I)
C                !Pause
C              Endif    
                
              IF(ALFADEM / RPDEM(I) .GT. ALFADEM_MAX)THEN
                 ALFADEM_MAX = ALFADEM / RPDEM(I)
              ENDIF                

              DIATAX = 0.0
              DIATAY = 0.0
              DIATAZ = 0.0
              N_SYM1 = 0
              
              DO K = 1, N_P2(I)   
                IF(NCONT_INFO(I,J) + I_PARTICLE .EQ. N_PP2(I,K))THEN
                  DIATAX = DELTAT_PP1X(I,K)
                  DIATAY = DELTAT_PP1Y(I,K)
                  DIATAZ = DELTAT_PP1Z(I,K)
                  N_SYM1 = 1
                ENDIF
              ENDDO 
                  
              N_P1(I) = N_P1(I) + 1                 
              If(N_P1(I) .GT. I_CPN) Then
                Print*, 'Error3 in softmodule.f: N_P1(I)'   
                !Pause
              Endif
              N_PP1(I,N_P1(I)) = NCONT_INFO(I,J) + I_PARTICLE
                
c----- unit vector                                                                                                                                            
              PNXDEM = PNXDEM / PNLDEM  
              PNYDEM = PNYDEM / PNLDEM 
              PNZDEM = PNZDEM / PNLDEM

c------ ralative velocity between particle i and j
              UP_IJ = UPDEM(I)
              VP_IJ = VPDEM(I)
              WP_IJ = WPDEM(I)
                !print*, UPDEM(I)   
c------ UNDEM                    
              UNDEM = UP_IJ * PNXDEM + VP_IJ * PNYDEM + 
     &                WP_IJ * PNZDEM   
                    
C----- kn and CNDEM 
              IF(CONMOD .EQ. 'M1')THEN
                  
                PKNDEM = PKN_M1PW
                CNDEM = 2.0 * SQRT( PMDEM(I) * PMDEM(I) / 
     &		          ( PMDEM(I) + PMDEM(I) ) * PKNDEM ) * 
     &                  ABS( LOG(EN_M1PP) ) / 
     &                  SQRT(PI_CON**2+(LOG(EN_M1PW))**2)  
                
              ELSEIF(CONMOD .EQ. 'M2')THEN 
                  
                PKNDEM = 4./3. / ( (1. - PNIUPDEM(I)**2)/EPDEM(I) +
     &				   (1. - PNIUPDEM(I)**2) / EPDEM(I) ) /  
     &				   ( ( RPDEM(I) + RPDEM(I) ) / 
     &                   ( RPDEM(I) * RPDEM(I) ) )**0.5
                    
                CNDEM = 2. * Sqrt( PKNDEM * ( PMDEM(I) * PMDEM(I) / 
     &		          ( PMDEM(I) + PMDEM(I) ) ) )
                
              ELSEIF(CONMOD .EQ. 'M3')THEN    
                  
                PKNDEM = 2.*PI_CON**2*PMDEM(I)/(20.0*TDEM)**2
                CNDEM = 2.0*0.3*Sqrt(PMDEM(I)*PKNDEM)  

              ENDIF
                                   
C----- fn(i,j)
              FNDEM_IJ =  - PKNDEM * ALFADEM - CNDEM * UNDEM

              FNXDEM_IJ = FNDEM_IJ * PNXDEM				            
			FNYDEM_IJ = FNDEM_IJ * PNYDEM			           
			FNZDEM_IJ = FNDEM_IJ * PNZDEM      

              FNDEM_IJ = ABS( FNDEM_IJ )  
                                    
C----- VT(i,j)
              VCTXDEM = UP_IJ - UNDEM * PNXDEM +
     &			      RPDEM(I) * (OMYDEM(I)*PNZDEM - OMZDEM(I)*PNYDEM)
                
              VCTYDEM = VP_IJ - UNDEM * PNYDEM -		
     &			      RPDEM(I) * (OMXDEM(I)*PNZDEM - OMZDEM(I)*PNXDEM) 
     			
			VCTZDEM = WP_IJ - UNDEM * PNZDEM +
     &			      RPDEM(I) * (OMXDEM(I)*PNYDEM - OMYDEM(I)*PNXDEM)
     
              VCTNLDEM = Sqrt(VCTXDEM**2+VCTYDEM**2+VCTZDEM**2) 

C----- DIATAX   
              IF(N_SYM1 .EQ. 0)THEN
                    
                DIATAX = VCTXDEM * DMIN1(
     &                   ALFADEM / Dmax1(UNDEM,acdem), TDEM )
                DIATAY = VCTYDEM * DMIN1(
     &                   ALFADEM / Dmax1(UNDEM,acdem), TDEM )
                DIATAZ = VCTZDEM * DMIN1(
     &                   ALFADEM / Dmax1(UNDEM,acdem), TDEM )
                   
              ELSE
                    
                DIATATN = DIATAX * PNXDEM + DIATAY * PNYDEM 
     &                                    + DIATAZ * PNZDEM  
                DIATAX = DIATAX - DIATATN * PNXDEM
     &                          + VCTXDEM * TDEM 
                DIATAY = DIATAY - DIATATN * PNYDEM
     &                          + VCTYDEM * TDEM
                DIATAZ = DIATAZ - DIATATN * PNZDEM
     &                          + VCTZDEM * TDEM   
                    
              ENDIF                      
                  
C----- kt and CTDEM  
              IF(CONMOD .EQ. 'M1')THEN
                        
                PKTDEM = RK_M1PW * PKNDEM
                CTDEM = RC_M1PW * CNDEM 
                
              ELSEIF(CONMOD .EQ. 'M2')THEN                    

                PKTDEM = 8. * ALFADEM**0.5 / ( (1.-PNIUPDEM(I)**2) /
     &				   GPDEM(I) + (1. - PNIUPDEM(I)**2) /GPDEM(I) ) /
     &				   ( (RPDEM(I) + RPDEM(I)) / (RPDEM(I) * RPDEM(I)) )**0.5

			  CTDEM = 2. * Sqrt( ( PKTDEM * PMDEM(I) * PMDEM(I) / (PMDEM(I)
     &				  + PMDEM(I) ) ) )
                    
              ELSEIF(CONMOD .EQ. 'M3')THEN
                    
                PKTDEM = PKNDEM / (2.0*(1.0+PNIUPDEM(I)))
                CTDEM = CNDEM / (2.0*(1.0+PNIUPDEM(I)))   
                    
              ENDIF
                
C----- ft(i,j) 
      		FTXDEM_IJ = - PKTDEM * DIATAX - CTDEM * VCTXDEM 
      		FTYDEM_IJ = - PKTDEM * DIATAY - CTDEM * VCTYDEM
              FTZDEM_IJ = - PKTDEM * DIATAZ - CTDEM * VCTZDEM
          
              FTDEM_IJ =  Sqrt( FTXDEM_IJ * FTXDEM_IJ + 
     &                          FTYDEM_IJ * FTYDEM_IJ +
     &                          FTZDEM_IJ * FTZDEM_IJ )           
   
C----- Whether slipe, if slipe then ft will equal to 'PMIUSPW * fn	'
              If ( FTDEM_IJ .GT. ( PMIUSPW * FNDEM_IJ ) ) Then
                   
                FTXDEM_IJ = PMIUSPW * FNDEM_IJ * FTXDEM_IJ / 
     &                      Dmax1(FTDEM_IJ,acdem)
			  FTYDEM_IJ = PMIUSPW * FNDEM_IJ * FTYDEM_IJ / 
     &                      Dmax1(FTDEM_IJ,acdem)
			  FTZDEM_IJ = PMIUSPW * FNDEM_IJ * FTZDEM_IJ / 
     &                      Dmax1(FTDEM_IJ,acdem)              
                   
                IF(N_SYM1 .EQ. 1)THEN
                  DELTAT_PP1X(I,N_P1(I)) = DIATAX - VCTXDEM * TDEM
     &                                            + DIATATN * PNXDEM 
                  DELTAT_PP1Y(I,N_P1(I)) = DIATAY - VCTYDEM * TDEM
     &                                            + DIATATN * PNYDEM 
                  DELTAT_PP1Z(I,N_P1(I)) = DIATAZ - VCTZDEM * TDEM  
     &                                            + DIATATN * PNZDEM 
                ENDIF   

              ELSE
                  
                   DELTAT_PP1X(I,N_P1(I)) = DIATAX
                   DELTAT_PP1Y(I,N_P1(I)) = DIATAY
                   DELTAT_PP1Z(I,N_P1(I)) = DIATAZ    
                   
              Endif
                                                  
			TPXDEM_IJ = TPXDEM_IJ + RPDEM(I) *  
     &                    ( PNYDEM * FTZDEM_IJ - PNZDEM * FTYDEM_IJ )  
                    
			TPYDEM_IJ = TPYDEM_IJ + RPDEM(I) *  
     &                    ( PNZDEM * FTXDEM_IJ - PNXDEM * FTZDEM_IJ )   
                    
              TPZDEM_IJ = TPZDEM_IJ + RPDEM(I) *  
     &                    ( PNXDEM * FTYDEM_IJ - PNYDEM * FTXDEM_IJ ) 

C----- rolling friction torque             
              IF(ROLLMOD .EQ. 'CON1')THEN  
                 
                TRXDEM_IJ = TRXDEM_IJ - RMIUSPW * FNDEM_IJ *
     &                                  OMXDEM(I) / OMDEM
               
                TRYDEM_IJ = TRYDEM_IJ - RMIUSPW * FNDEM_IJ * 
     &                                  OMYDEM(I) / OMDEM    
               
                TRZDEM_IJ = TRZDEM_IJ - RMIUSPW * FNDEM_IJ * 
     &                                  OMZDEM(I) / OMDEM 
               
              ELSEIF(ROLLMOD .EQ. 'CON2')THEN
                                
                TRXDEM_IJ = TRXDEM_IJ -  RMIUSPW * FNDEM_IJ * 
     &                                   OMXDEM(I) / OMDEM
                TRYDEM_IJ = TRYDEM_IJ -  RMIUSPW * FNDEM_IJ * 
     &                                   OMYDEM(I) / OMDEM
                TRZDEM_IJ = TRZDEM_IJ -  RMIUSPW * FNDEM_IJ * 
     &                                   OMZDEM(I) / OMDEM
               
              ELSEIF(ROLLMOD .EQ. 'VIS1')THEN
                 
                VOMX = RPDEM(I) * ( OMYDEM(I) * PNZDEM - 
     &                              OMZDEM(I) * PNYDEM )
               
                VOMY = RPDEM(I) * ( OMZDEM(I) * PNXDEM - 
     &                              OMXDEM(I) * PNZDEM )
               
                VOMZ = RPDEM(I) * ( OMXDEM(I) * PNYDEM - 
     &                              OMYDEM(I) * PNXDEM )
                  
                VOM = SQRT( VOMX**2 + VOMY**2 + VOMZ**2 )

                TRXDEM_IJ = TRXDEM_IJ - RMIUSPW * VOM *
     &                      FNDEM_IJ * OMXDEM(I) / OMDEM
               
                TRYDEM_IJ = TRYDEM_IJ - RMIUSPW * VOM *
     &                      FNDEM_IJ * OMYDEM(I) / OMDEM    
               
                TRZDEM_IJ = TRZDEM_IJ - RMIUSPW * VOM *
     &                      FNDEM_IJ * OMZDEM(I) / OMDEM 
                              
              ELSEIF(ROLLMOD .EQ. 'VIS2')THEN      
                 
                VOMX = RPDEM(I) * ( OMYDEM(I) * PNZDEM - 
     &                                OMZDEM(I) * PNYDEM )
               
                VOMY = RPDEM(I) * ( OMZDEM(I) * PNXDEM - 
     &                                OMXDEM(I) * PNZDEM )
               
                VOMZ = RPDEM(I) * ( OMXDEM(I) * PNYDEM - 
     &                                OMYDEM(I) * PNXDEM )

                TRXDEM_IJ = TRXDEM_IJ -  RMIUSPW * VOMX * FNDEM_IJ 
                TRYDEM_IJ = TRYDEM_IJ -  RMIUSPW * VOMY * FNDEM_IJ
                TRZDEM_IJ = TRZDEM_IJ -  RMIUSPW * VOMZ * FNDEM_IJ  
                
              ELSEIF(ROLLMOD .EQ. 'NONE')THEN
                  
                TRXDEM_IJ = 0.0 
                TRYDEM_IJ = 0.0
                TRZDEM_IJ = 0.0    
               
              ENDIF                  
                
              FNXDEM(I) = FNXDEM(I) + FNXDEM_IJ                 
              FNYDEM(I) = FNYDEM(I) + FNYDEM_IJ  
              FNZDEM(I) = FNZDEM(I) + FNZDEM_IJ
          
              FTXDEM(I) = FTXDEM(I) + FTXDEM_IJ
              FTYDEM(I) = FTYDEM(I) + FTYDEM_IJ
              FTZDEM(I) = FTZDEM(I) + FTZDEM_IJ                    
                              
            ENDIF  
          ENDDO
        ENDIF
        
C---------------------------------------------------------------------------------C
C             bounded rolling friction torque, chenjun, 20200326                  c
c---------------------------------------------------------------------------------c     
        IF(ROLLBOUND .EQ. 'NONE')THEN        
             
          TPXDEM(I) = TPXDEM_IJ + TRXDEM_IJ
          TPYDEM(I) = TPYDEM_IJ + TRYDEM_IJ
          TPZDEM(I) = TPZDEM_IJ + TRZDEM_IJ           
     
        ELSEIF(ROLLBOUND .EQ. 'BOD1')THEN 

          TRX0 = - PIPDEM(I) * OMXDEM(I) / TDEM 
          TRY0 = - PIPDEM(I) * OMYDEM(I) / TDEM
          TRZ0 = - PIPDEM(I) * OMZDEM(I) / TDEM
           
          TPXDEM(I) = TPXDEM_IJ + SIGN( DMIN1( ABS(TRXDEM_IJ) ,
     &                ABS(TRX0) ), TRX0 )
          TPYDEM(I) = TPYDEM_IJ + SIGN( DMIN1( ABS(TRYDEM_IJ) ,
     &                ABS(TRY0) ), TRY0 )
          TPZDEM(I) = TPZDEM_IJ + SIGN( DMIN1( ABS(TRZDEM_IJ) ,
     &                ABS(TRZ0) ), TRZ0 ) 
      
        ELSEIF(ROLLBOUND .EQ. 'BOD2')THEN
             
          TRX0 = - PIPDEM(I) * OMXDEM(I) / TDEM - TPXDEM_IJ
          TRY0 = - PIPDEM(I) * OMYDEM(I) / TDEM - TPYDEM_IJ
          TRZ0 = - PIPDEM(I) * OMZDEM(I) / TDEM - TPZDEM_IJ
           
          TPXDEM(I) = TPXDEM_IJ + SIGN( DMIN1( ABS(TRXDEM_IJ) ,
     &                ABS(TRX0) ), TRX0 )
          TPYDEM(I) = TPYDEM_IJ + SIGN( DMIN1( ABS(TRYDEM_IJ) ,
     &                ABS(TRY0) ), TRY0 )
          TPZDEM(I) = TPZDEM_IJ + SIGN( DMIN1( ABS(TRZDEM_IJ) ,
     &                ABS(TRZ0) ), TRZ0 )      
  
        ENDIF

C---------------------------------------------------------------------------------C
C            assigning variable, chenjun, 20200326                                c
c---------------------------------------------------------------------------------c         
        N_P2(I) = N_P1(I)
        
        IF(N_P1(I) .NE. 0)THEN
          DO J = 1, N_P1(I) 
            N_PP2(I,J) = N_PP1(I,J)  
          ENDDO                        
        ENDIF

        DO J = N_P1(I)+1, I_CPN
          DELTAT_PP1X(I,J) = 0.0
          DELTAT_PP1Y(I,J) = 0.0
          DELTAT_PP1Z(I,J) = 0.0               
        ENDDO   
        
       Endif  
      Enddo
!$OMP END DO      
!$OMP END PARALLEL 
      WRITE(*,'(A14,F8.4)'),'ALFADEM_MAX=',ALFADEM_MAX
      Return
      End