c-----------------------------------------------------------------------------------c
c                             SUBROUTINE PROGRAM                                    #
C                             DESINGED BY ChenJun                                   #
C                                                 SHANGHAI JIAO TONG UNIVERSITY     #
C                                                 SHANGHAI, CHINA                   #
c-----------------------------------------------------------------------------------#
c 子程序功能：
c 1------颗粒接触力和力矩计算，基于背景网格碰撞接触搜索算法     
c-----------------------------------------------------------------------------------#
      
      Subroutine SOFTMODULE2
      Include './Include/OCERM_INF'
c      Parameter (ISCHEME = 2, N_SCHEME1 = 0)      
      Double Precision DIJDEM,PNXDEM,PNYDEM,PNZDEM,PNLDEM
      Double Precision UP_IJ,VP_IJ,WP_IJ,UNDEM,PKNDEM,CNDEM,ALFADEM
      Double Precision FNXDEM_IJ, FNYDEM_IJ, FNZDEM_IJ, FNDEM_IJ
      Double Precision VCTXDEM,VCTYDEM,VCTZDEM,VCTNLDEM
      Double Precision DIATAX,DIATAY,DIATAZ,PKTDEM,CTDEM
      Double Precision FTXDEM_IJ,FTYDEM_IJ,FTZDEM_IJ,FTDEM_IJ
      Double Precision TPXDEM_IJ,TPYDEM_IJ,TPZDEM_IJ,OMDEM    
      Double Precision RSTARDEM,MSTARDEM,ESTARDEM,GSTARDEM            
      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,L1,L2,L3,L4,M,N,I_CELL,
!$OMP&         N_TEMP,N_TEMP1,DIJDEM,PNXDEM,PNYDEM,PNZDEM,
!$OMP&         PNLDEM,UP_IJ,VP_IJ,WP_IJ,UNDEM,PKNDEM,CNDEM,ALFADEM,
!$OMP&         FNXDEM_IJ, FNYDEM_IJ, FNZDEM_IJ, FNDEM_IJ,
!$OMP&         VCTXDEM,VCTYDEM,VCTZDEM,VCTNLDEM,DIATAX,DIATAY,DIATAZ,
!$OMP&         PKTDEM,CTDEM,FTXDEM_IJ,FTYDEM_IJ,FTZDEM_IJ,FTDEM_IJ,
!$OMP&         TPXDEM_IJ,TPYDEM_IJ,TPZDEM_IJ,OMDEM,DIATATN,N_SYM1,TEMP,
!$OMP&         DOMX,DOMY,DOMZ,DOM,TRXDEM_IJ,TRYDEM_IJ,TRZDEM_IJ,
!$OMP&         TRX0,TRY0,TRZ0,VOMX,VOMY,VOMZ,VOM,
!$OMP&         RSTARDEM,MSTARDEM,ESTARDEM,GSTARDEM,
!$OMP&         ALFADEM_MAX,I_ALFADEM_MAX,J_ALFADEM_MAX) 
C!$OMP&         REDUCTION(+:FNXDEM,FNYDEM,FNZDEM,
C!$OMP&         FTXDEM,FTYDEM,FTZDEM,C!$OMP&         TPXDEM,TPYDEM,TPZDEM,
C!$OMP&         TRXDEM,TRYDEM,TRZDEM)

      ALFADEM_MAX = 0.0   
      I_ALFADEM_MAX = 1
      J_ALFADEM_MAX = 1
      
C----------------------------------------------------------------------C
C          低层级颗粒与同层级和高层级颗粒碰撞计算                 
C----------------------------------------------------------------------C
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
        TRXDEM(I) = 0.
        TRYDEM(I) = 0.
        TRZDEM(I) = 0.
        N_P1(I) = 0
        NP_HL(I)=0
        NP_LL(I)=0
           
        OMDEM =  DMAX1( Sqrt(OMXDEM(I) * OMXDEM(I) + 
     &                       OMYDEM(I) * OMYDEM(I) +
     &                       OMZDEM(I) * OMZDEM(I)), acdem)        
   
        DO L4 = I_PSLEV(I), N_PSLEV
        
        DO L3 = MAX(KP(I,L4)-1,1), MIN(KP(I,L4)+1,NBZ(L4))
            
         DO L2 = MAX(JP(I,L4)-1,1), MIN(JP(I,L4)+1,NBY(L4))
             
          DO L1 = MAX(IP(I,L4)-1,1), MIN(IP(I,L4)+1,NBX(L4))   
           
           N_TEMP1 = NBCP_INDEX(L4) + (L3-1)*NBX(L4)*NBY(L4) + 
     &            NBX(L4)*(L2-1) + L1   !颗粒在N层级所属网格单元在NBCP数组中的存储位置   
           
           N_TEMP = NBCP( N_TEMP1 )          
           
           Do N = 1, N_TEMP            
            
            !J = IBCP( (N_TEMP1 - 1) * IBCPMAX + N )
            J = IBCP( IBCP_INDEX(N_TEMP1) + N )
            
            If( I .NE. J ) Then
C----- normal contact displacement 
             IF(J==0)THEN
                  PRINT*,I,J
                  PRINT*,I_PSLEV(I)
                  PRINT*,L1,L2,L3,N
              PAUSE    
             ENDIF
             
             ALFADEM = RPDEM(I) + RPDEM(J) - Sqrt(  
     &                 ( XPDEM(I) - XPDEM(J) )**2 + 
     &                 ( YPDEM(I) - YPDEM(J) )**2 + 
     &                 ( ZPDEM(I) - ZPDEM(J) )**2 ) 
            
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
                 I_ALFADEM_MAX = I
                 J_ALFADEM_MAX = J
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
     &                WP_IJ * PNZDEM   
                    
C----- kn and CNDEM  
              IF(CONMOD .EQ. 'M1')THEN
                        
               PKNDEM = PKN_M1PP
               EN_M1PP = DMAX1(EN_M1PP,1.E-6)
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
               
               MSTARDEM = PMDEM(I) * PMDEM(J) / ( PMDEM(I) + PMDEM(J) )
                  
               PKNDEM = 2. * PI_CON**2 * MSTARDEM / (20.0*TDEM)**2
               CNDEM = 2.0 * 0.3 * Sqrt( MSTARDEM * PKNDEM )
                    
              ELSEIF(CONMOD .EQ. 'M4')THEN
               
               ESTARDEM = 1. / ( (1. - PNIUPDEM(I)**2)/EPDEM(I) +   
     &                           (1. - PNIUPDEM(J)**2)/EPDEM(J) )
                                 
               RSTARDEM = 1. / ( 1./RPDEM(I) + 1./RPDEM(J) )  
               
               MSTARDEM = 1. / ( 1./PMDEM(I) + 1./PMDEM(J) )
                  
               PKNDEM = 4. / 3. * ESTARDEM * SQRT(RSTARDEM * ALFADEM)  
               
               EN_M1PP = DMAX1(EN_M1PP,1.E-6)
               CNDEM = SQRT(5.) * ABS( LOG(EN_M1PP) ) / 
     &                 SQRT( PI_CON**2 + ( LOG(EN_M1PP) )**2 ) *       
     &                 SQRT( MSTARDEM * PKNDEM )               
                    
              ENDIF                        
                                       
C----- fn(i,j)
              FNDEM_IJ =  - PKNDEM * ALFADEM - CNDEM * UNDEM 
                    
              FNXDEM_IJ = FNDEM_IJ * PNXDEM				            
              FNYDEM_IJ = FNDEM_IJ * PNYDEM			           
              FNZDEM_IJ = FNDEM_IJ * PNZDEM    
             
              FNDEM_IJ = ABS( FNDEM_IJ )  		
                    
C----- VT(i,j)
              VCTXDEM = UP_IJ - UNDEM * PNXDEM +
     &			      RPDEM(I) * (OMYDEM(I)*PNZDEM - OMZDEM(I)*PNYDEM) +
     &			      RPDEM(J) * (OMYDEM(J)*PNZDEM - OMZDEM(J)*PNYDEM) 
     	
		    VCTYDEM = VP_IJ - UNDEM * PNYDEM -		
     &			      RPDEM(I) * (OMXDEM(I)*PNZDEM - OMZDEM(I)*PNXDEM) -
     &			      RPDEM(J) * (OMXDEM(J)*PNZDEM - OMZDEM(J)*PNXDEM) 
     			
		    VCTZDEM = WP_IJ - UNDEM * PNZDEM +
     &			      RPDEM(I) * (OMXDEM(I)*PNYDEM - OMYDEM(I)*PNXDEM) +
     &			      RPDEM(J) * (OMXDEM(J)*PNYDEM - OMYDEM(J)*PNXDEM)
     
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

              ELSEIF(CONMOD .EQ. 'M4')THEN   
               
               GSTARDEM = 1. / ( (2. - PNIUPDEM(I))/GPDEM(I) +   
     &                           (2. - PNIUPDEM(J))/GPDEM(J) )
                  
               PKTDEM = 8. * GSTARDEM * SQRT(RSTARDEM * ALFADEM)
               CTDEM = SQRT(10./3.) * ABS( LOG(EN_M1PP) ) / 
     &                 SQRT( PI_CON**2 + ( LOG(EN_M1PP) )**2 ) *  
     &                 SQRT( MSTARDEM * PKTDEM )               
               
              ENDIF                    
                      
C----- ft(i,j) 
              FTXDEM_IJ = - PKTDEM * DIATAX - CTDEM * VCTXDEM 
              FTYDEM_IJ = - PKTDEM * DIATAY - CTDEM * VCTYDEM
              FTZDEM_IJ = - PKTDEM * DIATAZ - CTDEM * VCTZDEM
                    
              FTDEM_IJ =  Sqrt( FTXDEM_IJ * FTXDEM_IJ + 
     &                          FTYDEM_IJ * FTYDEM_IJ +
     &                          FTZDEM_IJ * FTZDEM_IJ )           
                    
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

C              TPXDEM_IJ = TPXDEM_IJ + RPDEM(I) *  
C     &                    ( PNYDEM * FTZDEM_IJ - PNZDEM * FTYDEM_IJ )   
                    
C              TPYDEM_IJ = TPYDEM_IJ + RPDEM(I) *  
C     &                    ( PNZDEM * FTXDEM_IJ - PNXDEM * FTZDEM_IJ )   
                   
C		    TPZDEM_IJ = TPZDEM_IJ + RPDEM(I) *  
C     &                    ( PNXDEM * FTYDEM_IJ - PNYDEM * FTXDEM_IJ )   
                   
              TPXDEM_IJ = RPDEM(I) *  
     &                    ( PNYDEM * FTZDEM_IJ - PNZDEM * FTYDEM_IJ )   
                    
              TPYDEM_IJ = RPDEM(I) *  
     &                    ( PNZDEM * FTXDEM_IJ - PNXDEM * FTZDEM_IJ )   
                   
		    TPZDEM_IJ = RPDEM(I) *  
     &                    ( PNXDEM * FTYDEM_IJ - PNYDEM * FTXDEM_IJ ) 

C----- rolling friction torque               
              IF(ROLLMOD .EQ. 'CON1')THEN
                 
C               TRXDEM_IJ = TRXDEM_IJ - RMIUS * FNDEM_IJ *
C     &                                 OMXDEM(I) / OMDEM
               
C               TRYDEM_IJ = TRYDEM_IJ - RMIUS * FNDEM_IJ * 
C     &                                 OMYDEM(I) / OMDEM    
               
C               TRZDEM_IJ = TRZDEM_IJ - RMIUS * FNDEM_IJ * 
C     &                                 OMZDEM(I) / OMDEM 
               
               TRXDEM_IJ = - RMIUS * FNDEM_IJ * OMXDEM(I) / OMDEM
               
               TRYDEM_IJ = - RMIUS * FNDEM_IJ * OMYDEM(I) / OMDEM    
               
               TRZDEM_IJ = - RMIUS * FNDEM_IJ * OMZDEM(I) / OMDEM       
               
              ELSEIF(ROLLMOD .EQ. 'CON2')THEN
                 
               DOMX = OMXDEM(I) - OMXDEM(J)
               DOMY = OMYDEM(I) - OMYDEM(J)
               DOMZ = OMZDEM(I) - OMZDEM(J)
               
               DOM =  DMAX1( Sqrt( DOMX**2 + DOMY**2 + DOMZ**2 ), acdem)
               
C               TRXDEM_IJ = TRXDEM_IJ -  RMIUS * FNDEM_IJ * DOMX / DOM
C               TRYDEM_IJ = TRYDEM_IJ -  RMIUS * FNDEM_IJ * DOMY / DOM
C               TRZDEM_IJ = TRZDEM_IJ -  RMIUS * FNDEM_IJ * DOMZ / DOM
               
               TRXDEM_IJ = -  RMIUS * FNDEM_IJ * DOMX / DOM
               TRYDEM_IJ = -  RMIUS * FNDEM_IJ * DOMY / DOM
               TRZDEM_IJ = -  RMIUS * FNDEM_IJ * DOMZ / DOM           
               
              ELSEIF(ROLLMOD .EQ. 'VIS1')THEN
                 
               VOMX = RPDEM(I) * ( OMYDEM(I)*PNZDEM-OMZDEM(I)*PNYDEM ) +
     &                RPDEM(J) * ( OMYDEM(J)*PNZDEM-OMZDEM(J)*PNYDEM )
               
               VOMY = RPDEM(I) * ( OMZDEM(I)*PNXDEM-OMXDEM(I)*PNZDEM ) +
     &                RPDEM(J) * ( OMZDEM(J)*PNXDEM-OMXDEM(J)*PNZDEM )
               
               VOMZ = RPDEM(I) * ( OMXDEM(I)*PNYDEM-OMYDEM(I)*PNXDEM ) +
     &                RPDEM(J) * ( OMXDEM(J)*PNYDEM-OMYDEM(J)*PNXDEM )
               
               VOM = SQRT( VOMX**2 + VOMY**2 + VOMZ**2 )
               
C               TRXDEM_IJ = TRXDEM_IJ - RMIUS * VOM *
C     &                     FNDEM_IJ * OMXDEM(I) / OMDEM
               
C               TRYDEM_IJ = TRYDEM_IJ - RMIUS * VOM *
C     &                     FNDEM_IJ * OMYDEM(I) / OMDEM    
               
C               TRZDEM_IJ = TRZDEM_IJ - RMIUS * VOM *
C     &                     FNDEM_IJ * OMZDEM(I) / OMDEM    
               
               TRXDEM_IJ = - RMIUS * VOM * FNDEM_IJ * OMXDEM(I) / OMDEM
               
               TRYDEM_IJ = - RMIUS * VOM * FNDEM_IJ * OMYDEM(I) / OMDEM
               
               TRZDEM_IJ = - RMIUS * VOM * FNDEM_IJ * OMZDEM(I) / OMDEM
               
              ELSEIF(ROLLMOD .EQ. 'VIS2')THEN      
                 
               VOMX = RPDEM(I) * ( OMYDEM(I)*PNZDEM-OMZDEM(I)*PNYDEM ) +
     &                RPDEM(J) * ( OMYDEM(J)*PNZDEM-OMZDEM(J)*PNYDEM )
               
               VOMY = RPDEM(I) * ( OMZDEM(I)*PNXDEM-OMXDEM(I)*PNZDEM ) +
     &                RPDEM(J) * ( OMZDEM(J)*PNXDEM-OMXDEM(J)*PNZDEM )
               
               VOMZ = RPDEM(I) * ( OMXDEM(I)*PNYDEM-OMYDEM(I)*PNXDEM ) +
     &                RPDEM(J) * ( OMXDEM(J)*PNYDEM-OMYDEM(J)*PNXDEM )  

C               TRXDEM_IJ = TRXDEM_IJ -  RMIUS * VOMX * FNDEM_IJ 
C               TRYDEM_IJ = TRYDEM_IJ -  RMIUS * VOMY * FNDEM_IJ
C               TRZDEM_IJ = TRZDEM_IJ -  RMIUS * VOMZ * FNDEM_IJ      

               TRXDEM_IJ = -  RMIUS * VOMX * FNDEM_IJ 
               TRYDEM_IJ = -  RMIUS * VOMY * FNDEM_IJ
               TRZDEM_IJ = -  RMIUS * VOMZ * FNDEM_IJ   
               
              ENDIF             
                    
              FNXDEM(I) = FNXDEM(I) + FNXDEM_IJ                 
              FNYDEM(I) = FNYDEM(I) + FNYDEM_IJ  
              FNZDEM(I) = FNZDEM(I) + FNZDEM_IJ
                    
              FTXDEM(I) = FTXDEM(I) + FTXDEM_IJ
              FTYDEM(I) = FTYDEM(I) + FTYDEM_IJ
              FTZDEM(I) = FTZDEM(I) + FTZDEM_IJ
              
              TPXDEM(I) = TPXDEM(I) + TPXDEM_IJ
              TPYDEM(I) = TPYDEM(I) + TPYDEM_IJ
              TPZDEM(I) = TPZDEM(I) + TPZDEM_IJ
              
              TRXDEM(I) = TRXDEM(I) + TRXDEM_IJ
              TRYDEM(I) = TRYDEM(I) + TRYDEM_IJ
              TRZDEM(I) = TRZDEM(I) + TRZDEM_IJ              
              
              IF(L4 .GT. I_PSLEV(I))THEN   
                  
                NP_HL(I) = NP_HL(I) + 1   !与低层级颗粒接触的高层级颗粒数
                
                IF(NP_HL(I) .GT. 10)THEN
                    
                  PRINT*,'Error of NP_HL(I) in SOFTMODULE2.F!'
                  WRITE(IUPRT,'(2A14)') 'I','NP_HL(I)'
                  WRITE(IUPRT,'(2I14)') I,NP_HL(I)
                  WRITE(IUPRT,'(3A14)') 'XPDEM', 'YPDEM', 'ZPDEM' 
                  WRITE(IUPRT,'(3F14.7)') XPDEM(I),YPDEM(I),ZPDEM(I)
                  PAUSE
                  
                ENDIF  
                    
                IP_HL(I,NP_HL(I)) = J   !与低层级颗粒接触的高层级颗粒编号
                
              ENDIF                          
                                   
             ENDIF
            ENDIF       
           ENDDO 
          ENDDO
         ENDDO
        ENDDO
        ENDDO
           
       ENDIF          
      ENDDO
!$OMP END DO 

!$OMP BARRIER

!$OMP MASTER
      Do I = 1 , I_PARTICLE
       IF(NSYM_CONT(I) .EQ. 1)THEN          
        DO J = 1, NP_HL(I)   !高层级颗粒间循环
            
          NP_LL(IP_HL(I,J)) = NP_LL(IP_HL(I,J)) + 1   !与高层级颗粒接触的低层级颗粒数

          IF(NP_LL(IP_HL(I,J)) .GT. NPT_LL(IP_HL(I,J)))THEN
                    
            PRINT*,'Error of NP_LL(IP_HL(I,J)) in SOFTMODULE2.F!'
            WRITE(IUPRT,'(A14)') 'NP_LL(IP_HL(I,J))'
            WRITE(IUPRT,'(I14)') NP_LL(IP_HL(I,J))
            WRITE(IUPRT,'(A14)') 'NPT_LL(I)' 
            WRITE(IUPRT,'(I14)') NPT_LL(I)
            PAUSE
                  
          ENDIF            
          
          N_TEMP = IP_LL_INDEX(IP_HL(I,J)) + NP_LL(IP_HL(I,J))
          
          IP_LL(N_TEMP) = I   !与高层级颗粒接触的低层级颗粒编号          
          
        ENDDO
       ENDIF 
      ENDDO
!$OMP END MASTER      

!$OMP BARRIER

!$OMP DO
      Do I = 1 , I_PARTICLE
          
       IF(NSYM_CONT(I) .EQ. 1)THEN 
           
C----------------------------------------------------------------------C
C          高层级颗粒与低层级颗粒碰撞计算                 
C----------------------------------------------------------------------C 

         DO N_TEMP = 1, NP_LL(I)
           
          J = IP_LL( IP_LL_INDEX(I) + N_TEMP )   
         
          If( I .NE. J ) Then
C----- normal contact displacement 
             IF(J==0)THEN
                  PRINT*,I,J
                  PRINT*,I_PSLEV(I)
                  PRINT*,L1,L2,L3,N
              PAUSE    
             ENDIF
             
             ALFADEM = RPDEM(I) + RPDEM(J) - Sqrt(  
     &                 ( XPDEM(I) - XPDEM(J) )**2 + 
     &                 ( YPDEM(I) - YPDEM(J) )**2 + 
     &                 ( ZPDEM(I) - ZPDEM(J) )**2 ) 
            
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
                 I_ALFADEM_MAX = I
                 J_ALFADEM_MAX = J
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
     &                WP_IJ * PNZDEM   
                    
C----- kn and CNDEM  
              IF(CONMOD .EQ. 'M1')THEN
                        
               PKNDEM = PKN_M1PP
               EN_M1PP = DMAX1(EN_M1PP,1.E-6)
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
               
               MSTARDEM = PMDEM(I) * PMDEM(J) / ( PMDEM(I) + PMDEM(J) )
                  
               PKNDEM = 2. * PI_CON**2 * MSTARDEM / (20.0*TDEM)**2
               CNDEM = 2.0 * 0.3 * Sqrt( MSTARDEM * PKNDEM )
                    
              ELSEIF(CONMOD .EQ. 'M4')THEN
               
               ESTARDEM = 1. / ( (1. - PNIUPDEM(I)**2)/EPDEM(I) +   
     &                           (1. - PNIUPDEM(J)**2)/EPDEM(J) )
                                 
               RSTARDEM = 1. / ( 1./RPDEM(I) + 1./RPDEM(J) )  
               
               MSTARDEM = 1. / ( 1./PMDEM(I) + 1./PMDEM(J) )
                  
               PKNDEM = 4. / 3. * ESTARDEM * SQRT(RSTARDEM * ALFADEM)  
               
               EN_M1PP = DMAX1(EN_M1PP,1.E-6)
               CNDEM = SQRT(5.) * ABS( LOG(EN_M1PP) ) / 
     &                 SQRT( PI_CON**2 + ( LOG(EN_M1PP) )**2 ) *       
     &                 SQRT( MSTARDEM * PKNDEM )               
                    
              ENDIF                        
                                       
C----- fn(i,j)
              FNDEM_IJ =  - PKNDEM * ALFADEM - CNDEM * UNDEM 
                    
              FNXDEM_IJ = FNDEM_IJ * PNXDEM				            
              FNYDEM_IJ = FNDEM_IJ * PNYDEM			           
              FNZDEM_IJ = FNDEM_IJ * PNZDEM    
             
              FNDEM_IJ = ABS( FNDEM_IJ )  		
                    
C----- VT(i,j)
              VCTXDEM = UP_IJ - UNDEM * PNXDEM +
     &			      RPDEM(I) * (OMYDEM(I)*PNZDEM - OMZDEM(I)*PNYDEM) +
     &			      RPDEM(J) * (OMYDEM(J)*PNZDEM - OMZDEM(J)*PNYDEM) 
     	
		    VCTYDEM = VP_IJ - UNDEM * PNYDEM -		
     &			      RPDEM(I) * (OMXDEM(I)*PNZDEM - OMZDEM(I)*PNXDEM) -
     &			      RPDEM(J) * (OMXDEM(J)*PNZDEM - OMZDEM(J)*PNXDEM) 
     			
		    VCTZDEM = WP_IJ - UNDEM * PNZDEM +
     &			      RPDEM(I) * (OMXDEM(I)*PNYDEM - OMYDEM(I)*PNXDEM) +
     &			      RPDEM(J) * (OMXDEM(J)*PNYDEM - OMYDEM(J)*PNXDEM)
     
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

              ELSEIF(CONMOD .EQ. 'M4')THEN   
               
               GSTARDEM = 1. / ( (2. - PNIUPDEM(I))/GPDEM(I) +   
     &                           (2. - PNIUPDEM(J))/GPDEM(J) )
                  
               PKTDEM = 8. * GSTARDEM * SQRT(RSTARDEM * ALFADEM)
               CTDEM = SQRT(10./3.) * ABS( LOG(EN_M1PP) ) / 
     &                 SQRT( PI_CON**2 + ( LOG(EN_M1PP) )**2 ) *  
     &                 SQRT( MSTARDEM * PKTDEM )               
               
              ENDIF                    
                      
C----- ft(i,j) 
              FTXDEM_IJ = - PKTDEM * DIATAX - CTDEM * VCTXDEM 
              FTYDEM_IJ = - PKTDEM * DIATAY - CTDEM * VCTYDEM
              FTZDEM_IJ = - PKTDEM * DIATAZ - CTDEM * VCTZDEM
                    
              FTDEM_IJ =  Sqrt( FTXDEM_IJ * FTXDEM_IJ + 
     &                          FTYDEM_IJ * FTYDEM_IJ +
     &                          FTZDEM_IJ * FTZDEM_IJ )           
                    
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

C              TPXDEM_IJ = TPXDEM_IJ + RPDEM(I) *  
C     &                    ( PNYDEM * FTZDEM_IJ - PNZDEM * FTYDEM_IJ )   
                    
C              TPYDEM_IJ = TPYDEM_IJ + RPDEM(I) *  
C     &                    ( PNZDEM * FTXDEM_IJ - PNXDEM * FTZDEM_IJ )   
                   
C		    TPZDEM_IJ = TPZDEM_IJ + RPDEM(I) *  
C     &                    ( PNXDEM * FTYDEM_IJ - PNYDEM * FTXDEM_IJ )   
                   
              TPXDEM_IJ = RPDEM(I) *  
     &                    ( PNYDEM * FTZDEM_IJ - PNZDEM * FTYDEM_IJ )   
                    
              TPYDEM_IJ = RPDEM(I) *  
     &                    ( PNZDEM * FTXDEM_IJ - PNXDEM * FTZDEM_IJ )   
                   
		    TPZDEM_IJ = RPDEM(I) *  
     &                    ( PNXDEM * FTYDEM_IJ - PNYDEM * FTXDEM_IJ ) 

C----- rolling friction torque               
              IF(ROLLMOD .EQ. 'CON1')THEN
                 
C               TRXDEM_IJ = TRXDEM_IJ - RMIUS * FNDEM_IJ *
C     &                                 OMXDEM(I) / OMDEM
               
C               TRYDEM_IJ = TRYDEM_IJ - RMIUS * FNDEM_IJ * 
C     &                                 OMYDEM(I) / OMDEM    
               
C               TRZDEM_IJ = TRZDEM_IJ - RMIUS * FNDEM_IJ * 
C     &                                 OMZDEM(I) / OMDEM 
               
               TRXDEM_IJ = - RMIUS * FNDEM_IJ * OMXDEM(I) / OMDEM
               
               TRYDEM_IJ = - RMIUS * FNDEM_IJ * OMYDEM(I) / OMDEM    
               
               TRZDEM_IJ = - RMIUS * FNDEM_IJ * OMZDEM(I) / OMDEM       
               
              ELSEIF(ROLLMOD .EQ. 'CON2')THEN
                 
               DOMX = OMXDEM(I) - OMXDEM(J)
               DOMY = OMYDEM(I) - OMYDEM(J)
               DOMZ = OMZDEM(I) - OMZDEM(J)
               
               DOM =  DMAX1( Sqrt( DOMX**2 + DOMY**2 + DOMZ**2 ), acdem)
               
C               TRXDEM_IJ = TRXDEM_IJ -  RMIUS * FNDEM_IJ * DOMX / DOM
C               TRYDEM_IJ = TRYDEM_IJ -  RMIUS * FNDEM_IJ * DOMY / DOM
C               TRZDEM_IJ = TRZDEM_IJ -  RMIUS * FNDEM_IJ * DOMZ / DOM
               
               TRXDEM_IJ = -  RMIUS * FNDEM_IJ * DOMX / DOM
               TRYDEM_IJ = -  RMIUS * FNDEM_IJ * DOMY / DOM
               TRZDEM_IJ = -  RMIUS * FNDEM_IJ * DOMZ / DOM           
               
              ELSEIF(ROLLMOD .EQ. 'VIS1')THEN
                 
               VOMX = RPDEM(I) * ( OMYDEM(I)*PNZDEM-OMZDEM(I)*PNYDEM ) +
     &                RPDEM(J) * ( OMYDEM(J)*PNZDEM-OMZDEM(J)*PNYDEM )
               
               VOMY = RPDEM(I) * ( OMZDEM(I)*PNXDEM-OMXDEM(I)*PNZDEM ) +
     &                RPDEM(J) * ( OMZDEM(J)*PNXDEM-OMXDEM(J)*PNZDEM )
               
               VOMZ = RPDEM(I) * ( OMXDEM(I)*PNYDEM-OMYDEM(I)*PNXDEM ) +
     &                RPDEM(J) * ( OMXDEM(J)*PNYDEM-OMYDEM(J)*PNXDEM )
               
               VOM = SQRT( VOMX**2 + VOMY**2 + VOMZ**2 )
               
C               TRXDEM_IJ = TRXDEM_IJ - RMIUS * VOM *
C     &                     FNDEM_IJ * OMXDEM(I) / OMDEM
               
C               TRYDEM_IJ = TRYDEM_IJ - RMIUS * VOM *
C     &                     FNDEM_IJ * OMYDEM(I) / OMDEM    
               
C               TRZDEM_IJ = TRZDEM_IJ - RMIUS * VOM *
C     &                     FNDEM_IJ * OMZDEM(I) / OMDEM    
               
               TRXDEM_IJ = - RMIUS * VOM * FNDEM_IJ * OMXDEM(I) / OMDEM
               
               TRYDEM_IJ = - RMIUS * VOM * FNDEM_IJ * OMYDEM(I) / OMDEM
               
               TRZDEM_IJ = - RMIUS * VOM * FNDEM_IJ * OMZDEM(I) / OMDEM
               
              ELSEIF(ROLLMOD .EQ. 'VIS2')THEN      
                 
               VOMX = RPDEM(I) * ( OMYDEM(I)*PNZDEM-OMZDEM(I)*PNYDEM ) +
     &                RPDEM(J) * ( OMYDEM(J)*PNZDEM-OMZDEM(J)*PNYDEM )
               
               VOMY = RPDEM(I) * ( OMZDEM(I)*PNXDEM-OMXDEM(I)*PNZDEM ) +
     &                RPDEM(J) * ( OMZDEM(J)*PNXDEM-OMXDEM(J)*PNZDEM )
               
               VOMZ = RPDEM(I) * ( OMXDEM(I)*PNYDEM-OMYDEM(I)*PNXDEM ) +
     &                RPDEM(J) * ( OMXDEM(J)*PNYDEM-OMYDEM(J)*PNXDEM )  

C               TRXDEM_IJ = TRXDEM_IJ -  RMIUS * VOMX * FNDEM_IJ 
C               TRYDEM_IJ = TRYDEM_IJ -  RMIUS * VOMY * FNDEM_IJ
C               TRZDEM_IJ = TRZDEM_IJ -  RMIUS * VOMZ * FNDEM_IJ      

               TRXDEM_IJ = -  RMIUS * VOMX * FNDEM_IJ 
               TRYDEM_IJ = -  RMIUS * VOMY * FNDEM_IJ
               TRZDEM_IJ = -  RMIUS * VOMZ * FNDEM_IJ   
               
              ENDIF             
                    
              FNXDEM(I) = FNXDEM(I) + FNXDEM_IJ                 
              FNYDEM(I) = FNYDEM(I) + FNYDEM_IJ  
              FNZDEM(I) = FNZDEM(I) + FNZDEM_IJ
                    
              FTXDEM(I) = FTXDEM(I) + FTXDEM_IJ
              FTYDEM(I) = FTYDEM(I) + FTYDEM_IJ
              FTZDEM(I) = FTZDEM(I) + FTZDEM_IJ
              
              TPXDEM(I) = TPXDEM(I) + TPXDEM_IJ
              TPYDEM(I) = TPYDEM(I) + TPYDEM_IJ
              TPZDEM(I) = TPZDEM(I) + TPZDEM_IJ
              
              TRXDEM(I) = TRXDEM(I) + TRXDEM_IJ
              TRYDEM(I) = TRYDEM(I) + TRYDEM_IJ
              TRZDEM(I) = TRZDEM(I) + TRZDEM_IJ              
                
             ENDIF
             
          ENDIF

         ENDDO          
                                     
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
                
C              ALFADEM = 2 * RPDEM(I)- PNLDEM 
              ALFADEM = DEL_DEM(I,J)  
                
C              If(ALFADEM .GT. 0.1*2*RPDEM(I)) Then
C                Print*, 'Error1 in softmodule.f: contact error!' 
C                print*, 'ALFADEM=',ALFADEM 
C                print*, '(X,Y,Z)=',XPDEM(I),YPDEM(I),ZPDEM(I)
C                !Pause
C              Endif    
                
              IF(ALFADEM / RPDEM(I) .GT. ALFADEM_MAX)THEN
                 ALFADEM_MAX = ALFADEM / RPDEM(I)
                 I_ALFADEM_MAX = I
                 J_ALFADEM_MAX = I
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
                Pause
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
                    
              ELSEIF(CONMOD .EQ. 'M4')THEN
               
               ESTARDEM = 1. / ( (1. - PNIUPDEM(I)**2)/EPDEM(I) +   
     &                           (1. - PNIUPDEM(I)**2)/EPDEM(I) )
                                 
               RSTARDEM = 1. / ( 1./RPDEM(I) + 1./RPDEM(I) )  
               
               MSTARDEM = 1. / ( 1./PMDEM(I) + 1./PMDEM(I) )
                  
               PKNDEM = 4. / 3. * ESTARDEM * SQRT(RSTARDEM * ALFADEM)  
               
               CNDEM = SQRT(5.) * ABS( LOG(EN_M1PW) ) / 
     &                 SQRT( PI_CON**2 + ( LOG(EN_M1PW) )**2 ) *       
     &                 SQRT( MSTARDEM * PKNDEM )   
               
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
 
              ELSEIF(CONMOD .EQ. 'M4')THEN   
               
               GSTARDEM = 1. / ( (2. - PNIUPDEM(I))/GPDEM(I) +   
     &                           (2. - PNIUPDEM(I))/GPDEM(I) )
                  
               PKTDEM = 8. * GSTARDEM * SQRT(RSTARDEM * ALFADEM)
               CTDEM = SQRT(10./3.) * ABS( LOG(EN_M1PW) ) / 
     &                 SQRT( PI_CON**2 + ( LOG(EN_M1PW) )**2 ) *  
     &                 SQRT( MSTARDEM * PKTDEM )                 
                
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
                                                  
C			TPXDEM_IJ = TPXDEM_IJ + RPDEM(I) *  
C     &                    ( PNYDEM * FTZDEM_IJ - PNZDEM * FTYDEM_IJ )  
                    
C			TPYDEM_IJ = TPYDEM_IJ + RPDEM(I) *  
C     &                    ( PNZDEM * FTXDEM_IJ - PNXDEM * FTZDEM_IJ )   
                    
C              TPZDEM_IJ = TPZDEM_IJ + RPDEM(I) *  
C     &                    ( PNXDEM * FTYDEM_IJ - PNYDEM * FTXDEM_IJ ) 
                    
			TPXDEM(I) = TPXDEM(I) + RPDEM(I) *  
     &                    ( PNYDEM * FTZDEM_IJ - PNZDEM * FTYDEM_IJ )  
                    
			TPYDEM(I) = TPYDEM(I) + RPDEM(I) *  
     &                    ( PNZDEM * FTXDEM_IJ - PNXDEM * FTZDEM_IJ )   
                    
              TPZDEM(I) = TPZDEM(I) + RPDEM(I) *  
     &                    ( PNXDEM * FTYDEM_IJ - PNYDEM * FTXDEM_IJ ) 

C----- rolling friction torque             
              IF(ROLLMOD .EQ. 'CON1')THEN  
                 
C                TRXDEM_IJ = TRXDEM_IJ - RMIUSPW * FNDEM_IJ *
C     &                                  OMXDEM(I) / OMDEM
               
C                TRYDEM_IJ = TRYDEM_IJ - RMIUSPW * FNDEM_IJ * 
C     &                                  OMYDEM(I) / OMDEM    
               
C                TRZDEM_IJ = TRZDEM_IJ - RMIUSPW * FNDEM_IJ * 
C     &                                  OMZDEM(I) / OMDEM 
               
                TRXDEM(I) = TRXDEM(I) - RMIUSPW * FNDEM_IJ *
     &                                  OMXDEM(I) / OMDEM
               
                TRYDEM(I) = TRYDEM(I) - RMIUSPW * FNDEM_IJ * 
     &                                  OMYDEM(I) / OMDEM    
               
                TRZDEM(I) = TRZDEM(I) - RMIUSPW * FNDEM_IJ * 
     &                                  OMZDEM(I) / OMDEM             
               
              ELSEIF(ROLLMOD .EQ. 'CON2')THEN
                                
C                TRXDEM_IJ = TRXDEM_IJ -  RMIUSPW * FNDEM_IJ * 
C     &                                   OMXDEM(I) / OMDEM
C                TRYDEM_IJ = TRYDEM_IJ -  RMIUSPW * FNDEM_IJ * 
C     &                                   OMYDEM(I) / OMDEM
C                TRZDEM_IJ = TRZDEM_IJ -  RMIUSPW * FNDEM_IJ * 
C     &                                   OMZDEM(I) / OMDEM
                                
                TRXDEM(I) = TRXDEM(I) -  RMIUSPW * FNDEM_IJ * 
     &                                   OMXDEM(I) / OMDEM
                TRYDEM(I) = TRYDEM(I) -  RMIUSPW * FNDEM_IJ * 
     &                                   OMYDEM(I) / OMDEM
                TRZDEM(I) = TRZDEM(I) -  RMIUSPW * FNDEM_IJ * 
     &                                   OMZDEM(I) / OMDEM             
               
              ELSEIF(ROLLMOD .EQ. 'VIS1')THEN
                 
                VOMX = RPDEM(I) * ( OMYDEM(I) * PNZDEM - 
     &                              OMZDEM(I) * PNYDEM )
               
                VOMY = RPDEM(I) * ( OMZDEM(I) * PNXDEM - 
     &                              OMXDEM(I) * PNZDEM )
               
                VOMZ = RPDEM(I) * ( OMXDEM(I) * PNYDEM - 
     &                              OMYDEM(I) * PNXDEM )
                  
                VOM = SQRT( VOMX**2 + VOMY**2 + VOMZ**2 )

C                TRXDEM_IJ = TRXDEM_IJ - RMIUSPW * VOM *
C     &                      FNDEM_IJ * OMXDEM(I) / OMDEM
               
C                TRYDEM_IJ = TRYDEM_IJ - RMIUSPW * VOM *
C     &                      FNDEM_IJ * OMYDEM(I) / OMDEM    
               
C                TRZDEM_IJ = TRZDEM_IJ - RMIUSPW * VOM *
C     &                      FNDEM_IJ * OMZDEM(I) / OMDEM 
               
                TRXDEM(I) = TRXDEM(I) - RMIUSPW * VOM *
     &                      FNDEM_IJ * OMXDEM(I) / OMDEM
               
                TRYDEM(I) = TRYDEM(I) - RMIUSPW * VOM *
     &                      FNDEM_IJ * OMYDEM(I) / OMDEM    
               
                TRZDEM(I) = TRZDEM(I) - RMIUSPW * VOM *
     &                      FNDEM_IJ * OMZDEM(I) / OMDEM                
                              
              ELSEIF(ROLLMOD .EQ. 'VIS2')THEN      
                 
                VOMX = RPDEM(I) * ( OMYDEM(I) * PNZDEM - 
     &                                OMZDEM(I) * PNYDEM )
               
                VOMY = RPDEM(I) * ( OMZDEM(I) * PNXDEM - 
     &                                OMXDEM(I) * PNZDEM )
               
                VOMZ = RPDEM(I) * ( OMXDEM(I) * PNYDEM - 
     &                                OMYDEM(I) * PNXDEM )

C                TRXDEM_IJ = TRXDEM_IJ -  RMIUSPW * VOMX * FNDEM_IJ 
C                TRYDEM_IJ = TRYDEM_IJ -  RMIUSPW * VOMY * FNDEM_IJ
C                TRZDEM_IJ = TRZDEM_IJ -  RMIUSPW * VOMZ * FNDEM_IJ  

                TRXDEM(I) = TRXDEM(I) -  RMIUSPW * VOMX * FNDEM_IJ 
                TRYDEM(I) = TRYDEM(I) -  RMIUSPW * VOMY * FNDEM_IJ
                TRZDEM(I) = TRZDEM(I) -  RMIUSPW * VOMZ * FNDEM_IJ 
                
              ELSEIF(ROLLMOD .EQ. 'NONE')THEN
                  
C                TRXDEM_IJ = 0.0 
C                TRYDEM_IJ = 0.0
C                TRZDEM_IJ = 0.0    
               
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
             
C          TPXDEM(I) = TPXDEM_IJ + TRXDEM_IJ
C          TPYDEM(I) = TPYDEM_IJ + TRYDEM_IJ
C          TPZDEM(I) = TPZDEM_IJ + TRZDEM_IJ
             
          TPXDEM(I) = TPXDEM(I) + TRXDEM(I)
          TPYDEM(I) = TPYDEM(I) + TRYDEM(I)
          TPZDEM(I) = TPZDEM(I) + TRZDEM(I)             
     
        ELSEIF(ROLLBOUND .EQ. 'BOD1')THEN 

          TRX0 = - PIPDEM(I) * OMXDEM(I) / TDEM 
          TRY0 = - PIPDEM(I) * OMYDEM(I) / TDEM
          TRZ0 = - PIPDEM(I) * OMZDEM(I) / TDEM
           
C          TPXDEM(I) = TPXDEM_IJ + SIGN( DMIN1( ABS(TRXDEM_IJ) ,
C     &                ABS(TRX0) ), TRX0 )
C          TPYDEM(I) = TPYDEM_IJ + SIGN( DMIN1( ABS(TRYDEM_IJ) ,
C     &                ABS(TRY0) ), TRY0 )
C          TPZDEM(I) = TPZDEM_IJ + SIGN( DMIN1( ABS(TRZDEM_IJ) ,
C     &                ABS(TRZ0) ), TRZ0 ) 
           
          TPXDEM(I) = TPXDEM(I) + SIGN( DMIN1( ABS(TRXDEM(I)) ,
     &                ABS(TRX0) ), TRX0 )
          TPYDEM(I) = TPYDEM(I) + SIGN( DMIN1( ABS(TRYDEM(I)) ,
     &                ABS(TRY0) ), TRY0 )
          TPZDEM(I) = TPZDEM(I) + SIGN( DMIN1( ABS(TRZDEM(I)) ,
     &                ABS(TRZ0) ), TRZ0 )            
      
        ELSEIF(ROLLBOUND .EQ. 'BOD2')THEN
             
C          TRX0 = - PIPDEM(I) * OMXDEM(I) / TDEM - TPXDEM_IJ
C          TRY0 = - PIPDEM(I) * OMYDEM(I) / TDEM - TPYDEM_IJ
C          TRZ0 = - PIPDEM(I) * OMZDEM(I) / TDEM - TPZDEM_IJ
             
          TRX0 = - PIPDEM(I) * OMXDEM(I) / TDEM - TPXDEM(I)
          TRY0 = - PIPDEM(I) * OMYDEM(I) / TDEM - TPYDEM(I)
          TRZ0 = - PIPDEM(I) * OMZDEM(I) / TDEM - TPZDEM(I)             
           
C          TPXDEM(I) = TPXDEM_IJ + SIGN( DMIN1( ABS(TRXDEM_IJ) ,
C     &                ABS(TRX0) ), TRX0 )
C          TPYDEM(I) = TPYDEM_IJ + SIGN( DMIN1( ABS(TRYDEM_IJ) ,
C     &                ABS(TRY0) ), TRY0 )
C          TPZDEM(I) = TPZDEM_IJ + SIGN( DMIN1( ABS(TRZDEM_IJ) ,
C     &                ABS(TRZ0) ), TRZ0 )    
           
          TPXDEM(I) = TPXDEM(I) + SIGN( DMIN1( ABS(TRXDEM(I)) ,
     &                ABS(TRX0) ), TRX0 )
          TPYDEM(I) = TPYDEM(I) + SIGN( DMIN1( ABS(TRYDEM(I)) ,
     &                ABS(TRY0) ), TRY0 )
          TPZDEM(I) = TPZDEM(I) + SIGN( DMIN1( ABS(TRZDEM(I)) ,
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

        !DO J = N_P1(I)+1, I_CPN   !删除？
        !  DELTAT_PP1X(I,J) = 0.0
        !  DELTAT_PP1Y(I,J) = 0.0
        !  DELTAT_PP1Z(I,J) = 0.0               
        !ENDDO 
       
       ENDIF

      ENDDO         
!$OMP END DO          

C!$OMP MASTER
c      WRITE(*,'(A14,F8.4)'),'ALFADEM_MAX =',ALFADEM_MAX
c      PRINT*,I_ALFADEM_MAX,J_ALFADEM_MAX
c      WRITE(*,'(3F8.4)'),XPDEM(I_ALFADEM_MAX),YPDEM(I_ALFADEM_MAX),
c     &       ZPDEM(I_ALFADEM_MAX)
c      WRITE(*,'(3F8.4)'),XPDEM(J_ALFADEM_MAX),YPDEM(J_ALFADEM_MAX),
c     &       ZPDEM(J_ALFADEM_MAX) 
C!$OMP END MASTER

!$OMP END PARALLEL 
    
      Return
      End
      
C----------------------------------------------------------------------C
C          存储颗粒所在层级背景网格及相邻网格颗粒数目及颗粒编号                 
C----------------------------------------------------------------------C
      !Do I = 1 , I_PARTICLE
      !    
      ! IF(NSYM_CONT(I) .EQ. 1)THEN  
      !     
      !  NP_CAND1(I) = 0   
      !
      !  DO L4 = I_PSLEV(I), N_PSLEV
      !  
      !  DO L3 = MAX(KP(I,L4)-1,1), MIN(KP(I,L4)+1,NBZ(L4))
      !      
      !   DO L2 = MAX(JP(I,L4)-1,1), MIN(JP(I,L4)+1,NBY(L4))
      !       
      !    DO L1 = MAX(IP(I,L4)-1,1), MIN(IP(I,L4)+1,NBX(L4))  
      !        
      !     IF(L4 .EQ. 1)THEN   
      !        N_TEMP = NBCP1(L1,L2,L3)
      !     ELSEIF(L4 .EQ. 2)THEN   
      !        N_TEMP = NBCP2(L1,L2,L3)
      !     ELSEIF(L4 .EQ. 3)THEN   
      !        N_TEMP = NBCP3(L1,L2,L3)   
      !     ELSEIF(L4 .EQ. 4)THEN   
      !        N_TEMP = NBCP4(L1,L2,L3)
      !     ELSEIF(L4 .EQ. 5)THEN   
      !        N_TEMP = NBCP5(L1,L2,L3) 
      !     ENDIF
      !     
      !     Do N = 1, N_TEMP
      !
      !      IF(L4 .EQ. 1)THEN   
      !        NP_CAND1(I) = NP_CAND1(I) + 1   
      !        IP_CAND1(I,NP_CAND1(I)) = IBCP1(L1,L2,L3,N)
      !      ELSEIF(L4 .EQ. 2)THEN   
      !        NP_CAND2(I) = NP_CAND2(I) + 1   
      !        IP_CAND2(I,NP_CAND2(I)) = IBCP2(L1,L2,L3,N)
      !      ELSEIF(L4 .EQ. 3)THEN   
      !        NP_CAND3(I) = NP_CAND3(I) + 1   
      !        IP_CAND3(I,NP_CAND3(I)) = IBCP3(L1,L2,L3,N)  
      !      ELSEIF(L4 .EQ. 4)THEN   
      !        NP_CAND4(I) = NP_CAND4(I) + 1   
      !        IP_CAND4(I,NP_CAND4(I)) = IBCP4(L1,L2,L3,N)
      !      ELSEIF(L4 .EQ. 5)THEN   
      !        NP_CAND5(I) = NP_CAND5(I) + 1   
      !        IP_CAND5(I,NP_CAND5(I)) = IBCP5(L1,L2,L3,N) 
      !      ENDIF           