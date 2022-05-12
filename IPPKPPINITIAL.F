c-----------------------------------------------------------------------------------c
c                             SUBROUTINE PROGRAM                                    #
C                             DESINGED BY ChenJun                                   #
C                                                 SHANGHAI JIAO TONG UNIVERSITY     #
C                                                 SHANGHAI, CHINA                   #
c-----------------------------------------------------------------------------------#
c 子程序功能：
c 1------初始计算颗粒流场网格寻址及颗粒存储      
c 1------初始计算背景网格寻址及颗粒存储       
c-----------------------------------------------------------------------------------#
      
      Subroutine IPPKPPINITIAL
      Include './Include/OCERM_INF'
      Parameter(NSCHEME=2)
      Double Precision X1DEM,X2DEM,Y1DEM,Y2DEM,DIJDEM,XPGDEM
      Double Precision R1,RTOL,ELFPDEM,DCPDEM,HPDEM
      
      CALL NORMAL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,
!$OMP&         X1DEM,X2DEM,Y1DEM,Y2DEM,DIJDEM,XPGDEM,
!$OMP&         R1,RTOL,ELFPDEM,DCPDEM,HPDEM)
!$OMP DO     
      Do I = 1 , I_PARTICLE        

c----horizontal gridding number        
        Do J = 1 , IJM         
            
          Do K = 1 , Cell_POLYGEN(J)       
            X1DEM = PXY(CELL_SIDE(J,K,3),1)
            Y1DEM = PXY(CELL_SIDE(J,K,3),2)
            X2DEM = PXY(CELL_SIDE(J,K,4),1)
            Y2DEM = PXY(CELL_SIDE(J,K,4),2)          
                    
            XPGDEM = ( XPDEM(I) - X1DEM ) * ( Y2DEM - Y1DEM ) -
     &               ( YPDEM(I) - Y1DEM ) * ( X2DEM - X1DEM )
                         
            If(XPGDEM .GT. 0.) Then
              Exit
            Endif    
          Enddo
                
          If(K .EQ. ( Cell_POLYGEN(J) + 1 ) ) Then     
            IPP(I) = J
            Exit
          Endif  

        Enddo
 
        If( J .EQ. IJM + 1 ) Then
            Print*, 'Error1 of IPP in IPPKPPINITIAL.F!'
            Print*, I,IPP(I),XPDEM(I),YPDEM(I),ZPDEM(I)           
            Pause
        Endif        

c----颗粒流场网格定位-垂向 
        IF(NSCHEME .EQ. 1)THEN
        !插值计算颗粒水平投影点水位、静水深和总水深
	  RTOL = 0.0
	  ELFPDEM = 0.0
        HPDEM = 0.0
        DCPDEM = 0.0
        
	  Do J = 1, CELL_POLYGEN(IPP(I))              
c	    R1 = Sqrt((XPDEM(I) - PXY( CELL_INDEX(IPP(I),J) , 1)) ** 2. +
c     &	          (YPDEM(I) - PXY( CELL_INDEX(IPP(I),J) , 2)) ** 2.)   !距离作为权重
	    R1 = (XPDEM(I) - PXY( CELL_INDEX(IPP(I),J) , 1)) ** 2. +
     &	     (YPDEM(I) - PXY( CELL_INDEX(IPP(I),J) , 2)) ** 2.   !距离平方作为权重  
          
          RTOL = RTOL + 1. / R1
          
          ELFPDEM = ELFPDEM + ELFV( CELL_INDEX(IPP(I),J) ) / R1
          
c          HPDEM = HPDEM + HP(CELL_INDEX(IPP(I),J)) / R1 
          
          DCPDEM = DCPDEM + ( HP(CELL_INDEX(IPP(I),J)) + 
     &             ELFV(CELL_INDEX(IPP(I),J)) ) / R1 
	  Enddo
	  ELFPDEM = ELFPDEM / RTOL
        HPDEM = HPDEM / RTOL
        DCPDEM = DCPDEM / RTOL    
      
        !SIGMA_P(I) = ( ZPDEM(I) - ELF(IPP(I)) ) / DC(IPP(I))
        SIGMA_P(I) = ( ZPDEM(I) - ELFPDEM ) / DCPDEM
        !SIGMA_P(I) = ( ZPDEM(I) - ELF(IPP(I)) ) / (ELF(IPP(I))+HPDEM) 
        
        ELSEIF(NSCHEME .EQ. 2)THEN

         HPDEM = ( DEMCOSB(IPP(I),1) * (XPDEM(I) - 
     &             PXY( CELL_INDEX(IPP(I),1) , 1)) +
     &             DEMCOSB(IPP(I),2) * (YPDEM(I) - 
     &             PXY( CELL_INDEX(IPP(I),1) , 2)) ) /
     &             DEMCOSB(IPP(I),3) + HP(CELL_INDEX(IPP(I),1))    
         
         ELFPDEM = - ( DEMCOSA(IPP(I),1) * (XPDEM(I) - 
     &               PXY( CELL_INDEX(IPP(I),1) , 1)) +
     &               DEMCOSA(IPP(I),2) * (YPDEM(I) - 
     &               PXY( CELL_INDEX(IPP(I),1) , 2)) ) /
     &               DEMCOSA(IPP(I),3) + ELFV(CELL_INDEX(IPP(I),1))  
         
         SIGMA_P(I) = ( ZPDEM(I) - ELFPDEM ) / ( HPDEM + ELFPDEM ) 
         
        ELSEIF(NSCHEME .EQ. 3)THEN       

         HPDEM = 0.5 * (
     &           ( DEMCOSB(IPP(I),1) * (XPDEM(I) - 
     &             PXY( CELL_INDEX(IPP(I),1) , 1)) +
     &             DEMCOSB(IPP(I),2) * (YPDEM(I) - 
     &             PXY( CELL_INDEX(IPP(I),1) , 2)) ) /
     &             DEMCOSB(IPP(I),3) + HP(CELL_INDEX(IPP(I),1))  
     &           + 
     &           ( DEMCOSB(IPP(I),4) * (XPDEM(I) - 
     &             PXY( CELL_INDEX(IPP(I),1) , 1)) +
     &             DEMCOSB(IPP(I),5) * (YPDEM(I) - 
     &             PXY( CELL_INDEX(IPP(I),1) , 2)) ) /
     &             DEMCOSB(IPP(I),6) + HP(CELL_INDEX(IPP(I),1)) 
     &           )
         
         ELFPDEM = 0.5 * ( ELFV(CELL_INDEX(IPP(I),1)) -
     &             ( DEMCOSA(IPP(I),1) * (XPDEM(I) - 
     &               PXY( CELL_INDEX(IPP(I),1) , 1)) +
     &               DEMCOSA(IPP(I),2) * (YPDEM(I) - 
     &               PXY( CELL_INDEX(IPP(I),1) , 2)) ) /
     &               DEMCOSA(IPP(I),3)    
     &            +  ELFV(CELL_INDEX(IPP(I),1)) -  
     &             ( DEMCOSA(IPP(I),4) * (XPDEM(I) - 
     &               PXY( CELL_INDEX(IPP(I),1) , 1)) +
     &               DEMCOSA(IPP(I),5) * (YPDEM(I) - 
     &               PXY( CELL_INDEX(IPP(I),1) , 2)) ) /
     &               DEMCOSA(IPP(I),6)  
     &             ) 
         
         SIGMA_P(I) = ( ZPDEM(I) - ELFPDEM ) / ( HPDEM + ELFPDEM )
         
        ENDIF
        
        IF(SIGMA_P(I) .LT. -1.0)THEN
          PRINT*,'Error2 of KPP in IPPKPPINITIAL.F!'
          PRINT*,'SIGMA_P(I) = ',SIGMA_P(I)
          PAUSE
        ELSEIF(SIGMA_P(I) .GT. 0.0)THEN  
          KPP(I) = 0
        ELSE            
          Do K = 1 , KBM
            If ( SIGMA_P(I) .LE. Z(K) .AND. 
     &        SIGMA_P(I) .GE. Z(K+1) ) Then
              KPP(I) = K           
              exit 
            Endif
          Enddo
        ENDIF 

C- 背景网格寻址       
        IF(PSMOD .EQ. 'MOD2')THEN

         IF(NSYM_CONT(I) .EQ. 1 .AND. I_PSLEV(I) .NE. 0)THEN          
            
         DO K = I_PSLEV(I), N_PSLEV  
            
          DO J = 1, NBX(K) 
              
           IF(XPDEM(I) .GE. XMV(J,K) .AND. XPDEM(I) .LT. XMV(J+1,K))THEN
              IP(I,K) = J
              EXIT
           ENDIF
           
          ENDDO
          
          IF(J .EQ. NBX(K)+1) THEN
            Write (IUPRT,*) 'Particle I= ',I,' is vanished in 
     &                       IPPKPPINITIAL.F!'    
            Write (IUPRT,*) 'XPDEM(I)=',XPDEM(I)
            Write (IUPRT,*) 'YPDEM(I)=',YPDEM(I)
            Write (IUPRT,*) 'ZPDEM(I)=',ZPDEM(I)
            NSYM_CONT(I) = 0
          ENDIF
          
          DO J = 1, NBY(K)
              
           IF(YPDEM(I) .GE. YMV(J,K) .AND. YPDEM(I) .LT. YMV(J+1,K))THEN
              JP(I,K) = J
              EXIT
           ENDIF
           
          ENDDO
          
          IF(J .EQ. NBY(K)+1) THEN
            Write (IUPRT,*) 'Particle I= ',I,' is vanished in 
     &                       IPPKPPINITIAL.F!'    
            Write (IUPRT,*) 'XPDEM(I)=',XPDEM(I)
            Write (IUPRT,*) 'YPDEM(I)=',YPDEM(I)
            Write (IUPRT,*) 'ZPDEM(I)=',ZPDEM(I)
            NSYM_CONT(I) = 0
          ENDIF
          
          DO J = 1, NBZ(K)
              
           IF(ZPDEM(I) .GE. ZMV(J,K) .AND. ZPDEM(I) .LT. ZMV(J+1,K))THEN
              KP(I,K) = J
              EXIT
           ENDIF
           
          ENDDO
          
          IF(J .EQ. NBZ(K)+1) THEN
            Write (IUPRT,*) 'Particle I= ',I,' is vanished in 
     &                       IPPKPPINITIAL.F!'    
            Write (IUPRT,*) 'XPDEM(I)=',XPDEM(I)
            Write (IUPRT,*) 'YPDEM(I)=',YPDEM(I)
            Write (IUPRT,*) 'ZPDEM(I)=',ZPDEM(I)
            NSYM_CONT(I) = 0
          ENDIF
         
         ENDDO 
          
         ENDIF
                
        ENDIF

        IIP(I) = IP(I,I_PSLEV(I)) 
        JJP(I) = JP(I,I_PSLEV(I)) 
        KKP(I) = KP(I,I_PSLEV(I))         
        
      Enddo
!$OMP END DO
!$OMP END PARALLEL 

      DO I = 1, NBL    
        NBCP( I ) = 0
      ENDDO 

      Do I = 1 , I_PARTICLE 
          
       IF(NSYM_CONT(I) .EQ. 1 .AND. I_PSLEV(I) .NE. 0)THEN    
         
         N = I_PSLEV(I)    
           
         N_TEMP = NBCP_INDEX(N) + (KP(I,N)-1)*NBX(N)*NBY(N) + 
     &            NBX(N)*(JP(I,N)-1) + IP(I,N)   !颗粒在N层级所属网格单元在NBCP数组中的存储位置
          
         NBCP( N_TEMP ) = NBCP( N_TEMP ) + 1   !颗粒在N层级所属网格单元颗粒数+1    
         
         IF(NBCP( N_TEMP ) .GT. IBCPMAX(N))THEN
            PRINT*,'Error4 of IBCPMAX',N,' in IPJPKPDEM.F!'
            WRITE(IUPRT,*)
            WRITE(IUPRT,'(3A14)') 'XPDEM', 'YPDEM', 'ZPDEM' 
            WRITE(IUPRT,'(3F14.7)') XPDEM(I),YPDEM(I),ZPDEM(I)
            PAUSE
         ENDIF          
         
         !IBCP( (N_TEMP - 1) * IBCPMAX + NBCP( N_TEMP ) ) = I   !将颗粒I单元编号保存到颗粒在所在层级所属网格单元
         IBCP( IBCP_INDEX(N_TEMP) + NBCP( N_TEMP ) ) = I   !将颗粒I单元编号保存到颗粒在所在层级所属网格单元         
           
       ENDIF
       
      ENDDO
      
      Return 
      End
