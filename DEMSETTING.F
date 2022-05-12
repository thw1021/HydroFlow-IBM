c-----------------------------------------------------------------------------------c
c                             SUBROUTINE PROGRAM                                    #
C                             DESINGED BY ChenJun                                   #
C                                                 SHANGHAI JIAO TONG UNIVERSITY     #
C                                                 SHANGHAI, CHINA                   #
c-----------------------------------------------------------------------------------#
c 子程序功能：
c 1------DEM变量初始化
c 2------颗粒初始信息读入（从文件）   
c 3------DEM计算设置读入（从文件）      
c-----------------------------------------------------------------------------------#
      
      Subroutine DEMSETTING
      Include './Include/OCERM_INF'

        CALL DEMSETREADING
      
        IF(I_TIMESAVE .EQ. 1)THEN
          Open(1121,file='./Results/DEM-CPU-TIME-TOTAL.DAT')
        ENDIF

c        Write (1121,'(9X8A6)') ,' NNDEM',' TOTAL','IPPKPP','IPJPKP',
c     &                 ' FCOUP',' BCDEM','  SOFT',' PMOVE','ACHIVE' 
        IF(I_TIMESAVE .EQ. 1)THEN      
          Open(1131,file='./Results/DEM-CPU-TIME-SELECT.DAT')
        ENDIF
        
c        Write (1131,'(9X8A6)') ,' NNDEM',' TOTAL','IPPKPP','IPJPKP',
c     &                 ' FCOUP',' BCDEM','  SOFT',' PMOVE','ACHIVE'
        
        !PRINT*,'INPUT NREAD:'
        !READ(*,*)NREAD
        
        I_READ_DEMSET = 0
        I_READ_PARTICLE_SOFT = 0
      
        NNDEM = 0        
        TTTDEM = 0.0
        N_TEMP1 = 0
        
        DO I = 1, N_PSLEV
          IF(I .EQ. 1)THEN
              
            NBX(I) = NBX1
            NBY(I) = NBY1
            NBZ(I) = NBZ1
            DELBM(I) = DELBM1
            IBCPMAX(I) = IBCPMAX1
            
C            NBL1 = NBX1 * NBY1 * NBZ1
            NBCP_INDEX(I) = 0
             
          ELSEIF(I .EQ. 2)THEN
              
            NBX(I) = NBX2
            NBY(I) = NBY2
            NBZ(I) = NBZ2
            DELBM(I) = DELBM2
            IBCPMAX(I) = IBCPMAX2
            
C            NBL2 = NBX2 * NBY2 * NBZ2
            NBCP_INDEX(I) = NBL1
            
          ELSEIF(I .EQ. 3)THEN
              
            NBX(I) = NBX3
            NBY(I) = NBY3
            NBZ(I) = NBZ3
            DELBM(I) = DELBM3
            IBCPMAX(I) = IBCPMAX3
            
C            NBL3 = NBX3 * NBY3 * NBZ3
            NBCP_INDEX(I) = NBL1 + NBL2
            
          ELSEIF(I .EQ. 4)THEN
              
            NBX(I) = NBX4
            NBY(I) = NBY4
            NBZ(I) = NBZ4
            DELBM(I) = DELBM4
            IBCPMAX(I) = IBCPMAX4
            
C            NBL4 = NBX4 * NBY4 * NBZ4
            NBCP_INDEX(I) = NBL1 + NBL2 + NBL3
            
          ELSEIF(I .EQ. 5)THEN
              
            NBX(I) = NBX5
            NBY(I) = NBY5
            NBZ(I) = NBZ5
            DELBM(I) = DELBM5
            IBCPMAX(I) = IBCPMAX5
            
C            NBL5 = NBX5 * NBY5 * NBZ5
            NBCP_INDEX(I) = NBL1 + NBL2 + NBL3 + NBL4
            
          ELSE
              
            PRINT*,'ERROR OF N_PSLEV IN OCERM_INF'  
            PAUSE
            
          ENDIF
          
          DO L = 1, NBZ(I)
                
           DO K = 1, NBY(I)
                  
            DO J = 1, NBX(I) 
         
             N_TEMP = NBCP_INDEX(I) + (L-1)*NBX(I)*NBY(I) + 
     &                NBX(I)*(K-1) + J   !根据层级和所在层I,J,K确定网格全局编号          
               
             IBCP_INDEX(N_TEMP) = N_TEMP1   !确定当前网格颗粒存储起始位置

             N_TEMP1 = N_TEMP1 + IBCPMAX(I)   !当前网格颗粒存储终止位置，作为下一网格颗粒存储起始位置
             
            ENDDO
              
           ENDDO
             
          ENDDO          
          
        ENDDO    
        
        !NBL = NBL1 + NBL2 + NBL3 + NBL4 + NBL5
          
        CALL PARTICLE_INFO

C--- 床底单元边方向向量
        DO I = 1, IJE
          
          N1 = IEND_EDGE(I,1)
          N2 = IEND_EDGE(I,2)
            
          VECTOR_BE(I,1) = PXY(N2,1) - PXY(N1,1) 
          VECTOR_BE(I,2) = PXY(N2,2) - PXY(N1,2) 
          VECTOR_BE(I,3) = - HP(N2) + HP(N1)
          VECTOR_BE(I,4) = VECTOR_BE(I,1)**2 + VECTOR_BE(I,2)**2 +
     &                     VECTOR_BE(I,3)**2
          
        ENDDO
          
C        DO I = 1, IJM
            
C         DO J = 1, CELL_POLYGEN(I)   
          
C          N1 = IEND_EDGE(I,1)
C          N2 = IEND_EDGE(I,2)
            
C          VECTOR_BE(I,1) = PXY(N2,1) - PXY(N1,1) 
C          VECTOR_BE(I,2) = PXY(N2,2) - PXY(N1,2) 
C          VECTOR_BE(I,3) = - HP(N2) + HP(N1)
C          VECTOR_BE(I,4) = VECTOR_BE(I,1)**2 + VECTOR_BE(I,2)**2 +
C     &                     VECTOR_BE(I,3)**2
          
C        ENDDO   
          
C---  床底非与周围网格不共面的单元搜索 
      DO I = 1, IJM
        
         NSYM_BC(I) = 0 
          
         IF(CELL_POLYGEN(I) .EQ. 3)THEN
           NLS = 3
         ELSEIF(CELL_POLYGEN(I) .EQ. 4)THEN
           NLS = 2
         ELSE
           PRINT*,'ERROR OF CELL_POLYGEN(I) IN IPPKPPINITIAL.F' 
           PAUSE
         ENDIF
          
         DO J = 2, INL(I,NLS) 
             
          K =  INE(I,J)           
          
          TEMP1 = SQRT(DEMCOSB(I,1)**2+DEMCOSB(I,2)**2+DEMCOSB(I,3)**2)
          TEMP2 = SQRT(DEMCOSB(K,1)**2+DEMCOSB(K,2)**2+DEMCOSB(K,3)**2)
          TEMP3 = DEMCOSB(I,1)*DEMCOSB(K,1)+DEMCOSB(I,2)*DEMCOSB(K,2)+
     &            DEMCOSB(I,3)*DEMCOSB(K,3)
          TEMP4 = TEMP3 / (TEMP1 * TEMP2)

          IF(ABS(ABS(TEMP4) - 1.0) .GT. 1.E-10)THEN
          
            NSYM_BC(I) = 1
          
          ENDIF

         
         ENDDO
         
       ENDDO

      !DO I = 1, IJM
      !    IF(NSYM_BC(I) == 1)PRINT*,I
      !ENDDO
      !PRINT*,CXY(1198,1),CXY(1198,2)                              
      
      End  