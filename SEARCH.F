c-----------------------------------------------------------------------------------c
c                             SUBROUTINE PROGRAM                                    #
C                             DESINGED BY ChenJun                                   #
C                                                 SHANGHAI JIAO TONG UNIVERSITY     #
C                                                 SHANGHAI, CHINA                   #
c-----------------------------------------------------------------------------------#
c 子程序功能：
c 1------计算简单边界侧壁法向量和床底坐标，识别并记录侧壁边界（左右上下）及内部边界各自的边数及边编号       
c 2------搜索水平网格邻单元编号
c 3------生成背景网格坐标                                  
c-----------------------------------------------------------------------------------c
      Subroutine SEARCH
      Include './Include/OCERM_INF'
      INTEGER N_SYM(IJM)
      !INE(IJM,(NL-1)*NL*2+1)   !网格I(第一层)和其周围网格(第N层)编号数组，(NL-1)*NL*2+1为每个网格周围N层最大网格数
      !INL(IJM,NL)：网格IJM周围n层实际网格数
      
C---计算简单边界侧壁法向量和床底坐标，识别并记录侧壁边界（左右上下）及内部边界各自的边数及边编号
      N_LEFT = 0
      N_RIGHT = 0
      N_SIDE1 = 0
      N_SIDE2 = 0
      N_INNER = 0
      
      DO K = 1, IJE
          DO I = 1, 4
              NE_INFO(I,K) = 0
          ENDDO
      ENDDO

      XMAXP = PXY(1,1)
      YMAXP = PXY(1,2)
      XMINP = PXY(1,1)
      YMINP = PXY(1,2)      
      
      DO J = 2, IJP
          XMAXP = MAX( XMAXP , PXY(J,1) )
          YMAXP = MAX( YMAXP , PXY(J,2) )
          XMINP = MIN( XMINP , PXY(J,1) )
          YMINP = MIN( YMINP , PXY(J,2) )   
      ENDDO
      
      ZMINP = - MAXVAL(HP)
      DO I = 1, IJM
          
        IF(ABS(ZMINP)-ABS(HP(I)) .GT. 1.E-5 .AND. 
     &     BCDMOD .EQ. 'MOD2')THEN
          PRINT*,'Checking the BCDMOD in DEM_SETTING.DAT'  
          PAUSE
        ENDIF

      ENDDO 
      
      DO K = 1, IJE
          If(CFM(K) .EQ. -1.0 .OR. CFM(K) .EQ. -2.0 
     &        .OR. CFM(K) .EQ. -3.0) Then
c---  normal vector              
              DX_EDGE = PXY(IEND_EDGE(K,1),1) - PXY(IEND_EDGE(K,2),1)
              DY_EDGE = PXY(IEND_EDGE(K,1),2) - PXY(IEND_EDGE(K,2),2)
              BE_NORMAL(K,1) =   DY_EDGE / SQRT(DX_EDGE**2 + DY_EDGE**2)
              BE_NORMAL(K,2) = - DX_EDGE / SQRT(DX_EDGE**2 + DY_EDGE**2)
              BE_NORMAL(K,3) =   0.0
              
c---  boundary edges              
              X_POINT1 = PXY(IEND_EDGE(K,1),1)
              Y_POINT1 = PXY(IEND_EDGE(K,1),2)
              X_POINT2 = PXY(IEND_EDGE(K,2),1)
              Y_POINT2 = PXY(IEND_EDGE(K,2),2)
              
              IF( X_POINT1 .LT. (XMINP + 0.0001) .AND.
     &            X_POINT2 .LT. (XMINP + 0.0001) ) THEN
                  N_LEFT = N_LEFT + 1                  
                  NE_INFO(1,N_LEFT) = K              
              ELSEIF( X_POINT1 .GT. (XMAXP - 0.0001) .AND.
     &            X_POINT2 .GT. (XMAXP - 0.0001) ) THEN
                  N_RIGHT = N_RIGHT + 1                  
                  NE_INFO(2,N_RIGHT) = K
              ELSEIF( Y_POINT1 .LT. (YMINP + 0.0001) .AND.
     &            Y_POINT2 .LT. (YMINP + 0.0001) ) THEN
                  N_SIDE1 = N_SIDE1 + 1                  
                  NE_INFO(3,N_SIDE1) = K
              ELSEIF( Y_POINT1 .GT. (YMAXP - 0.0001) .AND.
     &            Y_POINT2 .GT. (YMAXP - 0.0001) ) THEN
                  N_SIDE2 = N_SIDE2 + 1                  
                  NE_INFO(4,N_SIDE2) = K                  
              ELSE
                  N_INNER = N_INNER + 1                  
                  NE_INFO(5,N_INNER) = K
              ENDIF              
          ENDIF
      ENDDO 
      NE_BCD(1) = N_LEFT
      NE_BCD(2) = N_RIGHT
      NE_BCD(3) = N_SIDE1
      NE_BCD(4) = N_SIDE2
      NE_BCD(5) = N_INNER            

C---搜索任一网格周围NL层水平网格，排序并记录其编号      
      DO I = 1, IJM
        INE(I,1) = I

        DO J = 1, IJM
          N_SYM(J) = 0
        ENDDO
        
        N_SYM(I) = 1
        NC1 = 1
        INL(I,1) = 1
        
        DO J = 2, NL
          
          !NC2 = 0  
            
          IF(J .EQ. 2)THEN
            IS = 1
            IE = 1
          ELSE
            IS = INL(I,J-2) + 1
            IE = INL(I,J-1)           
          ENDIF

          DO K = IS, IE 
            !PRINT*,K,INE(I,K)  
            DO L = 1, CELL_POLYGEN(INE(I,K)) 
              IC = CELL_SIDE(INE(I,K), L, 2)
              IF(IC .NE. -999)THEN
               IF(N_SYM(IC) .EQ. 0)THEN   
                NC1 = NC1 + 1  
                INE(I,NC1) = IC 
                !print*,NC1,INE(I,NC1)
                N_SYM(IC) = 1
               ENDIF
              ENDIF
            ENDDO                                    
          ENDDO
          INL(I,J) = NC1
          !print*,I,J,INL(I,J)
        ENDDO        
        !IF(I .EQ. 1)PRINT*,INL(I,NL),(INE(1,K),K=1,INL(I,NL))
      ENDDO
      !PAUSE

C---背景网格生成      
      IF(PSMOD .EQ. 'MOD2')THEN
          
       DO J = 1, N_PSLEV   
      
        DO I = 1, NBX(J)+1
          XMV(I,J) = XBMIN + DFLOAT(I-1) * DELBM(J) * 1.0001
        ENDDO

        DO I = 1, NBY(J)+1  
          YMV(I,J) = YBMIN + DFLOAT(I-1) * DELBM(J) *  1.0001 
        ENDDO
      
        DO I = 1, NBZ(J)+1  
          ZMV(I,J) = ZBMIN + DFLOAT(I-1) * DELBM(J) * 1.0001 
        ENDDO
        
       ENDDO 
        
      ENDIF
      
      Return 
      End
C---------------------------------------------------------