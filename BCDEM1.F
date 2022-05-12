c-----------------------------------------------------------------------------------c
c                             SUBROUTINE PROGRAM                                    #
C                             DESINGED BY ChenJun                                   #
C                                                 SHANGHAI JIAO TONG UNIVERSITY     #
C                                                 SHANGHAI, CHINA                   #
c-----------------------------------------------------------------------------------#
c 子程序功能：
c 1------简单侧壁边界和简单床底边界接触虚拟颗粒生成     
c-----------------------------------------------------------------------------------#
      
      Subroutine BCDEM1
      Include './Include/OCERM_INF'
      INTEGER I,J,N_C,K,L,NNNN
      DOUBLE PRECISION A_BC,B_BC,C_BC,ABC_BC,DX_BC,DY_BC,DZ_BC,DIS_BC
      DOUBLE PRECISION S_BC,Z_N,T_BC,X2_BC,Y2_BC,Z2_BC,XN1,YN1,XN2,YN2
      DOUBLE PRECISION XPGDEM,DISTANCE
      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,N_C,K,L,NNNN,A_BC,B_BC,C_BC,
!$OMP&         ABC_BC,DX_BC,DY_BC,DZ_BC,DIS_BC,S_BC,Z_N,T_BC,X2_BC,Y2_BC,
!$OMP&         Z2_BC,XN1,YN1,XN2,YN2,XPGDEM,DISTANCE)     
      
!$OMP DO      
      Do I = 1 , I_PARTICLE

      IF(NSYM_CONT(I) .EQ. 1 .AND. NSYM_MOVE(I) .EQ. 1 )THEN          
          
        NP_BCD(I) = 0         
c--- Left, right and side boundary                                        

        IF(XPDEM(I) .LT. XMINP + RPDEM(I))THEN
                
            NP_BCD(I) = NP_BCD(I) + 1   ! saving contact boundary numbers
            NCONT_INFO(I,NP_BCD(I)) = 1   ! saving contact boundary type (1:LEFT, 2:RIGHT, 3:SIDE1, 4:SIDE2, 5:INNER, 6:BOTTOM, 7:FREE SURFACE)                

            BCDEM_INFO(I,3*NP_BCD(I)-2) = XMINP - RPDEM(I)   !  coordinate of symmetrical particle
            BCDEM_INFO(I,3*NP_BCD(I)-1) = YPDEM(I)
            BCDEM_INFO(I,3*NP_BCD(I)) = ZPDEM(I)  
            
            DEL_DEM(I,NP_BCD(I)) = RPDEM(I) - (XPDEM(I) - XMINP)
            
        ENDIF              
              
        IF(XPDEM(I) .GT. XMAXP - RPDEM(I))THEN
               
            NP_BCD(I) = NP_BCD(I) + 1
            NCONT_INFO(I,NP_BCD(I)) = 2
              
            BCDEM_INFO(I,3*NP_BCD(I)-2) = XMAXP + RPDEM(I)
            BCDEM_INFO(I,3*NP_BCD(I)-1) = YPDEM(I)
            BCDEM_INFO(I,3*NP_BCD(I)) = ZPDEM(I) 
            
            DEL_DEM(I,NP_BCD(I)) = RPDEM(I) - (XMAXP - XPDEM(I))
            
        ENDIF              
              
        IF(YPDEM(I) .LT. YMINP + RPDEM(I))THEN
                
            NP_BCD(I) = NP_BCD(I) + 1   
            NCONT_INFO(I,NP_BCD(I)) = 3   
              
            BCDEM_INFO(I,3*NP_BCD(I)-2) = XPDEM(I)   
            BCDEM_INFO(I,3*NP_BCD(I)-1) = YMINP - RPDEM(I)
            BCDEM_INFO(I,3*NP_BCD(I)) = ZPDEM(I) 
            
            DEL_DEM(I,NP_BCD(I)) = RPDEM(I) - (YPDEM(I) - YMINP)
            
        ENDIF
              
        IF(YPDEM(I) .GT. YMAXP - RPDEM(I))THEN
               
            NP_BCD(I) = NP_BCD(I) + 1
            NCONT_INFO(I,NP_BCD(I)) = 4
              
            BCDEM_INFO(I,3*NP_BCD(I)-2) = XPDEM(I)
            BCDEM_INFO(I,3*NP_BCD(I)-1) = YMAXP + RPDEM(I)
            BCDEM_INFO(I,3*NP_BCD(I)) = ZPDEM(I)    
            
            DEL_DEM(I,NP_BCD(I)) = RPDEM(I) - (YMAXP - YPDEM(I))
            
        ENDIF
        
c--- Cylinder boundary      
        IF(N_CYLINDER .NE. 0)THEN
          DO J = 1, N_CYLINDER  
            DIS_BC = SQRT( (XPDEM(I)-X_CYLINDER(J))**2 + 
     &                (YPDEM(I)-Y_CYLINDER(J))**2 +
     &                (XPDEM(I)-Z_CYLINDER(J))**2 )  
              
            IF( DIS_BC .LT. (RPDEM(I) + R_CYLINDER(J)) ) THEN

               NP_BCD(I) = NP_BCD(I) + 1
               NCONT_INFO(I,NP_BCD(I)) = 5
               
               BCDEM_INFO(I,3*NP_BCD(I)-2) = X_CYLINDER(J)
               BCDEM_INFO(I,3*NP_BCD(I)-1) = Y_CYLINDER(J)
               BCDEM_INFO(I,3*NP_BCD(I)) = ZPDEM(I)

               DEL_DEM(I,NP_BCD(I)) = RPDEM(I) + R_CYLINDER(J) - 
     &                                DIS_BC
               
            ENDIF
          ENDDO
        ENDIF
        
c--- Bottom boundary
        IF(ZPDEM(I) .LT. ZMINP + RPDEM(I))THEN
            
            NP_BCD(I) = NP_BCD(I) + 1
            NCONT_INFO(I,NP_BCD(I)) = 6
               
            BCDEM_INFO(I,3*NP_BCD(I)-2) = XPDEM(I)
            BCDEM_INFO(I,3*NP_BCD(I)-1) = YPDEM(I)
            BCDEM_INFO(I,3*NP_BCD(I)) = ZMINP - RPDEM(I)  
            
            DEL_DEM(I,NP_BCD(I)) = RPDEM(I) - (ZPDEM(I) - ZMINP)
            
        ENDIF

       ENDIF 
        
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      !print*,'NCONT_INFO(1,1)=',NCONT_INFO(1,1)
      Return
      END