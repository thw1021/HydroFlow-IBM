c-----------------------------------------------------------------------------------c
c                             SUBROUTINE PROGRAM                                    #
C                             DESINGED BY ChenJun                                   #
C                                                 SHANGHAI JIAO TONG UNIVERSITY     #
C                                                 SHANGHAI, CHINA                   #
c-----------------------------------------------------------------------------------#
c 子程序功能：
c 1------简单侧壁边界和复杂床底边界接触虚拟颗粒生成     
c-----------------------------------------------------------------------------------#
      
      Subroutine BCDEM
      Include './Include/OCERM_INF'
      PARAMETER(N_SCHEME = 1)
      INTEGER I,J,N_C,K,L,NNNN,N_TEMP1
      DOUBLE PRECISION A_BC,B_BC,C_BC,ABC_BC,DX_BC,DY_BC,DZ_BC,DIS_BC
      DOUBLE PRECISION S_BC,Z_N,T_BC,X2_BC,Y2_BC,Z2_BC
      DOUBLE PRECISION XPGDEM,VX,VY,VZ
      DOUBLE PRECISION DISTANCE,DISX,DISY,DISZ
      DOUBLE PRECISION CONST_A,CONST_B,CONST_C,DELTA
      DOUBLE PRECISION T_INT1,T_INT2,DEL_AVE
      DOUBLE PRECISION XN1(10),YN1(10),XN2(10),YN2(10)
      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,N_C,K,L,NNNN,N_TEMP1,VX,VY,
!$OMP&         VZ,A_BC,B_BC,C_BC,ABC_BC,DX_BC,DY_BC,DZ_BC,DIS_BC,S_BC,
!$OMP&         Z_N,T_BC,X2_BC,Y2_BC,Z2_BC,XN1,YN1,XN2,YN2,XPGDEM,
!$OMP&         DISTANCE,DISX,DISY,DISZ,CONST_A,CONST_B,CONST_C,DELTA,
!$OMP&         T_INT1,T_INT2,DEL_AVE)    
      
!$OMP DO      
      Do I = 1 , I_PARTICLE
          
       IF( NSYM_CONT(I) .EQ. 1 .AND. NSYM_MOVE(I) .EQ. 1 )THEN    
          
        NP_BCD(I) = 0         
c--- Left, right and side boundary        
        DO J = 1, 4
            IF(NE_BCD(J) .NE. 0)THEN
                A_BC = BE_NORMAL(NE_INFO(J,1),1)   ! normal vector
                B_BC = BE_NORMAL(NE_INFO(J,1),2)
                C_BC = BE_NORMAL(NE_INFO(J,1),3)
                ABC_BC = A_BC**2+B_BC**2+C_BC**2
                  
                DX_BC = PXY(IEND_EDGE(NE_INFO(J,1),1),1)-XPDEM(I)   ! vector from any point of boundary face to particle center
                DY_BC = PXY(IEND_EDGE(NE_INFO(J,1),1),2)-YPDEM(I)
                DZ_BC = 0.0
                  
                DIS_BC = ABS(A_BC*DX_BC+B_BC*DY_BC+C_BC*DZ_BC)/
     &                              SQRT(ABC_BC)   ! distance between paticle center and boundary face                                 
                
                IF(DIS_BC .LT. RPDEM(I))THEN
                    NP_BCD(I) = NP_BCD(I) + 1   ! saving contact boundary numbers
                    NCONT_INFO(I,NP_BCD(I)) = J   ! saving contact boundary type (1:LEFT, 2:RIGHT, 3:SIDE1, 4:SIDE2, 5:INNER, 6:BOTTOM, 7:FREE SURFACE)
                      
 
                    T_BC = (A_BC*DX_BC+B_BC*DY_BC+C_BC*DZ_BC)/ABC_BC
                        
                    X2_BC = XPDEM(I)+A_BC*T_BC   ! intersection point coordinate
                    Y2_BC = YPDEM(I)+B_BC*T_BC
                    Z2_BC = ZPDEM(I)+C_BC*T_BC
                        
                    BCDEM_INFO(I,3*NP_BCD(I)-2) = X2_BC+RPDEM(I)/
     &                                    DIS_BC*(X2_BC-XPDEM(I))   !  coordinate of symmetrical particle
                    BCDEM_INFO(I,3*NP_BCD(I)-1) = Y2_BC+RPDEM(I)/
     &                                    DIS_BC*(Y2_BC-YPDEM(I))
                    BCDEM_INFO(I,3*NP_BCD(I)) = Z2_BC+RPDEM(I)/
     &                                    DIS_BC*(Z2_BC-ZPDEM(I))  
                    
                    DEL_DEM(I,NP_BCD(I)) = RPDEM(I) - DIS_BC
                    
                ENDIF
            ENDIF            
        ENDDO
        
c--- Cylinder boundary      
        IF(N_CYLINDER .NE. 0)THEN
          DO J = 1, N_CYLINDER  
              
            DIS_BC = SQRT( (XPDEM(I)-X_CYLINDER(J))**2 + 
     &                     (YPDEM(I)-Y_CYLINDER(J))**2 +
     &                     (XPDEM(I)-Z_CYLINDER(J))**2 )              
              
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
        N_TEMP1 = 0
        
        DO J = 1, IJM
          A_BC = DEMCOSB(J,1)
          B_BC = DEMCOSB(J,2)
          C_BC = DEMCOSB(J,3)

          ABC_BC = A_BC**2+B_BC**2+C_BC**2
          
          N_C = CELL_INDEX(J,1)
          Z_N = Z(KB) * (HP(N_C) + ELFV(N_C)) + ELFV(N_C)
          
          DX_BC = PXY(N_C,1)-XPDEM(I)
          DY_BC = PXY(N_C,2)-YPDEM(I)                              
          DZ_BC = Z_N - ZPDEM(I)
                  
          DIS_BC = ABS(A_BC*DX_BC+B_BC*DY_BC+C_BC*DZ_BC)/
     &                              SQRT(ABC_BC) 
          
          IF(DIS_BC .LT. RPDEM(I))THEN
              
            T_BC = (A_BC*DX_BC+B_BC*DY_BC+C_BC*DZ_BC)/ABC_BC   ! intermidiate variable
                      
            X2_BC = XPDEM(I)+A_BC*T_BC      ! 交点坐标
            Y2_BC = YPDEM(I)+B_BC*T_BC
            Z2_BC = ZPDEM(I)+C_BC*T_BC
                      
            DO K =  1, CELL_POLYGEN(J)   ! 判断交点是否位于单元内    
                XN1(K) = PXY(CELL_SIDE(J,K,3),1)
                YN1(K) = PXY(CELL_SIDE(J,K,3),2)
                XN2(K) = PXY(CELL_SIDE(J,K,4),1)
                YN2(K) = PXY(CELL_SIDE(J,K,4),2)           
      
                XPGDEM = ( X2_BC - XN1(K) ) * ( YN2(K) - YN1(K) ) - 
     &                   ( Y2_BC - YN1(K) ) * ( XN2(K) - XN1(K) )
                
                If(XPGDEM .GT. 0.) Then
                    Exit
                Endif
            ENDDO
            
c            DO L = 1, CELL_POLYGEN(J)   ! judging whether the diatance between cell vertex and particle center is smaller than particel radius
c                N_C = CELL_INDEX(J,L)
c                Z_N = Z(KB) * (HP(N_C) + ELFV(N_C)) + ELFV(N_C)
c                DISTANCE=SQRT( (XPDEM(I)-PXY(N_C,1))**2 + 
c     &                         (YPDEM(I)-PXY(N_C,2))**2 +
c     &                         (ZPDEM(I)- Z_N)**2 )
              
c                IF(DISTANCE .LT. RPDEM(I))THEN            
c                    Exit
c                Endif             
c            ENDDO
            
c            If(K .EQ. ( CELL_POLYGEN(J) + 1 ) .OR.
c     &         L .NE. ( CELL_POLYGEN(J) + 1 )) Then             

             If(K .EQ. ( CELL_POLYGEN(J) + 1 )) Then            
            
                NP_BCD(I) = NP_BCD(I) + 1
                NCONT_INFO(I,NP_BCD(I)) = 6
                N_TEMP1 = N_TEMP1 + 1
                          
                BCDEM_INFO(I,3*NP_BCD(I)-2) = X2_BC+RPDEM(I)/
     &                                DIS_BC*(X2_BC-XPDEM(I))
                BCDEM_INFO(I,3*NP_BCD(I)-1) = Y2_BC+RPDEM(I)/
     &                                DIS_BC*(Y2_BC-YPDEM(I))
                BCDEM_INFO(I,3*NP_BCD(I)) = Z2_BC+RPDEM(I)/
     &                                DIS_BC*(Z2_BC-ZPDEM(I))
                
                DEL_DEM(I,NP_BCD(I)) = RPDEM(I) - DIS_BC
                
                !Exit                
            ELSE

             IF(N_SCHEME .EQ. 1 .AND. NSYM_BC(J) .EQ. 1)THEN
                 
              DO K =  1, CELL_POLYGEN(J)  
                  
               VX = VECTOR_BE(CELL_SIDE(J,K,1),1)
               VY = VECTOR_BE(CELL_SIDE(J,K,1),2)
               VZ = VECTOR_BE(CELL_SIDE(J,K,1),3)
 
               !XN1(K) = PXY(CELL_SIDE(J,K,3),1)
               !YN1(K) = PXY(CELL_SIDE(J,K,3),2) 
               XN1(K) = PXY(IEND_EDGE(CELL_SIDE(J,K,1),1),1)
               YN1(K) = PXY(IEND_EDGE(CELL_SIDE(J,K,1),1),2)
               
               DISX = XN1(K) - XPDEM(I)
               DISY = YN1(K) - YPDEM(I)
               DISZ = - HP( IEND_EDGE(CELL_SIDE(J,K,1),1) ) - ZPDEM(I)
               !print*,VX,VY,VZ
               !print*,DISX,DISY,DISZ
               CONST_A = VECTOR_BE(CELL_SIDE(J,K,1),4)
               
               CONST_B = 2. * ( VX * DISX + VY * DISY + VZ * DISZ )
               
               CONST_C = DISX**2 + DISY**2 + DISZ**2 - RPDEM(I)**2
               
               DELTA = CONST_B**2 - 4. * CONST_A * CONST_C            
               
               IF(DELTA .GE. 0.)THEN
                   
                 T_INT1 = ( - CONST_B + SQRT(DELTA) ) / ( 2. * CONST_A )
                 T_INT2 = ( - CONST_B - SQRT(DELTA) ) / ( 2. * CONST_A )

               ENDIF
               
              ENDDO
              
               
              IF(K .NE. ( CELL_POLYGEN(J) + 1 )) Then            
            
                NP_BCD(I) = NP_BCD(I) + 1
                NCONT_INFO(I,NP_BCD(I)) = 6
                !PRINT*,'NP_BCD(I)-222=',NP_BCD(I)
                N_TEMP1 = N_TEMP1 + 1
                          
                BCDEM_INFO(I,3*NP_BCD(I)-2) = X2_BC+RPDEM(I)/
     &                                DIS_BC*(X2_BC-XPDEM(I))
                BCDEM_INFO(I,3*NP_BCD(I)-1) = Y2_BC+RPDEM(I)/
     &                                DIS_BC*(Y2_BC-YPDEM(I))
                BCDEM_INFO(I,3*NP_BCD(I)) = Z2_BC+RPDEM(I)/
     &                                DIS_BC*(Z2_BC-ZPDEM(I))
                !PRINT*,'222',BCDEM_INFO
                !PRINT*,I,J
                DEL_DEM(I,NP_BCD(I)) = DMIN1( RPDEM(I) - DIS_BC , 
     &                                        0.02 * RPDEM(I) )
c                print*,'P1:',XPDEM(I),YPDEM(I),ZPDEM(I)
                !PRINT*,J
                !print*,'P2:',BCDEM_INFO(I,3*NP_BCD(I)-2),
c     &          BCDEM_INFO(I,3*NP_BCD(I)-1),
c     &          BCDEM_INFO(I,3*NP_BCD(I))
c                print*,'DEL_DEM:',DEL_DEM(I,NP_BCD(I))
                !pause
C                IF(DEL_DEM(NP_BCD(I)) .LT. 0.04 * RPDEM(I))THEN
                    
C                 N_TEMP2 = N_TEMP2 + 1 
            
C                 DELTOT_DEM = DELTOT_DEM + DEL_DEM(NP_BCD(I))
                
C                ENDIF
                
              ENDIF
                                          
             ENDIF            
                             
            ENDIF
             
          ENDIF
          
       ENDDO

c--- 多虚拟颗粒重叠量校正
        IF(N_TEMP1 .EQ. 0)THEN   !与床底不接触(大多数情况)
          
         DEL_HIS(I) = 0.   !历史重叠量   
         
        ELSEIF(N_TEMP1 .EQ. 1)THEN   !与床底接触且只生产一个虚拟颗粒(部分情况) 
            
         DEL_HIS(I) = DEL_DEM(I,NP_BCD(I))
         
        ELSE   !与床底接触且生成多个虚拟颗粒(少数情况)          
          
c         DEL_AVE = SUM( ( DEL_DEM(I,K) , K = NP_BCD(I)-N_TEMP1+1,
c     &                NP_BCD(I) ) ) / N_TEMP1   !多个虚拟颗粒重叠量平均值 
         DEL_AVE = 0. 
         DO K = NP_BCD(I)-N_TEMP1+1, NP_BCD(I)
           !PRINT*,DEL_DEM(I,K)  
           DEL_AVE = DEL_AVE + DEL_DEM(I,K) 
           
         ENDDO
         
         DEL_AVE = DEL_AVE / N_TEMP1   !多个虚拟颗粒重叠量平均值          
          
         !IF(DEL_AVE .GT. 1.2 * DEL_HIS(I))THEN   !若超过1.2倍历史重叠量，取平均值作为重叠量
         !  
         ! DEL_AVE = 0.5 * ( DEL_AVE + DEL_HIS(I) )   
         !  
         !ENDIF
         
         DEL_HIS(I) = DEL_AVE
          
         DO J = NP_BCD(I)-N_TEMP1+1 , NP_BCD(I)
              
          DEL_DEM(I,J) = DEL_AVE / N_TEMP1   !将重叠量平均分配到各虚拟颗粒
            
         ENDDO           
         
        ENDIF         
       
c--- Free surface
        NNNN = 0 
        IF(NNNN .EQ. 1)THEN  
         DO J = 1, IJM
          A_BC = DEMCOSA(J,1)
          B_BC = DEMCOSA(J,2)
          C_BC = DEMCOSA(J,3)
          ABC_BC = A_BC**2+B_BC**2+C_BC**2

          N_C = CELL_INDEX(J,1)
          Z_N = Z(1) * (HP(N_C) + ELFV(N_C)) + ELFV(N_C)           
          DX_BC = CXY(J,1)-XPDEM(I)
          DY_BC = CXY(J,2)-YPDEM(I)         
          DZ_BC = Z_N - ZPDEM(I)          
                  
          DIS_BC = ABS(A_BC*DX_BC+B_BC*DY_BC+C_BC*DZ_BC)/
     &                              SQRT(ABC_BC) 
          
          IF(DIS_BC .LT. RPDEM(I))THEN
              
            T_BC = (A_BC*DX_BC+B_BC*DY_BC+C_BC*DZ_BC)/ABC_BC
                      
            X2_BC = XPDEM(I)-A_BC*T_BC
            Y2_BC = YPDEM(I)-B_BC*T_BC
            Z2_BC = ZPDEM(I)-C_BC*T_BC
                      
            DO K =  1, CELL_POLYGEN(J)    
c                XN1 = PXY(CELL_SIDE(J,K,3),1)
c                YN1 = PXY(CELL_SIDE(J,K,3),2)
c                XN2 = PXY(CELL_SIDE(J,K,4),1)
c                YN2 = PXY(CELL_SIDE(J,K,4),2)           
      
c                XPGDEM = ( X2_BC - XN1 ) * ( YN2 - YN1 ) - 
c     &                   ( Y2_BC - YN1 ) * ( XN2 - XN1 )

                XN1(K) = PXY(CELL_SIDE(J,K,3),1)
                YN1(K) = PXY(CELL_SIDE(J,K,3),2)
                XN2(K) = PXY(CELL_SIDE(J,K,4),1)
                YN2(K) = PXY(CELL_SIDE(J,K,4),2)           
      
                XPGDEM = ( X2_BC - XN1(K) ) * ( YN2(K) - YN1(K) ) - 
     &                   ( Y2_BC - YN1(K) ) * ( XN2(K) - XN1(K) )
                
                If(XPGDEM .GT. 0.) Then
                    Exit
                Endif
            ENDDO
            
            DO L = 1, CELL_POLYGEN(J)
                N_C = CELL_INDEX(J,L)
                Z_N = Z(1) * (HP(N_C) + ELFV(N_C)) + ELFV(N_C)
                DISTANCE=SQRT( (XPDEM(I)-PXY(N_C,1))**2 + 
     &                         (YPDEM(I)-PXY(N_C,2))**2 +
     &                         (ZPDEM(I)- Z_N)**2 )
              
                IF(DISTANCE .LT. RPDEM(I))THEN            
                    Exit
                Endif             
            ENDDO
            
            If(K .EQ. ( CELL_POLYGEN(J) + 1 ) .OR.
     &         L .NE. ( CELL_POLYGEN(J) + 1 )) Then             

                NP_BCD(I) = NP_BCD(I) + 1
                NCONT_INFO(I,NP_BCD(I)) = 7
                
C                BCDEM_INFO(I,3*NP_BCD(I)-2) = X2_BC+RPDEM(I)/
C     &                                DIS_BC*(X2_BC-XPDEM(I))
C                BCDEM_INFO(I,3*NP_BCD(I)-1) = Y2_BC+RPDEM(I)/
C     &                                DIS_BC*(Y2_BC-YPDEM(I))
C                BCDEM_INFO(I,3*NP_BCD(I)) = Z2_BC+RPDEM(I)/
C     &                                DIS_BC*(Z2_BC-ZPDEM(I))
                Exit
            ENDIF
          ENDIF                  
         ENDDO
        ENDIF

        IF(NP_BCD(I) .GT. 7)THEN
          PRINT*,'ERROR IN BCDEM.F: NP_BCD(I) .GT. 7'
          PRINT*,'I: ',I
          PRINT*,'NP_BCD(I): ',NP_BCD(I)
          PRINT*,'RPDEM(I): ',RPDEM(I)
          PRINT*,'XPDEM(I): ',XPDEM(I)
          PRINT*,'YPDEM(I): ',YPDEM(I)
          PRINT*,'ZPDEM(I): ',ZPDEM(I)
          PRINT*,BCDEM_INFO(I,:)
          PAUSE
        ENDIF    
        
       ENDIF
        
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      !print*,'NCONT_INFO(1,1)=',NCONT_INFO(1,1)
      Return
      END
      
      
      
    !            IF( ( T_INT1 .GT. 0. .AND. T_INT1 .LT. 1. ) .OR.
    !&               ( T_INT2 .GT. 0. .AND. T_INT2 .LT. 1. ) )THEN
    !            IF(1 == 0)THEN      
    !            WRITE(IUPRT,*),
    !            WRITE(IUPRT,*),'*****************************'
    !            WRITE(IUPRT,*),'NNDEM:',NNDEM
    !            WRITE(IUPRT,*),'I,J,K:'
    !            WRITE(IUPRT,*),I,J,CELL_SIDE(J,K,1)
    !            WRITE(IUPRT,*),'P(X,Y,Z):' 
    !            WRITE(IUPRT,*),XPDEM(I),YPDEM(I),ZPDEM(I)
    !            WRITE(IUPRT,*),'RP(I):' 
    !            WRITE(IUPRT,*),RPDEM(I)
    !            WRITE(IUPRT,*),'N1(X,Y,Z):'
    !            WRITE(IUPRT,*),XN1(K),YN1(K),
    !&                         - HP( IEND_EDGE(CELL_SIDE(J,K,1),1) )
    !            XN2(K) = PXY(IEND_EDGE(CELL_SIDE(J,K,1),2),1)
    !            YN2(K) = PXY(IEND_EDGE(CELL_SIDE(J,K,1),2),2)
    !            WRITE(IUPRT,*),'N2(X,Y,Z):'
    !            WRITE(IUPRT,*),XN2(K),YN2(K),
    !&                         - HP( IEND_EDGE(CELL_SIDE(J,K,1),1) )
    !            WRITE(IUPRT,*),'VX,VY,VZ:'
    !            WRITE(IUPRT,*),VX,VY,VZ
    !            WRITE(IUPRT,*),'DISX,DISY,DISZ:'
    !            WRITE(IUPRT,*),DISX,DISY,DISZ  
    !            WRITE(IUPRT,*),'CONST_A,CONST_B,CONST_C:'
    !            WRITE(IUPRT,*),CONST_A,CONST_B,CONST_C
    !            WRITE(IUPRT,*),'DELTA:',DELTA
    !            WRITE(IUPRT,*),'T_INT1,T_INT2:',T_INT1,T_INT2  
    !            WRITE(IUPRT,*),'*****************************'
    !            WRITE(IUPRT,*),
    !            ENDIF  
    !            EXIT
    !               
    !            ENDIF      