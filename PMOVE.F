c-----------------------------------------------------------------------------------c
c                             SUBROUTINE PROGRAM                                    #
C                             DESINGED BY ChenJun                                   #
C                                                 SHANGHAI JIAO TONG UNIVERSITY     #
C                                                 SHANGHAI, CHINA                   #
c-----------------------------------------------------------------------------------#
c 子程序功能：
c 1------颗粒加速度和角加速度计算
c 2------颗粒位移、速度和旋转角速度更新   
c 3------流场时步内所有DEM时步颗粒速度累加求和      
c-----------------------------------------------------------------------------------#
      
      Subroutine PMOVE
      Include './Include/OCERM_INF'
      Double Precision FMDEM, F_GRAVITY, F_BUOYANCY 
      Double Precision A_X, A_Y, A_Z
      Double Precision DOMXDEM, DOMYDEM, DOMZDEM, VRATIO
      DIMENSION A_X(I_PARTICLE), A_Y(I_PARTICLE), A_Z(I_PARTICLE)
             
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,FMDEM,F_GRAVITY,F_BUOYANCY,
!$OMP&                       DOMXDEM,DOMYDEM,DOMZDEM,VRATIO)      
!$OMP DO
      Do I = 1 , I_PARTICLE 

       IF(NSYM_CONT(I) .EQ. 1 .AND. NSYM_MOVE(I) .EQ. 1)THEN          
          
        IF(I_GRA .EQ. 1)THEN  
          F_GRAVITY = - PMDEM(I) * GRAV   !gravity
        ELSE
          F_GRAVITY = 0.0 
        ENDIF
        
C        IF(I_PARTICLE .GT. 1 .AND. I_CV_PORE .EQ. 1 .AND. 
C     &      KPP(I) .NE. 0)THEN
C          CV_PORE = C_VIRTUAL * ( 1.0 + 2.0 * ( 1.0 - 
C     &              PORE(IPP(I),KPP(I)) ) ) / PORE(IPP(I),KPP(I))
C        ELSE
C          CV_PORE = C_VIRTUAL  
C        ENDIF
        
C        IF(KPP(I) .NE. 0 .AND. I_FTP .NE. 0)THEN
C          FMDEM = PMDEM(I) + CV_PORE(I) * RMEAN(IPP(I),KPP(I)) * 
C     &            V_PARTICLE(I)   !total mass = particle mass + fictitious mass
CC          F_BUOYANCY = RMEAN(IPP(I),KPP(I)) * V_PARTICLE(I) * GRAV   !buoyancy
C        ELSE
C          FMDEM = PMDEM(I)
C        ENDIF        
        
        IF(I_FTP .NE. 0)THEN   !the force that fluid acting on particle is neglect            
         IF(MOD_FA .EQ. 'IMPLICIT')THEN            
          FMDEM = PMDEM(I) + CV_PORE(I) * RMEAN(IPP(I),KPP(I)) * 
     &            V_PARTICLE(I)   !DEM隐式处理附加质量力项,total mass = particle mass + fictitious mass
         ELSE
          FMDEM = PMDEM(I)
         ENDIF
        ELSE          
          FDXDEM(I) = 0.0
          FDYDEM(I) = 0.0
          FDZDEM(I) = 0.0
          
          IF(MOD_FA .EQ. 'IMPLICIT')THEN   !DEM隐式处理附加质量力项   
           FAFXDEM(I) = 0.0   
           FAFYDEM(I) = 0.0 
           FAFZDEM(I) = 0.0  
          ELSE   !DEM显式处理附加质量力项 
           FAXDEM(I) = 0.0   
           FAYDEM(I) = 0.0 
           FAZDEM(I) = 0.0                       
          ENDIF           
          
          F_PG_X(I) = 0.0
          F_PG_Y(I) = 0.0
          F_PG_Z(I) = 0.0
          FMDEM = PMDEM(I)          
        ENDIF                  
        
        IF(I_TRAN .EQ. 1) THEN
         
         !velocity accelaration   
C         A_X(I) = (FNXDEM(I) + FTXDEM(I) + FDXDEM(I) + F_PG_X(I)) / 
C     &             FMDEM   
C         A_Y(I) = (FNYDEM(I) + FTYDEM(I) + FDYDEM(I) + F_PG_Y(I)) / 
C     &             FMDEM
C         A_Z(I) = (FNZDEM(I) + FTZDEM(I) + FDZDEM(I) + F_PG_Z(I) + 
C     &             F_GRAVITY) / FMDEM       
C!    &             F_GRAVITY + F_BUOYANCY) / FMDEM 
         
         IF(MOD_FA .EQ. 'IMPLICIT')THEN   !DEM隐式处理附加质量力项   
          A_X(I) = (FNXDEM(I) + FTXDEM(I) + FDXDEM(I) + F_PG_X(I) +  
     &              FAFXDEM(I) ) / FMDEM   
          A_Y(I) = (FNYDEM(I) + FTYDEM(I) + FDYDEM(I) + F_PG_Y(I) +  
     &              FAFYDEM(I) ) / FMDEM
          A_Z(I) = (FNZDEM(I) + FTZDEM(I) + FDZDEM(I) + F_PG_Z(I) + 
     &              FAFZDEM(I) + F_GRAVITY) / FMDEM       
!    &             F_GRAVITY + F_BUOYANCY) / FMDEM 
         ELSE   !DEM显式处理附加质量力项         
          A_X(I) = (FNXDEM(I) + FTXDEM(I) + FDXDEM(I) + F_PG_X(I) +  
     &              FAXDEM(I) ) / FMDEM   
          A_Y(I) = (FNYDEM(I) + FTYDEM(I) + FDYDEM(I) + F_PG_Y(I) +  
     &              FAYDEM(I) ) / FMDEM
          A_Z(I) = (FNZDEM(I) + FTZDEM(I) + FDZDEM(I) + F_PG_Z(I) + 
     &              FAZDEM(I) + F_GRAVITY) / FMDEM   
!    &             F_GRAVITY + F_BUOYANCY) / FMDEM                    
         ENDIF              
            
         XPDEM(I) = XPDEM(I) + UPDEM(I) * TDEM + 0.5 * A_X(I) * TDEM**2   !new position
         YPDEM(I) = YPDEM(I) + VPDEM(I) * TDEM + 0.5 * A_Y(I) * TDEM**2
         ZPDEM(I) = ZPDEM(I) + WPDEM(I) * TDEM + 0.5 * A_Z(I) * TDEM**2
        
         UPDEM(I) = UPDEM(I) + A_X(I) * TDEM   !new velocity
         VPDEM(I) = VPDEM(I) + A_Y(I) * TDEM
         WPDEM(I) = WPDEM(I) + A_Z(I) * TDEM          
         
        ENDIF

        IF(I_FCOUPLING .EQ. 1)THEN
          UPSUMDEM(I) = UPSUMDEM(I) + UPDEM(I)
          VPSUMDEM(I) = VPSUMDEM(I) + VPDEM(I)
          WPSUMDEM(I) = WPSUMDEM(I) + WPDEM(I)  
        ENDIF 
         
        IF(I_ROT .EQ. 1) THEN
            
 	   DOMXDEM = TPXDEM(I) / PIPDEM(I) * TDEM   !angular accelaration
         DOMYDEM = TPYDEM(I) / PIPDEM(I) * TDEM
         DOMZDEM = TPZDEM(I) / PIPDEM(I) * TDEM
         
         OMXDEM(I) = OMXDEM(I) + DOMXDEM   !new angular velocity
         OMYDEM(I) = OMYDEM(I) + DOMYDEM
         OMZDEM(I) = OMZDEM(I) + DOMZDEM   
         
        ENDIF
        
       ENDIF
       
      Enddo
!$OMP END DO
!$OMP END PARALLEL     
C      write(*,'(A16,F8.4)'),'MAX_VELOCITY = ',DMAX1( DMAX1( 
C     &      maxval(ABS(UPDEM)),maxval(ABS(VPDEM)) ),maxval(ABS(WPDEM)) )
      IF(1 .EQ. 0)THEN       

         UMAX_TEMP = 0.  
         VMAX_TEMP = 0. 
         WMAX_TEMP = 0.           
          
       Do I = 1 , I_PARTICLE 
        IF(NSYM_CONT(I) .EQ. 1 .AND. NSYM_MOVE(I) .EQ. 1)THEN 
         UMAX_TEMP = DMAX1(UMAX_TEMP,ABS(UPDEM(I)))   
         VMAX_TEMP = DMAX1(VMAX_TEMP,ABS(VPDEM(I))) 
         WMAX_TEMP = DMAX1(WMAX_TEMP,ABS(WPDEM(I))) 
        ENDIF
       ENDDO
       
       write(*,'(A16,F8.4)'),'U_MAX = ', UMAX_TEMP     
       write(*,'(A16,F8.4)'),'V_MAX = ', VMAX_TEMP  
       write(*,'(A16,F8.4)'),'W_MAX = ', WMAX_TEMP  
       
      ENDIF      
      
      Return
      End
      
c         VRATIO = (1.6 * RPDEM(I) / TDEM) / DMAX1( 1.0E-20, 
c     &            SQRT(UPDEM(I)**2 + VPDEM(I)**2 + WPDEM(I)**2) )
c         IF(VRATIO .LE. 1.0)THEN
c          PRINT*,'warning1 in PMOVE.F: the particle velocity is 
c     &            too high'

c          WRITE(IUPRT,*)'warning1 in PMOVE.F: the particle velocity is 
c     &            too high'
c          WRITE(IUPRT,*)'I=',I,'X=',XPDEM(I),'Y=',YPDEM(I),'Z=',ZPDEM(I)
c          WRITE(IUPRT,*)'I=',I,'U=',UPDEM(I),'V=',VPDEM(I),'W=',WPDEM(I)
          
c          UPDEM(I) = VRATIO * UPDEM(I)
c          VPDEM(I) = VRATIO * VPDEM(I)
c          WPDEM(I) = VRATIO * WPDEM(I)
c         ENDIF      
