      Subroutine PARTICLE_INFO    
        Include './Include/OCERM_INF'
        
C--- 初始化        
        DO I = 1, I_PARTICLE
            
          N_P2(I) = 0
          
          N_SYNFA(I) = 0
          UFPDEM_PRE(I) = 0.
          VFPDEM_PRE(I) = 0.
          WFPDEM_PRE(I) = 0.     
          
          UFDEM_PRE(I) = 0.
          VFDEM_PRE(I) = 0.
          WFDEM_PRE(I) = 0.           
          
          DO K = 1, I_CPN
            N_PP2(I,K) = 0
            DELTAT_PP1X(I,K) = 0.0
            DELTAT_PP1Y(I,K) = 0.0
            DELTAT_PP1Z(I,K) = 0.0
          ENDDO
            
          UPSUMDEM(I) = 0.0
          VPSUMDEM(I) = 0.0
          WPSUMDEM(I) = 0.0 
          
          DO J = 1, 7
            DEL_DEM(I,J) = 0.
          ENDDO
          
          DEL_HIS(I) = 0.
          
        ENDDO        
        
        Open(1001,file='./Include/PARTICLE_SOFT.DAT')
	  Read(1001,*)
	  
          	  
        Do I = 1 , NPARTICLE
		 Read(1001,'(15E28.20)') XPDEM(I),YPDEM(I),ZPDEM(I),UPDEM(I),VPDEM(I),
     &                 WPDEM(I),OMXDEM(I),OMYDEM(I),OMZDEM(I),RPMEAN(I),
     &			     PMDEM(I),RPDEM(I),PNIUPDEM(I),EPDEM(I),GPDEM(I)
        
        Enddo
        
        IF(NPARTICLE1 .NE. 0)THEN
            
          Open(1101,file='./Include/PARTICLE_SOFT1.DAT')
	    Read(1101,*)            
          Do I = NPARTICLE + 1 , I_PARTICLE
		 Read(1101,'(15E28.20)') XPDEM(I),YPDEM(I),ZPDEM(I),UPDEM(I),VPDEM(I),
     &                 WPDEM(I),OMXDEM(I),OMYDEM(I),OMZDEM(I),RPMEAN(I),
     &			     PMDEM(I),RPDEM(I),PNIUPDEM(I),EPDEM(I),GPDEM(I)
          Enddo 
          Close(1101)
          
        ENDIF        

        !是否随机分配颗粒（是：IP_RANDOM = 1；否：其他）
        
        If(DEMRESTART .NE. 'COLD START') IP_RANDOM = 0
        
        IF(IP_RANDOM .EQ. 1)THEN
        
        DO K = 1, 100   !随机分配颗粒
        
         DO I = 1, I_PARTICLE
      
          call RANDOM_NUMBER(x)
          J=floor(x*FLOAT(I_PARTICLE))+1
          
          
          TEMP = XPDEM(I)
          XPDEM(I) = XPDEM(J)
          XPDEM(J) = TEMP
          
          TEMP = YPDEM(I)
          YPDEM(I) = YPDEM(J)
          YPDEM(J) = TEMP
          
          TEMP = ZPDEM(I)
          ZPDEM(I) = ZPDEM(J)
          ZPDEM(J) = TEMP
          
          TEMP = UPDEM(I)
          UPDEM(I) = UPDEM(J)
          UPDEM(J) = TEMP
          
          TEMP = VPDEM(I)
          VPDEM(I) = VPDEM(J)
          VPDEM(J) = TEMP
          
          TEMP = WPDEM(I)
          WPDEM(I) = WPDEM(J)
          WPDEM(J) = TEMP
          
          TEMP = OMXDEM(I)
          OMXDEM(I) = OMXDEM(J)
          OMXDEM(J) = TEMP
          
          TEMP = OMYDEM(I)
          OMYDEM(I) = OMYDEM(J)
          OMYDEM(J) = TEMP
          
          TEMP = OMZDEM(I)
          OMZDEM(I) = OMZDEM(J)
          OMZDEM(J) = TEMP
          
          TEMP = RPMEAN(I)
          RPMEAN(I) = RPMEAN(J)
          RPMEAN(J) = TEMP
          
          TEMP = PMDEM(I)
          PMDEM(I) = PMDEM(J)
          PMDEM(J) = TEMP
          
          TEMP = RPDEM(I)
          RPDEM(I) = RPDEM(J)
          RPDEM(J) = TEMP
          
          TEMP = PNIUPDEM(I)
          PNIUPDEM(I) = PNIUPDEM(J)
          PNIUPDEM(J) = TEMP
          
          TEMP = EPDEM(I)
          EPDEM(I) = EPDEM(J)
          EPDEM(J) = TEMP
          
          TEMP = GPDEM(I)
          GPDEM(I) = GPDEM(J)
          GPDEM(J) = TEMP
          
         ENDDO  
        ENDDO
        
        ENDIF
        
        !J=0
        !DO I = 1, 10222
        !  IF(RPDEM(I) .GT. 0.0006)THEN
        !    J = J + 1
        !    PRINT*,J
        !  ENDIF
        !ENDDO      

        !PRINT*,'INPUT RHOP_CONST'
        !READ(*,*)RHOP_CONST

        DO I = 1, I_PARTICLE
        
           V_PARTICLE(I) = ( 4. / 3. ) *  ACOS(-1.0) * RPDEM(I)**3
           
           IF(RHO_PARTICLE .LT. 0.000001)THEN
             PMDEM(I) = RPMEAN(I) * V_PARTICLE(I)  
           ELSE
             PMDEM(I) = RHO_PARTICLE * V_PARTICLE(I)
             !print*, YPDEM(I),zPDEM(I)
             !PMDEM(I) = RPMEAN(I) * V_PARTICLE(I)     
           ENDIF
           
           PIPDEM(I) = 2. / 5. * PMDEM(I) * RPDEM(I)**2
           
           NSYM_CONT(I) = 1
           DPDEM(I) = 2.0 * RPDEM(I)
           UPDEM(I)=0.
           VPDEM(I)=0.
           WPDEM(I)=0.
           OMXDEM(I)=0.
           OMYDEM(I)=0.
           OMZDEM(I)=0.
           
           NSYM_MOVE(I) = 1
c           IF(XPDEM(I) .LT. 0.06 .OR. XPDEM(I) .GT. 0.09 .OR. 
c     &        YPDEM(I) .GT. 0.03 .OR. ZPDEM(I) .LT. -0.1025)THEN
              
c              NSYM_MOVE(I) = 0
              
c           ENDIF 
 
        ENDDO
        
        DPMAX = 2 * MAXVAL(RPDEM)
        
        Close(1001)
        
        IF(IP_RANDOM .EQ. 1)THEN
          Open(1102,file='./Include/PARTICLE_SOFT_RANDOM.DAT')
	    WRITE(1102,*)            
          Do I = 1 , I_PARTICLE
		 WRITE(1102,'(15E28.20)')XPDEM(I),YPDEM(I),ZPDEM(I),UPDEM(I),VPDEM(I),
     &                 WPDEM(I),OMXDEM(I),OMYDEM(I),OMZDEM(I),RPMEAN(I),
     &			     PMDEM(I),RPDEM(I),PNIUPDEM(I),EPDEM(I),GPDEM(I)
          Enddo 
          Close(1102)  
        ENDIF
        
        N_TOTAL = 0
        
        DO I = 1, I_PARTICLE
C         IF(RPDEM(I) .GT. 0.0015 .OR. RPDEM(I) .LT. 0.00075)THEN 
             
          DO J = 1, N_PSLEV
            
            I_PSLEV(I) = J  
              
            IF(2.0 * RPDEM(I) .LE. DELBM(J)*1.000001)THEN
              
              I_PSLEV(I) = J
              
              EXIT
              
            ENDIF
                        
          ENDDO
 
          IF(J .EQ. N_PSLEV+1)THEN
            PRINT*,'ERROR OF DELBM IN DEMSETTING.F!'  
            PRINT*,2.0 * RPDEM(I),DELBM(J-1)
            PAUSE
          ENDIF          
          
C         ELSE
             
C          I_PSLEV(I) = 0
          
C         ENDIF
          
          IF(I_PSLEV(I) .GT. 1)THEN
              
              IP_LL_INDEX(I) = N_TOTAL
              
              NPT_LL(I) = 5 * 
     &          (INT((RPDEM(I)+MINVAL(RPDEM))/MINVAL(RPDEM)))**2
              
              N_TOTAL = N_TOTAL + NPT_LL(I)
  
              
              IF(N_TOTAL .GT. NP_CONTACT)THEN
                  
                PRINT*,'ERROR OF NP_CONTACT IN PARTICLE_INFO!'  
                PAUSE
                
              ENDIF
              
          ENDIF
          
        ENDDO              
      
      End
      
 
