c-----------------------------------------------------------------------------------c
c                             SUBROUTINE PROGRAM                                    #
C                             DESINGED BY ChenJun                                   #
C                                                 SHANGHAI JIAO TONG UNIVERSITY     #
C                                                 SHANGHAI, CHINA                   #
c-----------------------------------------------------------------------------------#
c 子程序功能：
c 1------颗粒背景网格寻址及颗粒存储       
c-----------------------------------------------------------------------------------#

      Subroutine IPJPKPDEM2
      Include './Include/OCERM_INF'
c      Double Precision 
      
!      IF(PSMOD .EQ. 'MOD2')THEN   
!       DO L = 1, N_PSLEV    
!        Do K = 1, NBZ(L)
!          Do J = 1, NBY(L)
!c!$OMP DO	              
!            Do I = 1, NBX(L)  
!              IF(L .EQ. 1)THEN
!                NBCP1(I,J,K) = 0
!              ELSEIF(L .EQ. 2)THEN
!                NBCP2(I,J,K) = 0
!              ELSEIF(L .EQ. 3)THEN
!                NBCP3(I,J,K) = 0
!              ELSEIF(L .EQ. 4)THEN
!                NBCP4(I,J,K) = 0
!              ELSEIF(L .EQ. 5)THEN
!                NBCP5(I,J,K) = 0
!              ENDIF
!            ENDDO
!c!$OMP END DO
!          ENDDO
!        ENDDO
!       ENDDO 
!      ENDIF       
      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,L,N)  
C!$OMP&         REDUCTION(+:NBCP1,NBCP2,NBCP3,NBCP4,NBCP5)
C!$OMP&         REDUCTION(.EQV.:IBCP1,IBCP2,IBCP3,IBCP4,IBCP5)

C      IF(PSMOD .EQ. 'MOD2')THEN  
!$OMP DO          
      DO I = 1, NBL    
        NBCP( I ) = 0
      ENDDO 
!$OMP END DO       
C      ENDIF   

!$OMP DO     
      Do I = 1 , I_PARTICLE 
       IF(NSYM_CONT(I) .EQ. 1 .AND. I_PSLEV(I) .NE. 0)THEN   
        N = I_PSLEV(I)   !颗粒所在层级  
          
        DO J = N, N_PSLEV   !识别颗粒在所在层级和更高层级背景网格上所属网格单元编号
            
          IF(XPDEM(I) .GE. XMV(IP(I,J),J) .AND. 
     &       XPDEM(I) .LT. XMV(IP(I,J)+1,J))THEN
            IP(I,J) = IP(I,J)
          ELSEIF( XPDEM(I) .GE. XMV( MAX(IP(I,J)-1,1),J ) .AND. 
     &       XPDEM(I) .LT. XMV(IP(I,J),J) )THEN
            IP(I,J) = IP(I,J) - 1 
          ELSEIF(XPDEM(I) .GE. XMV(IP(I,J)+1,J) .AND. 
     &       XPDEM(I) .LT. XMV( MIN(IP(I,J)+2,NBX(J)+1),J ) )THEN
            IP(I,J) = IP(I,J) + 1     
          ELSE
            Write (IUPRT,*) 'Particle I= ',I,' is vanished in 
     &                       IPJPKPDEM.F'    
            Write (IUPRT,*) 'XPDEM(I)=',XPDEM(I)
            Write (IUPRT,*) 'YPDEM(I)=',YPDEM(I)
            Write (IUPRT,*) 'ZPDEM(I)=',ZPDEM(I)
            NSYM_CONT(I) = 0
          ENDIF

          IF(YPDEM(I) .GE. YMV(JP(I,J),J) .AND. 
     &       YPDEM(I) .LT. YMV(JP(I,J)+1,J))THEN
            JP(I,J) = JP(I,J)
          ELSEIF(YPDEM(I) .GE. YMV( MAX(JP(I,J)-1,1),J ) .AND. 
     &       YPDEM(I) .LT. YMV(JP(I,J),J))THEN
            JP(I,J) = JP(I,J) - 1 
          ELSEIF(YPDEM(I) .GE. YMV(JP(I,J)+1,J) .AND. 
     &       YPDEM(I) .LT. YMV( MIN(JP(I,J)+2,NBY(J)+1),J ) )THEN
            JP(I,J) = JP(I,J) + 1     
          ELSE
            Write (IUPRT,*) 'Particle I= ',I,' is vanished in 
     &                       IPJPKPDEM.F'    
            Write (IUPRT,*) 'XPDEM(I)=',XPDEM(I)
            Write (IUPRT,*) 'YPDEM(I)=',YPDEM(I)
            Write (IUPRT,*) 'ZPDEM(I)=',ZPDEM(I)
            NSYM_CONT(I) = 0
          ENDIF          

          IF(ZPDEM(I) .GE. ZMV(KP(I,J),J) .AND. 
     &       ZPDEM(I) .LT. ZMV(KP(I,J)+1,J))THEN
            KP(I,J) = KP(I,J)
          ELSEIF(ZPDEM(I) .GE. ZMV( MAX(KP(I,J)-1,1),J ) .AND. 
     &       ZPDEM(I) .LT. ZMV(KP(I,J),J))THEN
            KP(I,J) = KP(I,J) - 1 
          ELSEIF(ZPDEM(I) .GE. ZMV(KP(I,J)+1,J) .AND. 
     &       ZPDEM(I) .LT. ZMV( MIN(KP(I,J)+2,NBZ(J)+1),J ) )THEN
            KP(I,J) = KP(I,J) + 1     
          ELSE
            Write (IUPRT,*) 'Particle I= ',I,' is vanished in 
     &                       IPJPKPDEM.F'     
            Write (IUPRT,*) 'XPDEM(I)=',XPDEM(I)
            Write (IUPRT,*) 'YPDEM(I)=',YPDEM(I)
            Write (IUPRT,*) 'ZPDEM(I)=',ZPDEM(I)
            NSYM_CONT(I) = 0
          ENDIF       

        ENDDO                                          
       ENDIF
      Enddo
      
!$OMP END DO
!$OMP END PARALLEL           

      Do I = 1 , I_PARTICLE 
          
       IF(NSYM_CONT(I) .EQ. 1 .AND. I_PSLEV(I) .NE. 0)THEN
           
         N = I_PSLEV(I)  
         
         N_TEMP = NBCP_INDEX(N) + (KP(I,N)-1)*NBX(N)*NBY(N) + 
     &            NBX(N)*(JP(I,N)-1) + IP(I,N)   !颗粒在N层级所属网格单元在NBCP数组中的存储位置
          
         NBCP( N_TEMP ) = NBCP( N_TEMP ) + 1   !颗粒在N层级所属网格单元颗粒数+1
         
         IF(NBCP( N_TEMP ) .GT. IBCPMAX(N))THEN
            PRINT*,'Error4 of IBCPMAX in IPJPKPDEM.F!'
            WRITE(IUPRT,*)
            WRITE(IUPRT,'(3A14)') 'XPDEM', 'YPDEM', 'ZPDEM' 
            WRITE(IUPRT,'(3F14.7)') XPDEM(I),YPDEM(I),ZPDEM(I)
            PAUSE
         ENDIF          
         
         IBCP( (N_TEMP - 1) * IBCPMAX + NBCP( N_TEMP ) ) = I   !将颗粒I单元编号保存到颗粒在所在层级所属网格单元
         
       ENDIF
       
      ENDDO
      
      Return 

5200  Format (80A1)        
      
      End        
C-------------------------------------------------------------------  
         ! Do I = 1 , I_PARTICLE 
         ! 
         !  DO L4 = I_PSLEV(I), N_PSLEV
         !
         !   DO L3 = MAX(KP(I,L4)-1,1), MIN(KP(I,L4)+1,NBZ(L4))
         !   
         !    DO L2 = MAX(JP(I,L4)-1,1), MIN(JP(I,L4)+1,NBY(L4))
         !    
         !     DO L1 = MAX(IP(I,L4)-1,1), MIN(IP(I,L4)+1,NBX(L4))    
         !     
         !      N_TEMP = NBCP_INDEX(L4) + (L3-1)*NBX*NBY + 
         !&              NBX*(L2-1) + L1
         ! 
         !      NP_CAND(I,2) = NP_CAND_TOTAL
         !  
         !      DO J = 1, NBCP( N_TEMP )   
         !   
         !       NP_CAND1(I) = NP_CAND1(I) + 1
         !   
         !       IP_CAND(NP_CAND1(I)) = J
         !   
         !       NP_CAND1(J) = NP_CAND1(I) + 1
         !   
         !       IP_CAND(NP_CAND1(J)) = J
         !   
         !       IF(L4 .GT. I_PSLEV(I))THEN
         !       
         !       
         !      
         !       NP_CAND(I,1) = NP_CAND(I,1) + 1
         !   
         !       NP_CAND_TOTAL = NP_CAND_TOTAL + 1
         !   
         !       IP_CAND(NP_CAND_TOTAL) = J
         !   
         !      ENDDO
         !  
         !     ENDDO
         ! 
         !    ENDDO
         !
         !   ENDDO
         !
         !  ENDDO
         !
         ! ENDDO
         !
         !
         !    
         !       NP_CAND (I_PARTICLE,2)：存储（1：候选颗粒个数；2： IP_CAND 中颗粒存储起始位置）
         !       IP_CAND (I_PARTICLE*20)：存储所有候选颗粒编号
  
      
C-------------------------------------------------------------------           
         !
         ! Do I = 1 , I_PARTICLE 
         ! 
         !  IF(NSYM_CONT(I) .EQ. 1 .AND. I_PSLEV(I) .NE. 0)THEN   
         !  
         !   N = I_PSLEV(I)   !颗粒所在层级
         !
         !   IF(N .EQ. 1)THEN
         !   
         !     NBCP1(IP(I,N),JP(I,N),KP(I,N)) = 
         !&      NBCP1(IP(I,N),JP(I,N),KP(I,N)) + 1
         !
         !     IBCP1(IP(I,N),JP(I,N),KP(I,N),
         !&      NBCP1(IP(I,N),JP(I,N),KP(I,N))) = I
         ! 
         !     IF(NBCP1(IP(I,N),JP(I,N),KP(I,N)) .GT. IBCPMAX)THEN
         !       PRINT*,'Error4 of IBCPMAX in IPJPKPDEM.F!'
         !       WRITE(IUPRT,*)
         !       WRITE(IUPRT,'(3A14)') 'XPDEM', 'YPDEM', 'ZPDEM' 
         !       WRITE(IUPRT,'(3F14.7)') XPDEM(I),YPDEM(I),ZPDEM(I)
         !       PAUSE
         !     ENDIF            
         ! 
         !   ELSEIF(N .EQ. 2)THEN
         !   
         !     NBCP2(IP(I,N),JP(I,N),KP(I,N)) = 
         !&      NBCP2(IP(I,N),JP(I,N),KP(I,N)) + 1
         !
         !     IBCP2(IP(I,N),JP(I,N),KP(I,N),
         !&      NBCP2(IP(I,N),JP(I,N),KP(I,N))) = I   
         ! 
         !     IF(NBCP2(IP(I,N),JP(I,N),KP(I,N)) .GT. IBCPMAX)THEN
         !       PRINT*,'Error4 of IBCPMAX in IPJPKPDEM.F!'
         !       WRITE(IUPRT,*)
         !       WRITE(IUPRT,'(3A14)') 'XPDEM', 'YPDEM', 'ZPDEM' 
         !       WRITE(IUPRT,'(3F14.7)') XPDEM(I),YPDEM(I),ZPDEM(I)        
         !       PAUSE
         !     ENDIF            
         ! 
         !   ELSEIF(N .EQ. 3)THEN
         !   
         !     NBCP3(IP(I,N),JP(I,N),KP(I,N)) = 
         !&      NBCP3(IP(I,N),JP(I,N),KP(I,N)) + 1
         !
         !     IBCP3(IP(I,N),JP(I,N),KP(I,N),
         !&      NBCP3(IP(I,N),JP(I,N),KP(I,N))) = I  
         ! 
         !     IF(NBCP3(IP(I,N),JP(I,N),KP(I,N)) .GT. IBCPMAX)THEN
         !       PRINT*,'Error4 of IBCPMAX in IPJPKPDEM.F!'
         !       WRITE(IUPRT,*)
         !       WRITE(IUPRT,'(3A14)') 'XPDEM', 'YPDEM', 'ZPDEM' 
         !       WRITE(IUPRT,'(3F14.7)') XPDEM(I),YPDEM(I),ZPDEM(I)
         !       PAUSE
         !     ENDIF            
         ! 
         !   ELSEIF(N .EQ. 4)THEN
         ! 
         !     NBCP4(IP(I,N),JP(I,N),KP(I,N)) = 
         !&      NBCP4(IP(I,N),JP(I,N),KP(I,N)) + 1
         !
         !     IBCP4(IP(I,N),JP(I,N),KP(I,N),
         !&      NBCP4(IP(I,N),JP(I,N),KP(I,N))) = I  
         ! 
         !     IF(NBCP4(IP(I,N),JP(I,N),KP(I,N)) .GT. IBCPMAX)THEN
         !       PRINT*,'Error4 of IBCPMAX in IPJPKPDEM.F!'
         !       WRITE(IUPRT,*)
         !       WRITE(IUPRT,'(3A14)') 'XPDEM', 'YPDEM', 'ZPDEM' 
         !       WRITE(IUPRT,'(3F14.7)') XPDEM(I),YPDEM(I),ZPDEM(I)
         !       PAUSE
         !     ENDIF            
         ! 
         !   ELSEIF(N .EQ. 5)THEN
         !
         !     NBCP5(IP(I,N),JP(I,N),KP(I,N)) = 
         !&      NBCP5(IP(I,N),JP(I,N),KP(I,N)) + 1
         !
         !     IBCP5(IP(I,N),JP(I,N),KP(I,N),
         !&      NBCP5(IP(I,N),JP(I,N),KP(I,N))) = I  
         ! 
         !     IF(NBCP5(IP(I,N),JP(I,N),KP(I,N)) .GT. IBCPMAX)THEN
         !       PRINT*,'Error4 of IBCPMAX in IPJPKPDEM.F!'
         !       WRITE(IUPRT,*)
         !       WRITE(IUPRT,'(3A14)') 'XPDEM', 'YPDEM', 'ZPDEM' 
         !       WRITE(IUPRT,'(3F14.7)') XPDEM(I),YPDEM(I),ZPDEM(I)
         !       PAUSE
         !     ENDIF            
         ! 
         !   ENDIF        
         !
         !  ENDIF
         !
         ! Enddo        

      
         !   IF(J .EQ. 1)THEN
         !  
         !     NBCP1(IP(I,J),JP(I,J),KP(I,J)) = 
         !&      NBCP1(IP(I,J),JP(I,J),KP(I,J)) + 1
         !
         !     IBCP1(IP(I,J),JP(I,J),KP(I,J),
         !&      NBCP1(IP(I,J),JP(I,J),KP(I,J))) = I
         ! 
         !     IF(NBCP1(IP(I,J),JP(I,J),KP(I,J)) .GT. IBCPMAX1)THEN
         !       PRINT*,'Error4 of IBCPMAX1 in IPJPKPDEM.F!'
         !       WRITE(IUPRT,*)
         !       WRITE(IUPRT,'(3A14)') 'XPDEM', 'YPDEM', 'ZPDEM' 
         !       WRITE(IUPRT,'(3F14.7)') XPDEM(I),YPDEM(I),ZPDEM(I)
         !       PAUSE
         !     ENDIF            
         ! 
         !   ELSEIF(J .EQ. 2)THEN
         !   
         !     NBCP2(IP(I,J),JP(I,J),KP(I,J)) = 
         !&      NBCP2(IP(I,J),JP(I,J),KP(I,J)) + 1
         !
         !     IBCP2(IP(I,J),JP(I,J),KP(I,J),
         !&      NBCP2(IP(I,J),JP(I,J),KP(I,J))) = I   
         ! 
         !     IF(NBCP2(IP(I,J),JP(I,J),KP(I,J)) .GT. IBCPMAX2)THEN
         !       PRINT*,'Error4 of IBCPMAX2 in IPJPKPDEM.F!'
         !       WRITE(IUPRT,*)
         !       WRITE(IUPRT,'(3A14)') 'XPDEM', 'YPDEM', 'ZPDEM' 
         !       WRITE(IUPRT,'(3F14.7)') XPDEM(I),YPDEM(I),ZPDEM(I)        
         !       PAUSE
         !     ENDIF            
         ! 
         !   ELSEIF(J .EQ. 3)THEN
         !   
         !     NBCP3(IP(I,J),JP(I,J),KP(I,J)) = 
         !&      NBCP3(IP(I,J),JP(I,J),KP(I,J)) + 1
         !
         !     IBCP3(IP(I,J),JP(I,J),KP(I,J),
         !&      NBCP3(IP(I,J),JP(I,J),KP(I,J))) = I  
         ! 
         !     IF(NBCP3(IP(I,J),JP(I,J),KP(I,J)) .GT. IBCPMAX3)THEN
         !       PRINT*,'Error4 of IBCPMAX3 in IPJPKPDEM.F!'
         !       WRITE(IUPRT,*)
         !       WRITE(IUPRT,'(3A14)') 'XPDEM', 'YPDEM', 'ZPDEM' 
         !       WRITE(IUPRT,'(3F14.7)') XPDEM(I),YPDEM(I),ZPDEM(I)
         !       PAUSE
         !     ENDIF            
         ! 
         !   ELSEIF(J .EQ. 4)THEN
         ! 
         !     NBCP4(IP(I,J),JP(I,J),KP(I,J)) = 
         !&      NBCP4(IP(I,J),JP(I,J),KP(I,J)) + 1
         !
         !     IBCP4(IP(I,J),JP(I,J),KP(I,J),
         !&      NBCP4(IP(I,J),JP(I,J),KP(I,J))) = I  
         ! 
         !     IF(NBCP4(IP(I,J),JP(I,J),KP(I,J)) .GT. IBCPMAX4)THEN
         !       PRINT*,'Error4 of IBCPMAX4 in IPJPKPDEM.F!'
         !       WRITE(IUPRT,*)
         !       WRITE(IUPRT,'(3A14)') 'XPDEM', 'YPDEM', 'ZPDEM' 
         !       WRITE(IUPRT,'(3F14.7)') XPDEM(I),YPDEM(I),ZPDEM(I)
         !       PAUSE
         !     ENDIF            
         ! 
         !   ELSEIF(J .EQ. 5)THEN
         !
         !     NBCP5(IP(I,J),JP(I,J),KP(I,J)) = 
         !&      NBCP5(IP(I,J),JP(I,J),KP(I,J)) + 1
         !
         !     IBCP5(IP(I,J),JP(I,J),KP(I,J),
         !&      NBCP5(IP(I,J),JP(I,J),KP(I,J))) = I  
         ! 
         !     IF(NBCP5(IP(I,J),JP(I,J),KP(I,J)) .GT. IBCPMAX5)THEN
         !       PRINT*,'Error4 of IBCPMAX5 in IPJPKPDEM.F!'
         !       WRITE(IUPRT,*)
         !       WRITE(IUPRT,'(3A14)') 'XPDEM', 'YPDEM', 'ZPDEM' 
         !       WRITE(IUPRT,'(3F14.7)') XPDEM(I),YPDEM(I),ZPDEM(I)
         !       PAUSE
         !     ENDIF            
         ! 
         !   ENDIF                    

      !IF(N_PSLEV .EQ. 1 .AND. SUM(NBCP1) .NE. I_PARTICLE)THEN
      !  PRINT*,'ERROR in IPJPKPDEM.F: SUM(NBCP1) .NE. I_PARTICLE!'
      !  PRINT*,'SUM(NBCP1)=',SUM(NBCP1)
      !  PRINT*,'I_PARTICLE=',I_PARTICLE
      !ELSEIF(N_PSLEV .EQ. 2 .AND. SUM(NBCP2) .NE. I_PARTICLE)THEN
      !  PRINT*,'ERROR in IPJPKPDEM.F: SUM(NBCP2) .NE. I_PARTICLE!'
      !  PRINT*,'SUM(NBCP2)=',SUM(NBCP2)
      !  PRINT*,'I_PARTICLE=',I_PARTICLE   
      !ELSEIF(N_PSLEV .EQ. 3 .AND. SUM(NBCP3) .NE. I_PARTICLE)THEN
      !  PRINT*,'ERROR in IPJPKPDEM.F: SUM(NBCP3) .NE. I_PARTICLE!'
      !  PRINT*,'SUM(NBCP3)=',SUM(NBCP3)
      !  PRINT*,'I_PARTICLE=',I_PARTICLE
      !ELSEIF(N_PSLEV .EQ. 4 .AND. SUM(NBCP4) .NE. I_PARTICLE)THEN
      !  PRINT*,'ERROR in IPJPKPDEM.F: SUM(NBCP4) .NE. I_PARTICLE!'
      !  PRINT*,'SUM(NBCP4)=',SUM(NBCP4)
      !  PRINT*,'I_PARTICLE=',I_PARTICLE        
      !ELSEIF(N_PSLEV .EQ. 5 .AND. SUM(NBCP5) .NE. I_PARTICLE)THEN
      !  PRINT*,'ERROR in IPJPKPDEM.F: SUM(NBCP5) .NE. I_PARTICLE!'
      !  PRINT*,'SUM(NBCP5)=',SUM(NBCP5)
      !  PRINT*,'I_PARTICLE=',I_PARTICLE
      !ENDIF
     
        