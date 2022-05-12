c-----------------------------------------------------------------------------------c
c                             SUBROUTINE PROGRAM                                    #
C                             DESINGED BY ChenJun                                   #
C                                                 SHANGHAI JIAO TONG UNIVERSITY     #
C                                                 SHANGHAI, CHINA                   #
c-----------------------------------------------------------------------------------#
c 子程序功能：
c 1------CFD-DEM耦合计算      
c-----------------------------------------------------------------------------------#
      Subroutine DEMM
      Include './Include/OCERM_INF'
      Integer IIDEM,I_PRINT,IJK(10)
      
      I_PRINT = 0

      IIDEM = Int(DTI / TDEM + 0.1)
      
      IF(BCDMOD .NE. 'MOD2') CALL NORMAL   !计算床底和自由表面法向量
                
      Do I = 1 , IIDEM  
          
        NNDEM = NNDEM + 1
     
        TTTDEM = TTTDEM + TDEM 
        
        Print*, '     NUM. OF DEM.TSTEP=',I     
        
        BEGIN1 = OMP_GET_WTIME()
        
        BEGIN2 = OMP_GET_WTIME()
        If( I_IPKPDEM .EQ. 1 .AND. 
     &     (I .EQ. 1 .OR. PSMOD .EQ. 'MOD1') ) THEN
          Call IPPKPPDEM   ! 颗粒流场网格寻址   
        ENDIF
        END2 = OMP_GET_WTIME()
        
        BEGIN3 = OMP_GET_WTIME()
        IF(I_IPJPKPDEM .EQ. 1)THEN           
          IF(PSMOD .EQ. 'MOD2' )THEN
            CALL IPJPKPDEM1 ! 颗粒背景网格寻址(颗粒地址改变存储)
          ELSEIF(PSMOD .EQ. 'MOD3' )THEN
            CALL IPJPKPDEM2 ! 颗粒背景网格寻址(颗粒地址每步存储) 
          ENDIF
        ENDIF          
        END3 = OMP_GET_WTIME()
        
        BEGIN4 = OMP_GET_WTIME()
        IF(I_FCOUPLING .EQ. 1 .AND. I .EQ. 1) Call FCOUPLING   ! 流体颗粒耦合作用力计算 
        !Call FCOUPLINGoriginal
        END4 = OMP_GET_WTIME()       
     
        BEGIN5 = OMP_GET_WTIME()
        If(I_BCDEM .EQ. 1) THEN
            IF(BCDMOD .EQ. 'MOD2') THEN
                Call BCDEM1   ! 简单边界接触识别及虚拟颗粒生成
            ELSEIF(BCDMOD .EQ. 'MOD1') THEN
                Call BCDEM   ! 复杂边界接触识别及虚拟颗粒生成
            ENDIF
        ENDIF
        END5 = OMP_GET_WTIME()
        
        !If(DEM .EQ. 'DEMHARD') Call HARDMODULE   ! 硬球模型        
        
        BEGIN6 = OMP_GET_WTIME()
        If(DEM .EQ. 'DEMSOFT' .AND. I_SOFTMODULE .EQ. 1) THEN   ! 颗粒间及颗粒与边界碰撞接触计算
            IF(PSMOD .EQ. 'MOD1') Call SOFTMODULE1     
            IF(PSMOD .EQ. 'MOD2') Call SOFTMODULE2
            IF(PSMOD .EQ. 'NONE') Call SOFTMODULE0
        ENDIF
        !If(DEM .EQ. 'DEMSOFT') Call SOFTMODULEoriginal   ! soft sphere model  
        END6 = OMP_GET_WTIME()

        BEGIN7 = OMP_GET_WTIME()
        IF(I_MOVE .EQ. 1 .OR. I_FCOUPLING .EQ. 1) Call PMOVE   ! 颗粒平移和旋转运动
        END7 = OMP_GET_WTIME()

        BEGIN8 = OMP_GET_WTIME()
        IF(I_ACHIVEDEM .EQ. 1) Call ACHIVEDEM ! DEM结果保存
        END8 = OMP_GET_WTIME()
        
        END1 = OMP_GET_WTIME()
        
        IF(NNDEM .LE. 1E6 .AND. I_TIMESAVE .EQ. 1)THEN
            CPUTDEM(NNDEM,1) = END1-BEGIN1   !总时间
            CPUTDEM(NNDEM,2) = END2-BEGIN2   !流场寻址时间
            CPUTDEM(NNDEM,3) = END3-BEGIN3   !背景寻址时间
            CPUTDEM(NNDEM,4) = END4-BEGIN4   !耦合作用时间
            CPUTDEM(NNDEM,5) = END5-BEGIN5   !边界计算时间
            CPUTDEM(NNDEM,6) = END6-BEGIN6   !碰撞计算时间
            CPUTDEM(NNDEM,7) = END7-BEGIN7   !运动计算时间
            CPUTDEM(NNDEM,8) = END8-BEGIN8   !结果保存时间 
            Write(1121,'(I14,8F14.8)')NNDEM,(CPUTDEM(I,K),K=1,8)
        ENDIF
        
        IF(I_PRINT .EQ. 1) THEN
            Print*,'     TOTAL TIME: ',END1-BEGIN1
            Print*,' IPPKPPDEM TIME: ',END2-BEGIN2
            Print*,' IPJPKPDEM TIME: ',END3-BEGIN3
            Print*,' FCOUPLING TIME: ',END4-BEGIN4 
            Print*,'     BCDEM TIME: ',END5-BEGIN5
            Print*,'SOFTMODULE TIME: ',END6-BEGIN6
            Print*,'     PMOVE TIME: ',END7-BEGIN7
            Print*,' ACHIVEDEM TIME: ',END8-BEGIN8
        ENDIF         
                    
       IF(Mod(NNDEM,NREAD) .EQ. 0)THEN   !重新读取DEM_SETTING文件
        
        IF(I_READ_DEMSET .NE. 3)THEN   
          
         PRINT*,'INPUT I_READ_DEMSET:'
         PRINT*,'1: READDING'
         PRINT*,'2: NOT READDING'
         PRINT*,'3: CANCEL READDING'
         READ(*,*)I_READ_DEMSET
        
         IF(I_READ_DEMSET .EQ. 1)THEN
          
          CALL DEMSETREADING
          
         ENDIF
        
        ENDIF
      
        IF(I_READ_PARTICLE_SOFT .NE. 3)THEN
          
         PRINT*,'INPUT I_READ_PARTICLE_SOFT:'
         PRINT*,'1: READDING'
         PRINT*,'2: NOT READDING'
         PRINT*,'3: CANCEL READDING'
         READ(*,*)I_READ_PARTICLE_SOFT
        
         IF(I_READ_PARTICLE_SOFT .EQ. 1)THEN
          
          CALL PARTICLE_INFO

          CALL IPPKPPINITIAL
          
         ENDIF
        
        ENDIF
       
       !PRINT*,'INPUT NREAD:'
       !READ(*,*)NREAD
       
       ENDIF      
      
      Enddo      
      
      Return
      End


