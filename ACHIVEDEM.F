c-----------------------------------------------------------------------------------c
c                             SUBROUTINE PROGRAM                                    #
C                             DESINGED BY ChenJun                                   #
C                                                 SHANGHAI JIAO TONG UNIVERSITY     #
C                                                 SHANGHAI, CHINA                   #
c-----------------------------------------------------------------------------------#
c 子程序功能：
c 1------DEM结果保存      
c-----------------------------------------------------------------------------------#
      Subroutine ACHIVEDEM
      Include './Include/OCERM_INF'
      Character*10 FN
      !Integer NDEM
      
      If ((Mod(NNDEM,ISDEM) .EQ. 0 ) .OR. (NNDEM .EQ. 1)) Then
C        If(Mod(NNDEM,10*ISDEM) .EQ. -1 )Then
C           Write(IDEM,'(27A16)') 
C     &        'XPDEM','YPDEM','ZPDEM',
C     &        'UPDEM','VPDEM','WPDEM',
C     &        'OMXDEM','OMYDEM','OMZDEM',          
C     &        'FNXDEM','FNYDEM','FNZDEM',
C     &        'FTXDEM','FTYDEM','FTZDEM' 
C        Endif
        
        Write (IDEM,*) 'ZONE T=', '"',TTTDEM,'"' 
      
        Do I = 1,I_PARTICLE	
         IF(NSYM_CONT(I) .EQ. 1)THEN   
          Write(IDEM,'(30E25.16)')  
     &        XPDEM(I),YPDEM(I),ZPDEM(I),
     &        UPDEM(I),VPDEM(I),WPDEM(I),
     &        OMXDEM(I),OMYDEM(I),OMZDEM(I), 
     &        FDXDEM(I),FDYDEM(I),FDZDEM(I),
     &        FAXDEM(I),FAYDEM(I),FAZDEM(I),        
     &        F_PG_X(I),F_PG_Y(I),F_PG_Z(I),          
     &        FNXDEM(I),FNYDEM(I),FNZDEM(I),
     &        FTXDEM(I),FTYDEM(I),FTZDEM(I),RPDEM(I)          
!     &        DELTAT_PP1X(I,1),DELTAT_PP1Y(I,1),DELTAT_PP1Z(I,1),     
!     &        DELTAT_PP1X(I,2),DELTAT_PP1Y(I,2),DELTAT_PP1Z(I,2)      
         ENDIF
        Enddo  

C        WRITE(FN,'(F8.4)') DFLOAT(NNDEM)*TDEM
C        Open(1235,File='./TEMPSAVE/'//FN//'DEM.DAT',
C     &	        Status='unknown')
             
C        Write (1235,*) 'TITLE="DEM results of particles" '
C        Write (1235,*) 'VARIBLES= "X","Y","Z","U","V","W" '   
C        Write (1235,*) 'ZONE T=', '"',TTTDEM,'"'
C        Do I = 1,I_PARTICLE	
C          IF(NSYM_CONT(I) .EQ. 1)THEN  
C          Write(1235,'(30E25.16)')  
C     &        XPDEM(I),YPDEM(I),ZPDEM(I),
C     &        UPDEM(I),VPDEM(I),WPDEM(I),
C     &        OMXDEM(I),OMYDEM(I),OMZDEM(I), 
C     &        FDXDEM(I),FDYDEM(I),FDZDEM(I),
C     &        FAXDEM(I),FAYDEM(I),FAZDEM(I),            
C     &        F_PG_X(I),F_PG_Y(I),F_PG_Z(I),          
C     &        FNXDEM(I),FNYDEM(I),FNZDEM(I),
C     &        FTXDEM(I),FTYDEM(I),FTZDEM(I),RPDEM(I)          
C!     &        DELTAT_PP1X(I,1),DELTAT_PP1Y(I,1),DELTAT_PP1Z(I,1),     
C!     &        DELTAT_PP1X(I,2),DELTAT_PP1Y(I,2),DELTAT_PP1Z(I,2)  
C          ENDIF
C        Enddo
        
C       CLOSE(1235)
       
       IF(NNDEM .LE. 1E6 .AND. NNDEM .GT. 1)THEN
         Write(1131,'(I14,8F14.8)')NNDEM-1,(CPUTDEM(NNDEM-1,K),K=1,8)
       ENDIF
       
	Endif  
    
c      N_POSTDEM = 0
c      !PRINT*,IDEM
c      IF(N_POSTDEM .EQ. 1 .AND. (Mod(NNDEM,10) .EQ. 0 )) THEN
c        Do I = 1,I_PARTICLE	                             
c          Write(IDEM2,'(25F16.8)') TTTDEM, 
c     &        XPDEM(I),YPDEM(I),ZPDEM(I),
c     &        UPDEM(I),VPDEM(I),WPDEM(I),
c     &        OMXDEM(I),OMYDEM(I),OMZDEM(I), 
c     &        FDXDEM(I),FDYDEM(I),FDZDEM(I),
c     &        F_PG_X(I),F_PG_Y(I),F_PG_Z(I),          
c     &        FNXDEM(I),FNYDEM(I),FNZDEM(I),
c     &        FTXDEM(I),FTYDEM(I),FTZDEM(I),
c     &        UPDEM(I)-OMYDEM(I)*RPDEM(I),TPYDEM(I)
c!     &        DELTAT_PP1X(I,1),DELTAT_PP1Y(I,1),DELTAT_PP1Z(I,1),     
c!     &        DELTAT_PP1X(I,2),DELTAT_PP1Y(I,2),DELTAT_PP1Z(I,2)      
c        Enddo
c      ENDIF
        !close(10)
        
c===================output the particle information=======================c
      I_POUT = 0
      IF(Mod(NNDEM,100000) .EQ. 0 .AND. I_POUT .EQ. 1)THEN
        OPEN(24,FILE='./Results/PARTICLE_SOFT-SELECT.DAT',
     &          STATUS='UNKNOWN') 
        OPEN(25,FILE='./Results/PARTICLE_SOFT-TOTAL.DAT',
     &          STATUS='UNKNOWN') 
        
   !     NTEMP1 = 0
   !     NTEMP2 = 0
   !     Do I = 1 , I_PARTICLE
   !       IF(NSYM_CONT(I) .EQ. 1)THEN  
   !           
   !        IF(ZPDEM(I) .LE. -0.04-RPDEM(I))THEN
		 ! NTEMP1 = NTEMP1 + 1
   !        ENDIF
   !
		 !NTEMP2 = NTEMP2 + 1
   !        
   !       ENDIF
   !     Enddo  
        
	  WRITE(24,*)
        WRITE(25,*)
        Do I = 1 , I_PARTICLE
          IF(NSYM_CONT(I) .EQ. 1)THEN             
           IF(ZPDEM(I) .LE. -0.04-RPDEM(I))THEN
		 WRITE(24,'(15E28.20)') XPDEM(I),YPDEM(I),ZPDEM(I),UPDEM(I),
     &        VPDEM(I),WPDEM(I),OMXDEM(I),OMYDEM(I),OMZDEM(I),RPMEAN(I),
     &		PMDEM(I),RPDEM(I),PNIUPDEM(I),EPDEM(I),GPDEM(I) 
           ENDIF

		 WRITE(25,'(15E28.20)') XPDEM(I),YPDEM(I),ZPDEM(I),UPDEM(I),
     &        VPDEM(I),WPDEM(I),OMXDEM(I),OMYDEM(I),OMZDEM(I),RPMEAN(I),
     &		PMDEM(I),RPDEM(I),PNIUPDEM(I),EPDEM(I),GPDEM(I) 
          ENDIF
        Enddo 
        close(24)
        close(25)
      ENDIF       
c=========================================================================c        
        	            
      Return
      End     
