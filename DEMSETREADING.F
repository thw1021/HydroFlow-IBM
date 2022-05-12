      Subroutine DEMSETREADING 
        Include './Include/OCERM_INF'
        Dimension COM(80)
C--- Reading DEM setting        
        Open(10010,file='./Include/DEM_SETTING.DAT')   
        
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(5I14)') I_IPKPDEM,I_BCDEM,I_SOFTMODULE,
     &                       I_FCOUPLING,I_ACHIVEDEM
        Write (IUPRT,'(5I14)') I_IPKPDEM,I_BCDEM,I_SOFTMODULE,
     &                       I_FCOUPLING,I_ACHIVEDEM        
        
        Read(10010,5200) (COM(I), I= 1, 80)
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(5I14)') I_IPJPKPDEM,I_MOVE,I_TRAN,I_ROT,I_GRA
        Write (IUPRT,'(5I14)') I_IPJPKPDEM,I_MOVE,I_TRAN,I_ROT,I_GRA  
        
        Read(10010,5200) (COM(I), I= 1, 80)
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(E14.6,2I14,4X,A10,I14)') TDEM,ISDEM,IFLOW,
     &                                       DEMRESTART,IDEMRES
        Write (IUPRT,'(E14.6,2I14,4X,A10,I14)') TDEM,ISDEM,IFLOW,
     &                                         DEMRESTART,IDEMRES
        
        Read(10010,5200) (COM(I), I= 1, 80)
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(2I14,10X,A4,10X,A4)') I_PTF,I_FTP,DRAGMOD,CDMOD
        Write (IUPRT,'(2I14,10X,A4,10X,A4)') I_PTF,I_FTP,DRAGMOD,CDMOD
        
        Read(10010,5200) (COM(I), I= 1, 80)
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(10X,A4,3F14.6)') MPTMOD,CFKDEM,CFWDEM1,CFWDEM2
        Write (IUPRT,'(10X,A4,3F14.6)') MPTMOD,CFKDEM,CFWDEM1,CFWDEM2 
        
        Read(10010,5200) (COM(I), I= 1, 80)
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(10X,A4,10X,A4,12X,A2)') ROLLMOD,ROLLBOUND,CONMOD
        Write (IUPRT,'(10X,A4,10X,A4,12X,A2)') ROLLMOD,ROLLBOUND,CONMOD

        Read(10010,5200) (COM(I), I= 1, 80)
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(4E14.6)') PKN_M1PP,PKN_M1PW,EN_M1PP,EN_M1PW
        Write (IUPRT,'(4E14.6)') PKN_M1PP,PKN_M1PW,EN_M1PP,EN_M1PW  
        
        Read(10010,5200) (COM(I), I= 1, 80)
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(4E14.6)') RK_M1PP,RK_M1PW,RC_M1PP,RC_M1PW
        Write (IUPRT,'(4E14.6)') RK_M1PP,RK_M1PW,RC_M1PP,RC_M1PW      

        Read(10010,5200) (COM(I), I= 1, 80)        
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(4E14.6)') PMIUS,PMIUSPW,RMIUS,RMIUSPW
        Write (IUPRT,'(4E14.6)') PMIUS,PMIUSPW,RMIUS,RMIUSPW          

        Read(10010,5200) (COM(I), I= 1, 80)        
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(3E14.6,I14,6X,A8)') C_VIRTUAL,RHO_FLUID,
     &                         RHO_PARTICLE,I_CV_PORE,MOD_FA 
        Write (IUPRT,'(3E14.6,I14,6X,A8)') C_VIRTUAL,RHO_FLUID,
     &                         RHO_PARTICLE,I_CV_PORE,MOD_FA      
        
c        IF(I_FTP .EQ. 0) THEN
c          C_VIRTUAL = 0.0
c        ENDIF        

        Read(10010,5200) (COM(I), I= 1, 80)        
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(3E14.6)') R_ED1,R_ED2,R_ED3 
        Write (IUPRT,'(3E14.6)') R_ED1,R_ED2,R_ED3         
        
        Read(10010,5200) (COM(I), I= 1, 80)        
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(3E14.6)') R_SIGMA1,R_SIGMA2,R_SIGMA3 
        Write (IUPRT,'(3E14.6)') R_SIGMA1,R_SIGMA2,R_SIGMA3      
        
        Read(10010,5200) (COM(I), I= 1, 80)        
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(10X,A4,2I14)') PSMOD,NSH,NSV 
        Write (IUPRT,'(10X,A4,2I14)') PSMOD,NSH,NSV         
        
        Read(10010,5200) (COM(I), I= 1, 80)        
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(10X,A4,3I14)') BCDMOD,NREAD,I_TIMESAVE,IP_RANDOM
        Write (IUPRT,'(10X,A4,3I14)') BCDMOD,NREAD,I_TIMESAVE,IP_RANDOM
        
        Read(10010,5200) (COM(I), I= 1, 80)        
	  Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(I14)') N_CYLINDER
        Write (IUPRT,'(I14)') N_CYLINDER        

        Read(10010,5200) (COM(I), I= 1, 80)        
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)
        Read(10010,'(4E14.6)') X_CYLINDER(I),Y_CYLINDER(I),
     &                Z_CYLINDER(I),R_CYLINDER(I) 
        Write (IUPRT,'(4E14.6)') X_CYLINDER(I),Y_CYLINDER(I),
     &                Z_CYLINDER(I),R_CYLINDER(I) 
              
        IF(N_CYLINDER .GT. 1)THEN
          DO I = 1, N_CYLINDER
              Read(10010,'(4E14.6)') X_CYLINDER(I),Y_CYLINDER(I),
     &                      Z_CYLINDER(I),R_CYLINDER(I)   
              Write (IUPRT,'(4E14.6)') X_CYLINDER(I),Y_CYLINDER(I),
     &                      Z_CYLINDER(I),R_CYLINDER(I) 
          ENDDO
        ENDIF                
        
        Read(10010,5200) (COM(I), I= 1, 80)
        Read(10010,5200) (COM(I), I= 1, 80)
        Write (IUPRT,5200) (COM(I), I= 1, 80)        
        
        Close(10010)      
        
      Return
      
5200  Format (80A1)      
      
      End
