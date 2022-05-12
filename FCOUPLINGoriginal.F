c-----------------------------------------------------------------------------------c
c                   SUBROUTINE PROGRAM                                              #
C                 VERSION 1.0 (13/04/2015)                                          #
C                 AUTHORIZED BY CHE HAIOU                                           #
C                               SHANGHAI JIAO TONG UNIVERSITY                       #
C                               SHANGHAI, CHINA                                     #
c-----------------------------------------------------------------------------------#
c                       DEM Force  between  particle  and  fluid                    #
c-----------------------------------------------------------------------------------c
      Subroutine FCOUPLINGoriginal
      Include './Include/OCERM_INF'   
      Double Precision UFPDEM,VFPDEM,WFPDEM
      Double Precision DFXDEM,DFYDEM,DFZDEM
      Double Precision ReDEM,CDDEM
      Integer N
          
      Dimension IC(20000,2), N_SYM(IJM) 
      Dimension KK(2) 
      INTEGER I, K, L, M, N_E, N_COUNT1, N_COUNT2, N_LOOP1,N_LOOP2
      INTEGER N_NEIGHBOUR0, N_NEIGHBOUR, N_START, N_END, LMAX
      INTEGER I_CELL, I_N1, I_N2, I_N3, I_N4 
      Double Precision DISCRI, DISCRI2, DISCRI3
      Double Precision Z_N1, Z_N2, Z_N3, Z_N4, Z_N5, Z_N6, Z_N7, Z_N8
      Double Precision DISN1, DISN2, DISN3, DISN4
      Double Precision DISN5, DISN6, DISN7, DISN8              
      
      Dimension W_GAUSS(10000), VP_MESH(IJM,KBM) 
      Dimension FPFX_DEM(IJM,KBM), FPFY_DEM(IJM,KBM), FPFZ_DEM(IJM,KBM)
      Double Precision SIGMA_GAUSS, V_MESH, Z_CELL, D_IJ, SUM_WGAUSS
      Double Precision SUMU_WGAUSS, SUMV_WGAUSS, SUMW_WGAUSS
      Double Precision FUN_PORE, CONS_TCHI
        
      Double Precision SUMPGX_WGAUSS,SUMPGY_WGAUSS,SUMPGZ_WGAUSS
      Double Precision PG_X,PG_Y,PG_Z 
      
      Dimension W_GAUSS2(10000),W_GAUSS3(10000)
      Double Precision SIGMA_GAUSS2,SUM_WGAUSS2
      Double Precision SIGMA_GAUSS3, SUM_WGAUSS3
      
      PARAMETER(N_SCHEME = 1)  !identifing neighbour cells of particle. (1:centrol, 2: vertix)
      PARAMETER(N_SYM_PG = 1)  !pressure gradient force(0:neglect, 1:include)  
         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(K,J)      
      DO K=1,KBM
!$OMP DO          
          DO J=1,IJM
              VP_MESH(J,K) = 0.0
              FPFX_DEM(J,K) = 0.0
              FPFY_DEM(J,K) = 0.0
              FPFZ_DEM(J,K) = 0.0
          ENDDO
!$OMP END DO          
      ENDDO 
!$OMP END PARALLEL 

!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(J)      
!!$OMP DO         
      DO J=1,10000
          W_GAUSS(J) = 0.0
          W_GAUSS2(J) = 0.0
          W_GAUSS3(J) = 0.0
      ENDDO
!!$OMP END DO           
!!$OMP END PARALLEL 
      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(IC, N_SYM, I, KK, K, L, M, 
!$OMP&         N_E, N_COUNT1, N_COUNT2, N_LOOP1, N_LOOP2, 
!$OMP&         N_NEIGHBOUR0, N_NEIGHBOUR, N_START, N_END, 
!$OMP&         LMAX, I_CELL, I_N1, I_N2, I_N3, I_N4, DISCRI, DISCRI2, DISCRI3 
!$OMP&         Z_N1, Z_N2, Z_N3, Z_N4, Z_N5, Z_N6, Z_N7, Z_N8, 
!$OMP&         DISN1, DISN2, DISN3, DISN4, DISN5, DISN6, DISN7, DISN8,         
!$OMP&         W_GAUSS, SIGMA_GAUSS, V_MESH, Z_CELL, D_IJ, SUM_WGAUSS, 
!$OMP&         SUMU_WGAUSS, SUMV_WGAUSS, SUMW_WGAUSS, UFPDEM, VFPDEM, WFPDEM,
!$OMP&         ReDEM, CDDEM, FUN_PORE, CONS_TCHI, SUMPGX_WGAUSS, 
!$OMP&         SUMPGY_WGAUSS, SUMPGZ_WGAUSS, PG_X, PG_Y, PG_Z, 
!$OMP&         W_GAUSS2, SIGMA_GAUSS2, SUM_WGAUSS2, 
!$OMP&         W_GAUSS3, SIGMA_GAUSS3, SUM_WGAUSS3) 
!$OMP DO 

      DO I = 1, I_PARTICLE
      KK(1) = -1
      KK(2) = 1   
      IF(KPP(I) .NE. 0)THEN    

c----------------------------------------------------------------------------c
C             identifing  neighbour cells                                    c
c----------------------------------------------------------------------------c           
      DISCRI = R_ED1 * 2 * RPDEM(I)        
      SIGMA_GAUSS = R_SIGMA1 * 2 * RPDEM(I)       
        
      DISCRI2 = R_ED2 * 2 * RPDEM(I)  
      SIGMA_GAUSS2 = R_SIGMA2 * 2 * RPDEM(I)  
                
      DISCRI3 = R_ED3 * 2 * RPDEM(I)  
      SIGMA_GAUSS3 = R_SIGMA3 * 2 * RPDEM(I)  
          
      IC(1,1) = IPP(I)
      IC(1,2) = KPP(I)
        
      N_COUNT2 = 1

      N_LOOP1 = 0
      N_LOOP2 = 0 

      N_NEIGHBOUR = 1     
              
      DO WHILE(N_COUNT2 .NE. 0)                               
  
      N_NEIGHBOUR0 = N_NEIGHBOUR

      N_LOOP2 = N_LOOP2 + 1

      LMAX = MIN(N_LOOP2, 2)        

      DO L = 1, LMAX          
  
      K = KPP(I) + (N_LOOP2 - 1) * KK(L)
 
      IF((K .GE. 1) .AND. (K .LE. KBM))THEN

      N_COUNT1 = 1
      N_START = 1
      N_END = 1
      
      DO J = 1, IJM
          N_SYM(J) = 0
      ENDDO
      
      IF(LMAX .EQ. 1)THEN
          N_SYM(IPP(I)) = 1
      ENDIF      
 
      DO WHILE(N_COUNT1 .NE. 0)
  
        N_LOOP1 = N_LOOP1+1
        N_COUNT1 = 0

        DO M = N_START, N_END 

          DO N_E = 1, CELL_POLYGEN(IC(M, 1))            
            
            I_CELL = CELL_SIDE(IC(M, 1), N_E, 2)

            IF (I_CELL .NE. -999)THEN
              
            IF(N_SCHEME .EQ. 1)THEN   
                
            Z_CELL = ZZ(K) * DC(I_CELL) + ELF(I_CELL)               
            D_IJ=SQRT( (XPDEM(I)-CXY(I_CELL,1))*(XPDEM(I)-CXY(I_CELL,1))
     &           + (YPDEM(I)-CXY(I_CELL,2))*(YPDEM(I)-CXY(I_CELL,2))
     &           + (ZPDEM(I)-Z_CELL)*(ZPDEM(I)-Z_CELL) )    

            IF((D_IJ .LT. DISCRI) .AND. (N_SYM(I_CELL) .EQ. 0))THEN
                N_COUNT1 = N_COUNT1 + 1              
                N_NEIGHBOUR = N_NEIGHBOUR + 1
                IC(N_NEIGHBOUR,1) = I_CELL
                IC(N_NEIGHBOUR,2) = K 
                N_SYM(I_CELL) = 1                                 
            ENDIF            
            
            ELSEIF(N_SCHEME .EQ. 2)THEN
                
            I_N1 = CELL_INDEX(I_CELL, 1)
            I_N2 = CELL_INDEX(I_CELL, 2)
            I_N3 = CELL_INDEX(I_CELL, 3)
            I_N4 = CELL_INDEX(I_CELL, CELL_POLYGEN(I_CELL))
            Z_N1 = Z(K) * (HP(I_N1) + ELFV(I_N1)) + ELFV(I_N1)
            Z_N2 = Z(K) * (HP(I_N2) + ELFV(I_N2)) + ELFV(I_N2)
            Z_N3 = Z(K) * (HP(I_N3) + ELFV(I_N3)) + ELFV(I_N3)
            Z_N4 = Z(K) * (HP(I_N4) + ELFV(I_N4)) + ELFV(I_N4)
            Z_N5 = Z(K+1) * (HP(I_N1) + ELFV(I_N1)) + ELFV(I_N1)
            Z_N6 = Z(K+1) * (HP(I_N2) + ELFV(I_N2)) + ELFV(I_N2)
            Z_N7 = Z(K+1) * (HP(I_N3) + ELFV(I_N3)) + ELFV(I_N3)
            Z_N8 = Z(K+1) * (HP(I_N4) + ELFV(I_N4)) + ELFV(I_N4)        

            DISN1=SQRT((XPDEM(I)-PXY(I_N1,1)) * (XPDEM(I)-PXY(I_N1,1))
     &             + (YPDEM(I)-PXY(I_N1,2)) * (YPDEM(I)-PXY(I_N1,2))
     &             + (ZPDEM(I)-Z_N1) * (ZPDEM(I)-Z_N1))   
            DISN2=SQRT((XPDEM(I)-PXY(I_N2,1)) * (XPDEM(I)-PXY(I_N2,1))
     &             + (YPDEM(I)-PXY(I_N2,2)) * (YPDEM(I)-PXY(I_N2,2))
     &             + (ZPDEM(I)-Z_N2) * (ZPDEM(I)-Z_N2)) 
            DISN3=SQRT((XPDEM(I)-PXY(I_N3,1)) * (XPDEM(I)-PXY(I_N3,1))
     &             + (YPDEM(I)-PXY(I_N3,2)) * (YPDEM(I)-PXY(I_N3,2))
     &             + (ZPDEM(I)-Z_N3) * (ZPDEM(I)-Z_N3))
            DISN4=SQRT((XPDEM(I)-PXY(I_N4,1)) * (XPDEM(I)-PXY(I_N4,1))
     &             + (YPDEM(I)-PXY(I_N4,2)) * (YPDEM(I)-PXY(I_N4,2))
     &             + (ZPDEM(I)-Z_N4) * (ZPDEM(I)-Z_N4))
            DISN5=SQRT((XPDEM(I)-PXY(I_N1,1)) * (XPDEM(I)-PXY(I_N1,1))
     &             + (YPDEM(I)-PXY(I_N1,2)) * (YPDEM(I)-PXY(I_N1,2))
     &             + (ZPDEM(I)-Z_N5) * (ZPDEM(I)-Z_N5))   
            DISN6=SQRT((XPDEM(I)-PXY(I_N2,1)) * (XPDEM(I)-PXY(I_N2,1))
     &             + (YPDEM(I)-PXY(I_N2,2)) * (YPDEM(I)-PXY(I_N2,2))
     &             + (ZPDEM(I)-Z_N6) * (ZPDEM(I)-Z_N6)) 
            DISN7=SQRT((XPDEM(I)-PXY(I_N3,1)) * (XPDEM(I)-PXY(I_N3,1))
     &             + (YPDEM(I)-PXY(I_N3,2)) * (YPDEM(I)-PXY(I_N3,2))
     &             + (ZPDEM(I)-Z_N7) * (ZPDEM(I)-Z_N7))
            DISN8=SQRT((XPDEM(I)-PXY(I_N4,1)) * (XPDEM(I)-PXY(I_N4,1))
     &             + (YPDEM(I)-PXY(I_N4,2)) * (YPDEM(I)-PXY(I_N4,2))
     &             + (ZPDEM(I)-Z_N8) * (ZPDEM(I)-Z_N8))            

            IF(((DISN1 .LT. DISCRI) .OR. (DISN2 .LT. DISCRI)
     &        .OR. (DISN3 .LT. DISCRI) .OR. (DISN4 .LT. DISCRI)
     &        .OR. (DISN5 .LT. DISCRI) .OR. (DISN6 .LT. DISCRI)
     &        .OR. (DISN7 .LT. DISCRI) .OR. (DISN8 .LT. DISCRI))    
     &        .AND. (N_SYM(I_CELL) .EQ. 0))THEN 
                  N_COUNT1 = N_COUNT1 + 1              
                  N_NEIGHBOUR = N_NEIGHBOUR + 1
                  IC(N_NEIGHBOUR,1) = I_CELL
                  IC(N_NEIGHBOUR,2) = K 
                  N_SYM(I_CELL) = 1                                 
            ENDIF 
            
            ELSE
                print*,'please check the value of N_SCHEME!'
            ENDIF 
                                    
            ENDIF 
          ENDDO                       
        ENDDO 
         
        N_START=N_END+1
        N_END=N_END+N_COUNT1   
      
      ENDDO 
      ENDIF        
      ENDDO 
      
      N_COUNT2 = N_NEIGHBOUR - N_NEIGHBOUR0
      
      ENDDO 

      IF(N_NEIGHBOUR .GT. 10000)THEN
          PRINT*,'neighbour cell numbers is greater than 10000'
          PAUSE
      ENDIF     
      
      SUM_WGAUSS = 0.0
      SUMU_WGAUSS = 0.0
      SUMV_WGAUSS = 0.0
      SUMW_WGAUSS = 0.0
      SUMPGX_WGAUSS = 0.0
      SUMPGY_WGAUSS = 0.0
      SUMPGZ_WGAUSS = 0.0
      SUM_WGAUSS2 = 0.0
      SUM_WGAUSS3 = 0.0
      DO J = 1, N_NEIGHBOUR
          V_MESH = AREA(IC(J,1))*DZ(IC(J,2))*DC(IC(J,1))
          Z_CELL = ZZ(IC(J,2)) * DC(IC(J,1)) + ELF(IC(J,1))           
          D_IJ=SQRT((XPDEM(I)-CXY(IC(J,1),1))*(XPDEM(I)-CXY(IC(J,1),1))
     &           + (YPDEM(I)-CXY(IC(J,1),2))*(YPDEM(I)-CXY(IC(J,1),2))
     &           + (ZPDEM(I)-Z_CELL)*(ZPDEM(I)-Z_CELL))
          
          W_GAUSS(J)=1.0/(SQRT(2.0 * ACOS(-1.0)) * SIGMA_GAUSS) * 
     &        EXP(-D_IJ*D_IJ/(2.0*SIGMA_GAUSS*SIGMA_GAUSS))*V_MESH
          SUM_WGAUSS=SUM_WGAUSS+W_GAUSS(J)
          SUMU_WGAUSS=SUMU_WGAUSS+UR(IC(J,1),IC(J,2))*W_GAUSS(J)
          SUMV_WGAUSS=SUMV_WGAUSS+VR(IC(J,1),IC(J,2))*W_GAUSS(J)
          SUMW_WGAUSS=SUMW_WGAUSS+WR(IC(J,1),IC(J,2))*W_GAUSS(J)        
          
          SUMPGX_WGAUSS=SUMPGX_WGAUSS+
     &                  GRADP(IC(J,1),IC(J,2),1)*W_GAUSS(J)
          SUMPGY_WGAUSS=SUMPGY_WGAUSS+
     &                  GRADP(IC(J,1),IC(J,2),2)*W_GAUSS(J)
          SUMPGZ_WGAUSS=SUMPGZ_WGAUSS+
     &                  GRADP(IC(J,1),IC(J,2),3)*W_GAUSS(J)
          
          IF( (D_IJ .LT. DISCRI2) .OR. (J .EQ. 1) )THEN
          W_GAUSS2(J)=1.0/(SQRT(2.0 * ACOS(-1.0)) * SIGMA_GAUSS2) * 
     &        EXP(-D_IJ*D_IJ/(2.0*SIGMA_GAUSS2*SIGMA_GAUSS2))*V_MESH
          SUM_WGAUSS2=SUM_WGAUSS2+W_GAUSS2(J)       
          ENDIF
          
          IF( (D_IJ .LT. DISCRI3) .OR. (J .EQ. 1) )THEN
          W_GAUSS3(J)=1.0/(SQRT(2.0 * ACOS(-1.0)) * SIGMA_GAUSS3) * 
     &        EXP(-D_IJ*D_IJ/(2.0*SIGMA_GAUSS3*SIGMA_GAUSS3))*V_MESH
          SUM_WGAUSS3=SUM_WGAUSS3+W_GAUSS3(J)       
          ENDIF
          
      ENDDO
      
      UFPDEM=SUMU_WGAUSS/SUM_WGAUSS-UPDEM(I)
      VFPDEM=SUMV_WGAUSS/SUM_WGAUSS-VPDEM(I)
      WFPDEM=SUMW_WGAUSS/SUM_WGAUSS-WPDEM(I)
      
      PG_X=SUMPGX_WGAUSS/SUM_WGAUSS
      PG_Y=SUMPGY_WGAUSS/SUM_WGAUSS
      PG_Z=SUMPGZ_WGAUSS/SUM_WGAUSS
                   
      IF(N_SYM_PG .EQ. 1 .AND. NSTEP .NE. 1) THEN
          F_PG_X(I) = - 1.00 * PG_X * V_PARTICLE(I)
          F_PG_Y(I) = - 1.00 * PG_Y * V_PARTICLE(I)
          F_PG_Z(I) = - 1.00 * PG_Z * V_PARTICLE(I)   
      ELSE
          F_PG_X(I) = 0.0
          F_PG_Y(I) = 0.0
          F_PG_Z(I) = RMEAN(IPP(I),KPP(I)) * V_PARTICLE(I) * GRAV
      ENDIF      
C--------------------------------------------- 
c---Reynolds number
      IF(I_PARTICLE .EQ. 1)THEN      
          ReDEM = (2. * RPDEM(I) ) *  
     &            Sqrt(UFPDEM**2+VFPDEM**2+WFPDEM**2) / UMOL
      ELSE
          ReDEM = (2. * RPDEM(I) )*  
     &            Sqrt(UFPDEM**2+VFPDEM**2+WFPDEM**2) / UMOL * 
     &            PORE(IPP(I),KPP(I))
      ENDIF
c------drag coefficient according to Re          
      IF(ReDEM .GT. 0.0)THEN
          !CDDEM = (0.63 + 4.8 / SQRT(ReDEM))**2  ! (by DI FELICE, 1994)
          CDDEM = 24.0/9.06**2*(9.06/SQRT(ReDEM)+1)**2  ! (by Abraham, 1994)
          !CDDEM = 24.0*(1.0+0.15*ReDEM**0.687)/ReDEM   ! (Schiller, 1935)
          CDDEM = DMAX1(CDDEM,0.44)           
      ELSE
          CDDEM = 0.0
      ENDIF     
     	
      IF(I_PARTICLE .EQ. 1)THEN
          FUN_PORE=1.0
      ELSE
          CONS_TCHI=3.7-0.65*EXP(-(1.5-LOG10(ReDEM))**2/2.0)
          FUN_PORE=PORE(IPP(I),KPP(I))**(1-CONS_TCHI)
      ENDIF      
      !FUN_PORE=1.0
	FDXDEM(I) = 0.5 * CDDEM * RMEAN(IPP(I),KPP(I)) *  
     &	  Acos(-1.) * RPDEM(I)**2 * UFPDEM * 
     &      Sqrt(UFPDEM**2+VFPDEM**2+WFPDEM**2) * FUN_PORE
C-------------------------------------------------	
	FDYDEM(I) = 0.5 * CDDEM * RMEAN(IPP(I),KPP(I)) *  
     &	  Acos(-1.) * RPDEM(I)**2 * VFPDEM * 
     &      Sqrt(UFPDEM**2+VFPDEM**2+WFPDEM**2) * FUN_PORE 
C------------------------------------------------		
	FDZDEM(I) = 0.5 * CDDEM * RMEAN(IPP(I),KPP(I)) *  
     &	  Acos(-1.) * RPDEM(I)**2 * WFPDEM * 
     &      Sqrt(UFPDEM**2+VFPDEM**2+WFPDEM**2) * FUN_PORE
C-------------------------------------------------  
      DO J = 1, N_NEIGHBOUR
        Z_CELL = ZZ(IC(J,2)) * DC(IC(J,1)) + ELF(IC(J,1))   
        D_IJ=SQRT((XPDEM(I)-CXY(IC(J,1),1))*(XPDEM(I)-CXY(IC(J,1),1))
     &           + (YPDEM(I)-CXY(IC(J,1),2))*(YPDEM(I)-CXY(IC(J,1),2))
     &           + (ZPDEM(I)-Z_CELL)*(ZPDEM(I)-Z_CELL))

        IF( (D_IJ .LT. DISCRI2) .OR. (J .EQ. 1) )THEN
        FPFX_DEM(IC(J,1),IC(J,2))=FPFX_DEM(IC(J,1),IC(J,2))+
     &            W_GAUSS2(J)/SUM_WGAUSS2*FDXDEM(I)
        FPFY_DEM(IC(J,1),IC(J,2))=FPFY_DEM(IC(J,1),IC(J,2))+
     &            W_GAUSS2(J)/SUM_WGAUSS2*FDYDEM(I)
        FPFZ_DEM(IC(J,1),IC(J,2))=FPFZ_DEM(IC(J,1),IC(J,2))+
     &            W_GAUSS2(J)/SUM_WGAUSS2*FDZDEM(I)
        ENDIF

        IF( (D_IJ .LT. DISCRI3) .OR. (J .EQ. 1) )THEN
        VP_MESH(IC(J,1),IC(J,2))=VP_MESH(IC(J,1),IC(J,2))+
     &            W_GAUSS3(J)/SUM_WGAUSS3*V_PARTICLE(I)
        ENDIF        
        
        ENDDO         
      ELSE
        FDXDEM(I) = 0.0
        FDYDEM(I) = 0.0
        FDZDEM(I) = 0.0
        F_PG_X(I) = 0.0
        F_PG_Y(I) = 0.0
        F_PG_Z(I) = 0.0
      ENDIF
      
      ENDDO           
!$OMP END DO
!$OMP END PARALLEL

c     Calculate the sum of forces that particles force on fluid (gridding)
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(J,K,V_MESH)      
      Do K = 1 , KBM
!$OMP DO        
        Do J = 1 , IJM
            V_MESH = AREA(J)*DZ(K)*DC(J)
            
            IF(POREMODULE  .EQ. 'INCLUDE' .AND. I_PTF .EQ. 1) Then
             PORE(J,K)=DMIN1(DMAX1(1.0-VP_MESH(J,K)/V_MESH,0.0001),1.0)
            ENDIF
                              
            DEMFORCEX(J,K) = -1.0/RMEAN(J,K)*DC(J)*FPFX_DEM(J,K)/V_MESH
            DEMFORCEY(J,K) = -1.0/RMEAN(J,K)*DC(J)*FPFY_DEM(J,K)/V_MESH
            DEMFORCEZ(J,K) = -1.0/RMEAN(J,K)*DC(J)*FPFZ_DEM(J,K)/V_MESH
         Enddo
!$OMP END DO                    
      Enddo
!$OMP END PARALLEL             
      Return
      End       
 