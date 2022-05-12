c-----------------------------------------------------------------------------------c
c                             SUBROUTINE PROGRAM                                    #
C                             DESINGED BY ChenJun                                   #
C                                                 SHANGHAI JIAO TONG UNIVERSITY     #
C                                                 SHANGHAI, CHINA                   #
c-----------------------------------------------------------------------------------#
c 子程序功能：
c 1------流体颗粒耦合作用计算      
c-----------------------------------------------------------------------------------#
      Subroutine FCOUPLING
      Include './Include/OCERM_INF'  
      
      Double Precision UFPDEM,VFPDEM,WFPDEM,URFPDEM
      Double Precision DFXDEM,DFYDEM,DFZDEM
      Double Precision ReDEM,CDDEM
      Integer N
          
      Dimension IC(30000,2) 
      Dimension KK(2) 
      INTEGER I, K, L, M, N_COUNT1, N_COUNT2, N_LOOP1
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
      
      Double Precision TEMP
      
      Dimension W_GAUSS2(30000),W_GAUSS3(30000)
      Double Precision SIGMA_GAUSS2,SUM_WGAUSS2
      Double Precision SIGMA_GAUSS3, SUM_WGAUSS3
      
      Double Precision SUMPORE_WGAUSS,POROSITY,ALPHA_SOLID,DENSITY_F
      Double Precision SUMTKE_WGAUSS,SUMTDISS_WGAUSS,TKE_P,TDISS_P
      Double Precision BELTA_DRAG
      
      PARAMETER(N_SCHEME = 1)  !identifing neighbour cells of particle. (1:centrol, 2: vertix)
      PARAMETER(N_SYM_PG = 1)  !pressure gradient force(0:neglect, 1:include)  
      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(IC, I, J, KK, K, L, M, 
!$OMP&         N_COUNT1, N_COUNT2, N_LOOP1, 
!$OMP&         N_NEIGHBOUR0, N_NEIGHBOUR, N_START, N_END, 
!$OMP&         LMAX, I_CELL, I_N1, I_N2, I_N3, I_N4, 
!$OMP&         DISCRI, DISCRI1, DISCRI2, DISCRI3 
!$OMP&         Z_N1, Z_N2, Z_N3, Z_N4, Z_N5, Z_N6, Z_N7, Z_N8, 
!$OMP&         DISN1, DISN2, DISN3, DISN4, DISN5, DISN6, DISN7, DISN8,         
!$OMP&         W_GAUSS, SIGMA_GAUSS, V_MESH, Z_CELL, D_IJ, SUM_WGAUSS, 
!$OMP&         SUMU_WGAUSS, SUMV_WGAUSS, SUMW_WGAUSS, 
!$OMP&         UFPDEM, VFPDEM, WFPDEM,URFPDEM,TEMP,
!$OMP&         ReDEM, CDDEM, FUN_PORE, CONS_TCHI, SUMPGX_WGAUSS, 
!$OMP&         SUMPGY_WGAUSS, SUMPGZ_WGAUSS, PG_X, PG_Y, PG_Z, 
!$OMP&         W_GAUSS2, SIGMA_GAUSS2, SUM_WGAUSS2, 
!$OMP&         W_GAUSS3, SIGMA_GAUSS3, SUM_WGAUSS3,
!$OMP&         SUMPORE_WGAUSS,POROSITY,ALPHA_SOLID,DENSITY_F,
!$OMP&         SUMTKE_WGAUSS,SUMTDISS_WGAUSS,TKE_P,TDISS_P,
!$OMP&         BELTA_DRAG,AUFPDEMX,AUFPDEMY,AUFPDEMZ,AUFPDEM,
!$OMP&         AUFDEMX,AUFDEMY,AUFDEMZ,AC,
!$OMP&         FTFX,FTFY,FTFZ)

      DO K=1,KBM
!$OMP DO          
          DO J=1,IJM
              VP_MESH(J,K) = 0.0
              FPFX_DEM(J,K) = 0.0
              FPFY_DEM(J,K) = 0.0
              FPFZ_DEM(J,K) = 0.0              
              STKEDEM(J,K) = 0.0
              STDISSDEM(J,K) = 0.0
          ENDDO
!$OMP END DO          
      ENDDO

!$OMP BARRIER

C!$OMP DO         
C      Do J = 1 , I_PARTICLE   !1个流场时步内N个DEM时间步颗粒速度平均 
C        UPAVEDEM(J) = UPSUMDEM(J) / DFLOAT(IIDEM)
C        VPAVEDEM(J) = VPSUMDEM(J) / DFLOAT(IIDEM)
C        WPAVEDEM(J) = WPSUMDEM(J) / DFLOAT(IIDEM)
        
C        UPSUMDEM(J) = 0.0   !颗粒速度累加值归零
C        VPSUMDEM(J) = 0.0
C        WPSUMDEM(J) = 0.0 
C      ENDDO     
C!$OMP END DO       

!$OMP DO 

      DO I = 1, I_PARTICLE
       
       IF(NSYM_CONT(I) .EQ. 1)THEN              

       UPAVEDEM(I) = UPSUMDEM(I) / DFLOAT(IIDEM)   !1个流场时步内N个DEM时间步颗粒速度平均
       VPAVEDEM(I) = VPSUMDEM(I) / DFLOAT(IIDEM)
       WPAVEDEM(I) = WPSUMDEM(I) / DFLOAT(IIDEM)
        
       UPSUMDEM(I) = 0.0   !颗粒速度累加值归零
       VPSUMDEM(I) = 0.0
       WPSUMDEM(I) = 0.0            
           
       DENSITY_F = RMEAN(IPP(I),KPP(I))
          
       KK(1) = -1
       KK(2) = 1   
       
       IF(KPP(I) .NE. 0)THEN    

c----------------------------------------------------------------------------c
C             identifing  neighbour cells                                    c
c----------------------------------------------------------------------------c           
        DISCRI1 = R_ED1 * DPDEM(I)       
        SIGMA_GAUSS = R_SIGMA1 * DPDEM(I)       
        
        DISCRI2 = R_ED2 * DPDEM(I)  
        SIGMA_GAUSS2 = R_SIGMA2 * DPDEM(I) 
                
        DISCRI3 = R_ED3 * DPDEM(I)  
        SIGMA_GAUSS3 = R_SIGMA3 * DPDEM(I)      
        
        DISCRI = DMAX1(DISCRI1,DMAX1(DISCRI2,DISCRI3))
        
        N_COUNT1 = 1

        N_LOOP1 = 0 

        N_NEIGHBOUR = 0     
              
        DO WHILE(N_COUNT1 .NE. 0)
            
         N_NEIGHBOUR0 = N_NEIGHBOUR

         N_LOOP1 = N_LOOP1 + 1

         LMAX = MIN(N_LOOP1, 2)        

         DO L = 1, LMAX          
  
          K = KPP(I) + (N_LOOP1 - 1) * KK(L)
 
          IF((K .GE. 1) .AND. (K .LE. KBM))THEN            
        
           DO J = 1, NL
            
            N_COUNT2 = 0   
               
            IF(J .EQ. 1)THEN
             N_START = 1
             N_END = 1
            ELSE
             N_START = N_END + 1
             !N_END = N_END + INL(IPP(I),J)  
             N_END = INL(IPP(I),J)
            ENDIF
            !print*,N_START,N_END
            DO M = N_START, N_END
             !print*,J,INL(IPP(I),J),M,  INE(IPP(I),M)
             I_CELL = INE(IPP(I),M)  
              
             IF(N_SCHEME .EQ. 1)THEN   
                
              Z_CELL = ZZ(K) * DC(I_CELL) + ELF(I_CELL)               
              D_IJ=SQRT( (XPDEM(I)-CXY(I_CELL,1)) *
     &                   (XPDEM(I)-CXY(I_CELL,1)) +
     &                   (YPDEM(I)-CXY(I_CELL,2)) * 
     &                   (YPDEM(I)-CXY(I_CELL,2)) +
     &                   (ZPDEM(I)-Z_CELL) * 
     &                   (ZPDEM(I)-Z_CELL) )    
              !PRINT*,D_IJ,DISCRI
              IF( (D_IJ .LT. DISCRI) .OR. 
     &            (K .EQ. KPP(I) .AND. I_CELL .EQ. IPP(I)) )THEN
               !print*,D_IJ,DISCRI   
               N_COUNT2 = N_COUNT2 + 1    
               N_NEIGHBOUR = N_NEIGHBOUR + 1
               IC(N_NEIGHBOUR,1) = I_CELL
               IC(N_NEIGHBOUR,2) = K                                 
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
     &                 + (YPDEM(I)-PXY(I_N1,2)) * (YPDEM(I)-PXY(I_N1,2))
     &                 + (ZPDEM(I)-Z_N1) * (ZPDEM(I)-Z_N1))   
              DISN2=SQRT((XPDEM(I)-PXY(I_N2,1)) * (XPDEM(I)-PXY(I_N2,1))
     &                 + (YPDEM(I)-PXY(I_N2,2)) * (YPDEM(I)-PXY(I_N2,2))
     &                 + (ZPDEM(I)-Z_N2) * (ZPDEM(I)-Z_N2)) 
              DISN3=SQRT((XPDEM(I)-PXY(I_N3,1)) * (XPDEM(I)-PXY(I_N3,1))
     &                 + (YPDEM(I)-PXY(I_N3,2)) * (YPDEM(I)-PXY(I_N3,2))
     &                 + (ZPDEM(I)-Z_N3) * (ZPDEM(I)-Z_N3))
              DISN4=SQRT((XPDEM(I)-PXY(I_N4,1)) * (XPDEM(I)-PXY(I_N4,1))
     &                 + (YPDEM(I)-PXY(I_N4,2)) * (YPDEM(I)-PXY(I_N4,2))
     &                 + (ZPDEM(I)-Z_N4) * (ZPDEM(I)-Z_N4))
              DISN5=SQRT((XPDEM(I)-PXY(I_N1,1)) * (XPDEM(I)-PXY(I_N1,1))
     &                 + (YPDEM(I)-PXY(I_N1,2)) * (YPDEM(I)-PXY(I_N1,2))
     &                 + (ZPDEM(I)-Z_N5) * (ZPDEM(I)-Z_N5))   
              DISN6=SQRT((XPDEM(I)-PXY(I_N2,1)) * (XPDEM(I)-PXY(I_N2,1))
     &                 + (YPDEM(I)-PXY(I_N2,2)) * (YPDEM(I)-PXY(I_N2,2))
     &                 + (ZPDEM(I)-Z_N6) * (ZPDEM(I)-Z_N6)) 
              DISN7=SQRT((XPDEM(I)-PXY(I_N3,1)) * (XPDEM(I)-PXY(I_N3,1))
     &                 + (YPDEM(I)-PXY(I_N3,2)) * (YPDEM(I)-PXY(I_N3,2))
     &                 + (ZPDEM(I)-Z_N7) * (ZPDEM(I)-Z_N7))
              DISN8=SQRT((XPDEM(I)-PXY(I_N4,1)) * (XPDEM(I)-PXY(I_N4,1))
     &                 + (YPDEM(I)-PXY(I_N4,2)) * (YPDEM(I)-PXY(I_N4,2))
     &                 + (ZPDEM(I)-Z_N8) * (ZPDEM(I)-Z_N8))            

              IF( ( (DISN1 .LT. DISCRI) .OR. (DISN2 .LT. DISCRI)
     &         .OR. (DISN3 .LT. DISCRI) .OR. (DISN4 .LT. DISCRI)
     &         .OR. (DISN5 .LT. DISCRI) .OR. (DISN6 .LT. DISCRI)
     &         .OR. (DISN7 .LT. DISCRI) .OR. (DISN8 .LT. DISCRI) )    
     &         .OR. (K .EQ. KPP(I) .AND. I_CELL .EQ. IPP(I)) )THEN 
                N_COUNT2 = N_COUNT2 + 1              
                N_NEIGHBOUR = N_NEIGHBOUR + 1
                IC(N_NEIGHBOUR,1) = I_CELL
                IC(N_NEIGHBOUR,2) = K                                 
              ENDIF 
            
             ELSE
                print*,'please check the value of N_SCHEME!'
             ENDIF

             IF(N_NEIGHBOUR .GE. 30000)THEN
               PRINT*,'neighbour cell numbers is greater than 30000'
               PAUSE
             ENDIF
             
            ENDDO
            !print*,N_COUNT2
            IF(N_COUNT2 .EQ. 0)THEN
             EXIT  
            ELSEIF(J .EQ. NL)THEN
             PRINT*,'Error1 of NL in FCOUPLING.F'
             PAUSE
            ENDIF
            
           ENDDO           
            
          ENDIF
          
         ENDDO
         
         N_COUNT1 = N_NEIGHBOUR - N_NEIGHBOUR0 

        ENDDO
       
        !print*,KPP(I),((IC(J,K),K=1,2),J=1,N_NEIGHBOUR)
        !pause  
      
        SUM_WGAUSS = 0.0
        
        SUMU_WGAUSS = 0.0
        SUMV_WGAUSS = 0.0
        SUMW_WGAUSS = 0.0
        
        SUMPGX_WGAUSS = 0.0
        SUMPGY_WGAUSS = 0.0
        SUMPGZ_WGAUSS = 0.0
        
        SUMPORE_WGAUSS = 0.0
        SUMTKE_WGAUSS = 0.0
        SUMTDISS_WGAUSS = 0.0        
        
        SUM_WGAUSS2 = 0.0
        SUM_WGAUSS3 = 0.0
        DO J = 1, N_NEIGHBOUR
          V_MESH = AREA(IC(J,1))*DZ(IC(J,2))*DC(IC(J,1))
          Z_CELL = ZZ(IC(J,2)) * DC(IC(J,1)) + ELF(IC(J,1))           
          D_IJ=SQRT((XPDEM(I)-CXY(IC(J,1),1))*(XPDEM(I)-CXY(IC(J,1),1))
     &           + (YPDEM(I)-CXY(IC(J,1),2))*(YPDEM(I)-CXY(IC(J,1),2))
     &           + (ZPDEM(I)-Z_CELL)*(ZPDEM(I)-Z_CELL))

          IF( (D_IJ .LT. DISCRI1) .OR. (J .EQ. 1) )THEN          
          W_GAUSS(J)=1.0/(SQRT(2.0 * ACOS(-1.0)) * SIGMA_GAUSS) * 
     &        EXP(-D_IJ*D_IJ/(2.0*SIGMA_GAUSS*SIGMA_GAUSS))*V_MESH
          IF(N_NEIGHBOUR .EQ. 1)W_GAUSS(J)=1.0
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
          
          SUMPORE_WGAUSS = SUMPORE_WGAUSS +
     &                     PORE(IC(J,1),IC(J,2)) * W_GAUSS(J)  
          
          SUMTKE_WGAUSS = SUMTKE_WGAUSS +
     &                    TKE(IC(J,1),IC(J,2)) * W_GAUSS(J) 
          SUMTDISS_WGAUSS = SUMTDISS_WGAUSS +
     &                    TDISS(IC(J,1),IC(J,2)) * W_GAUSS(J)         
          ENDIF
          
          IF( (D_IJ .LT. DISCRI2) .OR. (J .EQ. 1) )THEN
          W_GAUSS2(J)=1.0/(SQRT(2.0 * ACOS(-1.0)) * SIGMA_GAUSS2) * 
     &        EXP(-D_IJ*D_IJ/(2.0*SIGMA_GAUSS2*SIGMA_GAUSS2))*V_MESH
          IF(N_NEIGHBOUR .EQ. 1)W_GAUSS2(J)=1.0
          SUM_WGAUSS2=SUM_WGAUSS2+W_GAUSS2(J)       
          ENDIF
          
          IF( (D_IJ .LT. DISCRI3) .OR. (J .EQ. 1) )THEN
          W_GAUSS3(J)=1.0/(SQRT(2.0 * ACOS(-1.0)) * SIGMA_GAUSS3) * 
     &        EXP(-D_IJ*D_IJ/(2.0*SIGMA_GAUSS3*SIGMA_GAUSS3))*V_MESH
          IF(N_NEIGHBOUR .EQ. 1)W_GAUSS3(J)=1.0
          SUM_WGAUSS3=SUM_WGAUSS3+W_GAUSS3(J)       
          ENDIF
          
        ENDDO        
        
        UFPDEM=SUMU_WGAUSS/SUM_WGAUSS-UPAVEDEM(I)
        VFPDEM=SUMV_WGAUSS/SUM_WGAUSS-VPAVEDEM(I)
        WFPDEM=SUMW_WGAUSS/SUM_WGAUSS-WPAVEDEM(I)
        URFPDEM = Sqrt(UFPDEM**2+VFPDEM**2+WFPDEM**2)
      
        PG_X=SUMPGX_WGAUSS/SUM_WGAUSS
        PG_Y=SUMPGY_WGAUSS/SUM_WGAUSS
        PG_Z=SUMPGZ_WGAUSS/SUM_WGAUSS
        
        POROSITY = DMIN1(DMAX1(SUMPORE_WGAUSS/SUM_WGAUSS,1.E-4),1.0)
        
        ALPHA_SOLID = DMAX1(1.0 - POROSITY,1.0E-10)
        
        TKE_P = SUMTKE_WGAUSS / SUM_WGAUSS
        TDISS_P = SUMTDISS_WGAUSS / SUM_WGAUSS        
                   
        IF(N_SYM_PG .EQ. 1 .AND. NSTEP .NE. 1) THEN
          F_PG_X(I) = - 1.00 * PG_X * V_PARTICLE(I)
          F_PG_Y(I) = - 1.00 * PG_Y * V_PARTICLE(I)
          F_PG_Z(I) = - 1.00 * PG_Z * V_PARTICLE(I)   
        ELSE
          F_PG_X(I) = 0.0
          F_PG_Y(I) = 0.0
          F_PG_Z(I) = DENSITY_F * V_PARTICLE(I) * GRAV
        ENDIF      
C--------------------------------------------- 
c---Reynolds number
        ReDEM = DPDEM(I) * URFPDEM / UMOL * POROSITY
   
c        IF(I_PARTICLE .EQ. 1)THEN      
c          ReDEM = (2. * RPDEM(I) ) * URFPDEM / UMOL
c        ELSE
c          ReDEM = (2. * RPDEM(I) ) * URFPDEM / UMOL * 
c     &            PORE(IPP(I),KPP(I))
c        ENDIF
c------drag coefficient according to Re          
        IF(ReDEM .GT. 0.0)THEN
            
          IF(CDMOD .EQ. 'MOD1') THEN
            CDDEM = (0.63 + 4.8 / SQRT(ReDEM))**2  ! (by DI FELICE, 1994)
          ELSEIF(CDMOD .EQ. 'MOD2') THEN
            CDDEM = 24.0/9.06**2*(9.06/SQRT(ReDEM)+1)**2  ! (by Abraham, 1994)
          ELSEIF(CDMOD .EQ. 'MOD3') THEN
            CDDEM = 24.0*(1.0+0.15*ReDEM**0.687)/ReDEM   ! (by Schiller, 1935)
          ENDIF
          
          CDDEM = DMAX1(CDDEM,0.44)   
          
        ELSE

          CDDEM = 0.0
          
        ENDIF     

        CONS_TCHI = 3.7 - 0.65*EXP(-(1.5-LOG10(ReDEM))**2/2.0)
        FUN_PORE = POROSITY**(1.0 - CONS_TCHI)
        
        IF(DRAGMOD .EQ. 'MOD1') THEN              
   	
          !IF(I_PARTICLE .EQ. 1)THEN
          !  FUN_PORE=1.0
          !ELSE
          !  CONS_TCHI=3.7-0.65*EXP(-(1.5-LOG10(ReDEM))**2/2.0)
          !  FUN_PORE=PORE(IPP(I),KPP(I))**(1-CONS_TCHI)
          !ENDIF      
          !FUN_PORE=1.0
        
          BELTA_DRAG = 0.75 * CDDEM * DENSITY_F * ALPHA_SOLID * 
     &                  URFPDEM / DPDEM(I) * FUN_PORE
          
        ELSEIF(DRAGMOD .EQ. 'MOD2')THEN
          
          IF(POROSITY .GT. 0.8)THEN

            BELTA_DRAG = 0.75 * CDDEM * DENSITY_F * ALPHA_SOLID *  
     &                   URFPDEM / DPDEM(I) * FUN_PORE 
            
          ELSEIF(POROSITY .LE. 0.8)THEN
              
            BELTA_DRAG = 150.0 * ALPHA_SOLID**2 / POROSITY * DENSITY_F
     &                   * UMOL / DPDEM(I)**2 + 1.75 * ALPHA_SOLID 
     &                   * DENSITY_F / DPDEM(I) * URFPDEM
                         
          ENDIF
          
        ENDIF
      
        TEMP = V_PARTICLE(I) / ALPHA_SOLID * BELTA_DRAG   
          
	  FDXDEM(I) = TEMP * UFPDEM

	  FDYDEM(I) = TEMP * VFPDEM 
	
	  FDZDEM(I) = TEMP * WFPDEM        
        
        IF(MPTMOD .EQ. 'MOD1')THEN
        
          TKEDEM(I) = TEMP * TKE_P
        
          TDISSDEM(I) = TEMP * TDISS_P
          
        ELSEIF(MPTMOD .EQ. 'MOD2')THEN
          
          TKEDEM(I) = TEMP * URFPDEM**2
        
          TDISSDEM(I) = TEMP * URFPDEM**2 * TDISS_P / 
     &                  DMAX1(TKE_P,1.0E-10)
          
        ENDIF

c        IF(I_PARTICLE .GT. 1 .AND. I_CV_PORE .EQ. 1 .AND. 
c     &      KPP(I) .NE. 0)THEN
c          CV_PORE = C_VIRTUAL * ( 1.0 + 2.0 * ( 1.0 - 
cC     &              PORE(IPP(I),KPP(I)) ) ) / PORE(IPP(I),KPP(I))
c     &              POROSITY ) ) / POROSITY     
c        ELSE
c          CV_PORE = C_VIRTUAL  
c        ENDIF
        
C---  附加质量力计算
        AUFPDEMX = ( UFPDEM - UFPDEM_PRE(I) ) / DTI
        AUFPDEMY = ( VFPDEM - VFPDEM_PRE(I) ) / DTI
        AUFPDEMZ = ( WFPDEM - WFPDEM_PRE(I) ) / DTI
        AUFPDEM = SQRT(AUFPDEMX**2+AUFPDEMY**2+AUFPDEMZ**2)    
        
        AC = URFPDEM**2 / DMAX1( (AUFPDEM * DPDEM(I)) , 1.E-10 )
        
        IF(C_VIRTUAL .GT. 0.001)THEN
          CV_PORE(I) = C_VIRTUAL
        ELSE
          CV_PORE(I) = 1.05 - 0.066 / ( AC**2 + 0.12 ) 
        ENDIF

C        IF(N_SYNFA(I) .EQ. 1)THEN
C         TEMP = CV_PORE * RMEAN(IPP(I),KPP(I)) * V_PARTICLE(I) / DTI 
C         FAXDEM(I) = TEMP * ( UFPDEM - UFPDEM_PRE(I) )
C         FAYDEM(I) = TEMP * ( VFPDEM - VFPDEM_PRE(I) )
C         FAZDEM(I) = TEMP * ( WFPDEM - WFPDEM_PRE(I) )
C        ELSE
C         N_SYNFA(I) = 1
C        ENDIF

        IF(N_SYNFA(I) .EQ. 1)THEN 
         TEMP = CV_PORE(I) * RMEAN(IPP(I),KPP(I)) * V_PARTICLE(I)   
         FAXDEM(I) = TEMP * AUFPDEMX
         FAYDEM(I) = TEMP * AUFPDEMY
         FAZDEM(I) = TEMP * AUFPDEMZ
         
         IF(MOD_FA .EQ. 'IMPLICIT')THEN   !DEM隐式处理附加质量力项时，计算附加质量力流速对时间导数部分DU/DT
          AUFDEMX = ( SUMU_WGAUSS/SUM_WGAUSS - UFDEM_PRE(I) ) / DTI
          AUFDEMY = ( SUMV_WGAUSS/SUM_WGAUSS - VFDEM_PRE(I) ) / DTI
          AUFDEMZ = ( SUMW_WGAUSS/SUM_WGAUSS - WFDEM_PRE(I) ) / DTI

          TEMP = CV_PORE(I) * RMEAN(IPP(I),KPP(I)) * V_PARTICLE(I)   
          FAFXDEM(I) = TEMP * AUFDEMX
          FAFYDEM(I) = TEMP * AUFDEMY
          FAFZDEM(I) = TEMP * AUFDEMZ
         ENDIF
         
        ELSE
         N_SYNFA(I) = 1
        ENDIF
      
        IF(MOD_FA .EQ. 'IMPLICIT')THEN   
         UFDEM_PRE(I) = SUMU_WGAUSS/SUM_WGAUSS
         VFDEM_PRE(I) = SUMV_WGAUSS/SUM_WGAUSS
         WFDEM_PRE(I) = SUMW_WGAUSS/SUM_WGAUSS  
        ENDIF

        UFPDEM_PRE(I) = UFPDEM
        VFPDEM_PRE(I) = VFPDEM
        WFPDEM_PRE(I) = WFPDEM    

C---  阻力+附加质量力   
        IF(1 .EQ. 0)THEN
         FTFX = FDXDEM(I) + FAXDEM(I)
         FTFY = FDYDEM(I) + FAYDEM(I)
         FTFZ = FDZDEM(I) + FAZDEM(I)
        ELSE
         FTFX = FDXDEM(I)
         FTFY = FDYDEM(I)
         FTFZ = FDZDEM(I)    
        ENDIF
         
C-------------------------------------------------  
        DO J = 1, N_NEIGHBOUR
            
         Z_CELL = ZZ(IC(J,2)) * DC(IC(J,1)) + ELF(IC(J,1))   
         D_IJ=SQRT((XPDEM(I)-CXY(IC(J,1),1))*(XPDEM(I)-CXY(IC(J,1),1))
     &           + (YPDEM(I)-CXY(IC(J,1),2))*(YPDEM(I)-CXY(IC(J,1),2))
     &           + (ZPDEM(I)-Z_CELL)*(ZPDEM(I)-Z_CELL))
c!$OMP CRITICAL
         IF( (D_IJ .LT. DISCRI2) .OR. (J .EQ. 1) )THEN
          FPFX_DEM(IC(J,1),IC(J,2))=FPFX_DEM(IC(J,1),IC(J,2))+
     &             W_GAUSS2(J)/SUM_WGAUSS2*FTFX
          FPFY_DEM(IC(J,1),IC(J,2))=FPFY_DEM(IC(J,1),IC(J,2))+
     &             W_GAUSS2(J)/SUM_WGAUSS2*FTFY
          FPFZ_DEM(IC(J,1),IC(J,2))=FPFZ_DEM(IC(J,1),IC(J,2))+
     &             W_GAUSS2(J)/SUM_WGAUSS2*FTFZ
          
          STKEDEM(IC(J,1),IC(J,2))=STKEDEM(IC(J,1),IC(J,2))+
     &             W_GAUSS2(J)/SUM_WGAUSS2*TKEDEM(I)
          STDISSDEM(IC(J,1),IC(J,2))=STDISSDEM(IC(J,1),IC(J,2))+
     &             W_GAUSS2(J)/SUM_WGAUSS2*TDISSDEM(I)          
         ENDIF

         IF( (D_IJ .LT. DISCRI3) .OR. (J .EQ. 1) )THEN
          VP_MESH(IC(J,1),IC(J,2))=VP_MESH(IC(J,1),IC(J,2))+
     &            W_GAUSS3(J)/SUM_WGAUSS3*V_PARTICLE(I)
         ENDIF        
c!$OMP END CRITICAL        
        ENDDO
         
       ELSE
          
        FDXDEM(I) = 0.0
        FDYDEM(I) = 0.0
        FDZDEM(I) = 0.0
        F_PG_X(I) = 0.0
        F_PG_Y(I) = 0.0
        F_PG_Z(I) = 0.0
        TKEDEM(I) = 0.0
        TDISSDEM(I) = 0.0      
        
        FAXDEM(I) = 0.
        FAYDEM(I) = 0.
        FAZDEM(I) = 0. 
        FAFXDEM(I) = 0.
        FAFYDEM(I) = 0.
        FAFZDEM(I) = 0.        
        N_SYNFA(I) = 0
        CV_PORE(I) = 0
        
       ENDIF

       ENDIF

      ENDDO           
!$OMP END DO
!$OMP END PARALLEL

      IF(1 .EQ. 0)THEN
      print*,maxval(abs(FDXDEM)),maxval(abs(FDYDEM)),maxval(abs(FDZDEM))
     &,maxloc(abs(FDXDEM))
      print*,maxval(abs(FAXDEM)),maxval(abs(FAYDEM)),maxval(abs(FAZDEM))
      print*,maxval(abs(FPFX_DEM)),maxval(abs(FPFY_DEM)),
     & maxval(abs(FPFZ_DEM))
      print*,abs(sum(FPFX_DEM)-sum(FDXDEM))/abs(sum(FDXDEM)),
     &       abs(sum(FPFy_DEM)-sum(FDYDEM))/abs(sum(FDYDEM)), 
     &       abs(sum(FPFz_DEM)-sum(FDZDEM))/abs(sum(FDZDEM))  
      print*,abs(sum(abs(VP_MESH))-sum(abs(V_PARTICLE)))/
     &       sum(abs(V_PARTICLE))
      ENDIF

c     Calculate the sum of forces that particles force on fluid (gridding)
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I,K,V_MESH)      
      Do K = 1 , KBM
!$OMP DO        
        Do I = 1 , IJM
          V_MESH = AREA(I)*DZ(K)*DC(I)
            
          IF(POREMODULE  .EQ. 'INCLUDE' .AND. I_PTF .EQ. 1) Then
            PORE(I,K)=DMIN1(DMAX1(1.0-VP_MESH(I,K)/V_MESH,0.0001),1.0)
          ENDIF
                              
          DEMFORCEX(I,K) = - DC(I)*FPFX_DEM(I,K)/RMEAN(I,K)/V_MESH
          DEMFORCEY(I,K) = - DC(I)*FPFY_DEM(I,K)/RMEAN(I,K)/V_MESH
          DEMFORCEZ(I,K) = - DC(I)*FPFZ_DEM(I,K)/RMEAN(I,K)/V_MESH
            
          STKEDEM(I,K) = CFKDEM * STKEDEM(I,K)/RMEAN(I,K)/V_MESH
          STDISSDEM(I,K) = (BF1(I,K)*CFWDEM1 + (1.-BF1(I,K))*CFWDEM2) * 
     &                      STDISSDEM(I,K) / RMEAN(I,K) / V_MESH
         Enddo
!$OMP END DO                    
      Enddo
!$OMP END PARALLEL        
c      WRITE(*,'(A13,F8.4)')'minvalPore = ',minval(PORE)
      Return
      End       
 
