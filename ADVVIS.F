C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (16/02/2011)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c      computes the advective,Coriolis,horizontal dispersive terms in the   #
c      momentum equation of u                                               #
c                                                                           #
c############################################################################
      Subroutine ADVVIS(CB1,CB2,CW1,CW2,CW3,SIGMA,VK,CV1,CDES,CT3,CT4,
     *                  CDES2)
      Include './Include/OCERM_INF'
      Include './Include/VORGEN_INF'

	Common/VISBLK/VISSOURCE(IJM,KB),VISMIN(IJM,KB)
	Parameter (TINT = 0.03, DESTYPE = 1, CW = 0.15)
C	Parameter (DWINTERFACE = 8.3E-3)
	Dimension HQ(IJE,KB),VQ(IJM,KB)
	Dimension TEMP(KBM)
	Dimension VAR_EDGE(IJE,KB),VAR_T(IJM,-1:KB+1),GRADX(IJM,KB),
     &	         GRADY(IJM,KB),GRADZ(IJM,KB)
	Dimension ELFX(IJM),ELFY(IJM),HX(IJM),HY(IJM)
c===========================================================================c
c                  initialiing arrays                                       c
c===========================================================================c
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(I,K) SHARED(HQ,VQ,VISF,VISSOURCE)
	Do K = 1, KB
!$OMP DO	
	   Do I = 1, IJM
           VISF(I,K) = 0.0
           VISSOURCE(I,K) = 0.0
	   Enddo
!$OMP END DO 
!$OMP DO	
	   Do I = 1, IJE
	      HQ(I,K) = 0.0
	   Enddo
!$OMP END DO 	   
!$OMP DO	
	   Do I = 1, IJM
	      VQ(I,K) = 0.0
	   Enddo
!$OMP END DO	   
	Enddo
!$OMP END PARALLEL	
      If(IWENOSCHEME .NE. 0) Then
         Do K = 1, KBM
            Do I = -1, NUM_GHOST, -1
               QGHOST(I,K) = VIS(INDEX_GHOST(I),K)
            Enddo
         Enddo   
	Endif
C============================================================================c
c            TVD schemes for the calculation of convective fluxes            c
c============================================================================c
C---- Variables on the midpoint of the cell surfaces
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,IL,IR)
!$OMP DO
      Do I = 1, IJM
	   Do K = 1, KBM
	      VAR_T(I,K) = VIS(I,K)
	   Enddo
	   VAR_T(I,0) = 2. * VIS(I,1) - VIS(I,2)
	   VAR_T(I,-1) = VAR_T(I,0)
	   VAR_T(I,KB) = 0.0
	   VAR_T(I,KB+1) = 0.0
	Enddo
!$OMP END DO
      Do K = 1, KBM
!$OMP DO
         Do I = 1, IJE
            VAR_EDGE(I,K) = 0.0
            IL = INDEX_EDGE(I,K,1)
            IR = INDEX_EDGE(I,K,2)
            If(CFM(I) .EQ. 1.0) Then
               VAR_EDGE(I,K) = (VIS(IL,K) + VIS(IR,K)) / 2.
            Endif
            If(CFM(I) .EQ. 0.0 .OR. CFM(I) .EQ. -1.0) Then
               If(ISLIP .EQ. 1) Then
                  If(IL .GT. 0) Then
                     If(CCM(IL) .EQ. 1.0) VAR_EDGE(I,K) = VIS(IL,K)
                  Endif   
                  If(IR .GT. 0) Then
                     If(CCM(IR) .EQ. 1.0) VAR_EDGE(I,K) = VIS(IR,K)
                  Endif
               Endif
            Endif
            If(CFM(I) .EQ. -3.0) Then
               If(IL .GT. 0) Then
                  If(CCM(IL) .EQ. 1.0) VAR_EDGE(I,K) = VIS(IL,K)
               Endif   
               If(IR .GT. 0) Then
                  If(CCM(IR) .EQ. 1.0) VAR_EDGE(I,K) = VIS(IR,K)
               Endif
            Endif
         Enddo
!$OMP END DO         
      Enddo
!$OMP DO 
      Do I = 1, IJM
         If(CCM(I) .EQ. 1.0) Then
            Do K = 2, KBM
               GRADZ(I,K) = (VIS(I,K-1) - VIS(I,K)) / DC(I) / DZZ(K-1)
            Enddo
            GRADZ(I,1) = GRADZ(I,2)
c            GRADZ(I,1) = 2. * GRADZ(I,2) - GRADZ(I,3)
            GRADZ(I,KB) = VIS(I,KBM) / (0.5 * DZ(KBM) * DC(I))
         Endif
      Enddo
!$OMP END DO      
!$OMP END PARALLEL	
	Call GRAD_XY(VAR_EDGE,GRADX,GRADY)
C	Call GRAD_Z(VAR_UD,GRADZ)
      Call TVDSCHEMEH(HQ,VIS,GRADX,GRADY,IH_TVD) 
      Call TVDSCHEMEV(VQ,VIS,VAR_T,GRADZ,IV_TVD)
c============================================================================c
C                Advection Descrization by 2nd Order TVD                     c
c============================================================================c
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP&         PRIVATE(I,J,K,XX,FV1,FV2,FV3,FT2,D1,D2,HMAX,HWN,DUP,
!$OMP&         SLOWERPART1,SLOWERPART2,SLOWER,DTIDLE,STIDLE,R,DIVXVIS,
!$OMP&         DIVYVIS,DIVZVIS,AAMF,FLUX1,FLUX2,FD,RDT,FT,FL,FE2,
!$OMP&         ALF1,FE1,FE,FB,FDT,FDITDLE,G,FW,ID,IS,ZDES,D2MAX,
!$OMP&         FRICTION_U,D1MAX,RD)
      If (ADVECT.EQ.'NON-LINEAR') Then
c----------------------------------------------------------------------------c
c                    horizontal advective terms                              c
c----------------------------------------------------------------------------c
	   Do K = 1, KBM
!$OMP DO	   
	      Do I = 1, IJM
	         If(CCM(I) .EQ. 1.0) Then
	            Do J = 1, CELL_POLYGEN(I)
		           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
					  VISF(I,K) = VISF(I,K) + DZ(K) *
     &				     HQ(CELL_SIDE(I,J,1),K) * CELL_CUV(I,J,6) *
     &				     (UN(CELL_SIDE(I,J,1),K) * CELL_CUV(I,J,7) +
     &                      VN(CELL_SIDE(I,J,1),K) * CELL_CUV(I,J,8))
	               Endif
	            Enddo
	         Endif
	      Enddo
!$OMP END DO NOWAIT	      
	   Enddo
!$OMP BARRIER
c----------------------------------------------------------------------------c
c                    vertical advective terms                                c
c----------------------------------------------------------------------------c
	   Do K = 1, KBM
!$OMP DO	   
	      Do I = 1, IJM
	         If(CCM(I) .EQ. 1.0) Then
	            VISF(I,K) = -VISF(I,K) - AREA(I) *
     &			          (VQ(I,K) * W(I,K) - VQ(I,K+1) * W(I,K+1))
	         Endif
	      Enddo
!$OMP END DO NOWAIT	      
	   Enddo
!$OMP BARRIER
	Endif
c============================================================================c
C                   Horizontal Diffusion Descrization by CS                  c
c============================================================================c
	Do K = 1, KBM
!$OMP DO	
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	        Do J = 1, CELL_POLYGEN(I)
	           FLUX1 = 0.0
	           FLUX2 = 0.0
	           If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) Then
	              AAMF = UMOL + 
     &		             (VIS(I,K) + VIS(CELL_SIDE(I,J,2),K)) / 2. 
		      FLUX1 = (DISCOE(I,J,1) - DISCOE(I,J,8)) * AAMF * 
     &			      (VIS(CELL_SIDE(I,J,2),K) - VIS(I,K))
		      FLUX2 = (DISCOE(I,J,7) - DISCOE(I,J,2)) * AAMF*
     &                        (VISV(CELL_SIDE(I,J,4),K) - 
     &                        VISV(CELL_SIDE(I,J,3),K))
 		      VISF(I,K) = VISF(I,K) + 
     &			          (FLUX1 + FLUX2) * DZ(K) / SIGMA
	            Endif
	            If(CFM(CELL_SIDE(I,J,1)) .EQ. 0.0 .OR.
     &		       CFM(CELL_SIDE(I,J,1)) .EQ. -1.0) Then
	               If(ISLIP .EQ. 0) Then
		          AAMF = UMOL + VIS(I,K)
	                  FLUX1 = -(DISCOE(I,J,1) - DISCOE(I,J,8)) * AAMF*
     &                            VIS(I,K)
			  VISF(I,K) = VISF(I,K) + FLUX1 * DZ(K) / SIGMA
	               Endif
	            Endif
		    If(CFM(CELL_SIDE(I,J,1)) .EQ. -2.0) Then  
                       AAMF = UMOL + VIS(I,K)
	               FLUX1 = -(DISCOE(I,J,1) - DISCOE(I,J,8)) * AAMF*
     &                       VIS(I,K)
		       VISF(I,K) = VISF(I,K) + FLUX1 * DZ(K) / SIGMA
	            Endif
	        Enddo
	      Endif
	   Enddo
!$OMP END DO NOWAIT	   
	Enddo
!$OMP BARRIER
C============================================================================C
c                       source and sink terms                                c
C============================================================================C
C       OPEN(1111,FILE='TEMP.DAT',STATUS='UNKNOWN')
!$OMP DO    
	Do I = 1, IJM
	   If(CCM(I) .EQ. 1.0) Then
	     D2MAX = 0.0
	     Do J = 1, CELL_POLYGEN(I)
	        If(CELL_CUV(I,J,6) .GT. D2MAX) D2MAX = CELL_CUV(I,J,6) 
           Enddo
	     Do K = 1, KBM
	        R = 0.0
	        XX = VIS(I,K) / UMOL
c-----        referring to the wall roughness
ccc	        XX = VIS(I,K) / UMOL + 0.5 * 30. * Z01(I) / 
ccc     &			 (DC(I) * (1. + ZZ(K)) + Z01(I))
		FV1 = XX ** 3. / (XX ** 3. + CV1 ** 3.)
ccc	        FV2 = 1. / (1. + XX / 5.) ** 3.
C	        FV3 = (1. + XX * FV1) *(3. * (1. + XX / 5.) + (XX / 5.) **2.)
C     &              / 5.0 / (1. + XX / 5.) ** 3.
	        FV2 = 1.0 - XX / (1. + XX * FV1)
	        FV3 = (1. + XX * FV1) * (1. - FV2) / (XX + 1.E-6)
	        FT2 = CT3 * Exp(-CT4 * XX ** 2.) 
ccc              FV2 = 1. - VIS(I,K) / (UMOL + VIS(I,K) * FV1) 
c	        ZSTAR = DC(I) * (1. + ZZ(K)) / 1.E-6 * 
c     &			    Sqrt(Sqrt(TBX(I)**2.+TBY(I)**2.)/RMEAN(I,KBM)) 
	        D1 = DC(I) * (1. + ZZ(K)) + Z01(I)
	        DUP = -1. * DC(I) * ZZ(K) 
       	      D1 = Dmin1(D1,D2D(I))
      	      D2 = CDES * Dmax1(D2MAX,DC(I)*DZ(K))
CCC	        D2 = CDES * Sqrt(4. * AREA(I) / ACOS(-1.0))!LES length-s cale DES97
	        HMAX = Dmax1(Sqrt(4. * AREA(I) / ACOS(-1.0)),DC(I)*DZ(K))
	        HWN = DC(I) * DZ(K)
	        SLOWERPART1 = 100.
			  SLOWERPART2 = (1. - 0.1355 / CW1 / VK ** 2. / 0.424 * 
     &                     FV2) / Dmax1(FV1,1.E-10)
	        SLOWER = Sqrt(Dmin1(SLOWERPART1,SLOWERPART2))

C              SLOWER = 1.0
c	        ZDES = 0.0
C	        D2 = CDES *   
C     &			 Dmin1(Dmax1(CW*D1,CW*HMAX,HWN),HMAX)   !LES length-scaleIDDES
	        If(DES .NE. 'NEGLECT') Then
                   If(DES .EQ. 'SADES  ') Goto 1
	           If(DES .EQ. 'SADDES ') Goto 2
	           If(DES .EQ. 'SAIDDES') Goto 3
                 If(DES .EQ. 'SAZDES ') Goto 4
C======================  MODE CHOICES  =========================================C
c-----           DES97 by Spalart et al.                                        C
1                Continue
                 DTIDLE = Dmin1(D1,D2*SLOWER)
			   STIDLE = 1.0 * TENSOR(I,K) +
     &                    VIS(I,K) / VK ** 2. / DTIDLE ** 2. * FV2
	           STIDLE = Dmax1(STIDLE, 0.0)
	           R = (VIS(I,K)+UMOL) / VK ** 2. / DTIDLE ** 2. / 
     *               Dmax1(STIDLE, 1.0E-10)
	           R = Dmin1(R,10.)
                   G = R + CW2 * R * (R ** 5. - 1.)
	           FW = G*((1. + CW3 ** 6.)/(G ** 6. + CW3 ** 6.)) ** (1./6.)
	           Goto 100
C===============================================================================C  
c-----           DDES by Spalart et al.                                         C
2                Continue
			   RD = (KM(I,K) + UMOL)/ VK ** 2. / D1 ** 2. / 
     &			    DMAX1(1.E-10,TENSOR(I,K))
C	           R = Dmin1(R,10.)
 	           FD = 1. - Dtanh((8. * RD) ** 3.)
C                 IF(I.EQ. 322) WRITE(1011,'(10F12.4)') DC(I)*ZZ(K),FD,RD
	           DTIDLE = D1 - FD * Dmax1(0.0, D1 - D2)
		       STIDLE = 1.0 * TENSOR(I,K) +             
     &                    VIS(I,K) / VK ** 2. / DTIDLE ** 2. * FV2
	           STIDLE = Dmax1(STIDLE, 0.0)
	           R = (VIS(I,K)+UMOL) / VK ** 2. / DTIDLE ** 2. / 
     *               Dmax1(STIDLE, 1.0E-10)
	           R = Dmin1(R,10.) 
                 G = R + CW2 * R * (R ** 5. - 1.)
	           FW = G*((1. + CW3 ** 6.)/(G ** 6. + CW3 ** 6.)) ** (1./6.)
	           Goto 100
C==============================================================================C  
c-----           IDDES by Spalart et al.                                       C                 
3                Continue
			     RDT = VIS(I,K) / VK ** 2. / D1 ** 2. / 
     &			   DMAX1(1.E-10,VORTICITY(I,K))
			     RDL = UMOL / VK ** 2. / D1 ** 2. / 
     &			   DMAX1(1.E-10,VORTICITY(I,K))
	           FT = Tanh((1.63**2.*RDT)**3.)
	           FL = Tanh((3.55**2.*RDL)**10.)
	           FE2 = 1.0 - Dmax1(FT,FL)
	           ALF1 = 0.25 - D1 / HMAX
	           If(ALF1 .GE. 0.0) Then
	              FE1 = 2. * Exp(-11.09*ALF1**2.)
	           Else
	              FE1 = 2. * Exp(-9.0*ALF1**2.)
	           Endif
	           FE = Dmax1((FE1-1.0), 0.0) * FE2
	           FB = Dmin1(2.*Exp(-9.*ALF1**2.),1.0)
	           FDT = 1. - Tanh((8.*RDT)**3.)
			     FDTIDLE = Dmax1(1.-FDT,FB)
	           DTIDLE = FDTIDLE * (1. + FE) * D1 + (1. - FDTIDLE) * D2
		         STIDLE = 1.0 * TENSOR(I,K) + 
     &                    VIS(I,K) / VK ** 2. / DTIDLE ** 2. * FV2
	           STIDLE = Dmax1(STIDLE, 1.E-10)	           
	           R = RDT
                 G = R + CW2 * R * (R ** 5. - 1.)
	           FW = G*((1. + CW3 ** 6.)/(G ** 6. + CW3 ** 6.)) ** (1./6.)
	           Goto 100
C=============================================================================C
c-----           ZDES
4                Continue
                 If(I_DES_MODE(I,K) .EQ. 0) Then
                    DTIDLE = D1
                    STIDLE = 1.0 * TENSOR(I,K) +
     &                    VIS(I,K) / VK ** 2. / DTIDLE ** 2. * FV2
	            STIDLE = Dmax1(STIDLE, 0.0)
	            R = (VIS(I,K)+UMOL) / VK ** 2. / DTIDLE ** 2. / 
     &                  Dmax1(STIDLE,1.0E-10)
           	    R = Dmin1(R,10.)
                    G = R + CW2 * R * (R ** 5. - 1.)
	            FW = G*((1.+CW3**6.)/(G**6.+CW3**6.))**(1./6.)
                 Endif  
C---------------  MODE I              
                 If(I_DES_MODE(I,K) .EQ. 1) Then
                    D2 = CDES*(AREA(I)*DC(I)*DZ(K))**(1./3.)
                    DTIDLE = Dmin1(D1,D2)
                    If(DTIDLE .EQ. D2) Then
                       FV1 = 1.0
                       FV2 = 0.0
                       FW = 1.0
C                   DTIDLE = CDES*(AREA(I)*DC(I)*DZ(K))**(1./3.)
                       STIDLE = 1.0 * TENSOR(I,K) +
     &                    VIS(I,K) / VK ** 2. / DTIDLE ** 2. * FV2
	                 STIDLE = Dmax1(STIDLE, 0.0)
c                       R = (VIS(I,K)) / VK ** 2. / DTIDLE ** 2. / 
c     *                     Dmax1(STIDLE, 1.0E-10)
c	                 R = Dmin1(R,10.)
c                       G = R + CW2 * R * (R ** 5. - 1.)
c	                 FW = G*((1.+CW3**6.)/(G**6.+CW3**6.)) ** (1./6.)
                    Else
                       STIDLE = 1.0 * TENSOR(I,K) +
     &                    VIS(I,K) / VK ** 2. / DTIDLE ** 2. * FV2
	                 STIDLE = Dmax1(STIDLE, 0.0)
                       R = (VIS(I,K)+UMOL) / VK ** 2. / DTIDLE ** 2. / 
     *                     Dmax1(STIDLE, 1.0E-10)
	                 R = Dmin1(R,10.)
                       G = R + CW2 * R * (R ** 5. - 1.)
	                 FW = G*((1.+CW3**6.)/(G**6.+CW3**6.)) ** (1./6.)
                    Endif
                 Endif
C---------------  MODE II                  
                 If(I_DES_MODE(I,K) .EQ. 2) Then
                    RD = (VIS(I,K)+UMOL)/ VK ** 2. / D1 ** 2. / 
     &			       Dmax1(1.E-10,TENSOR(I,K)) 
                    FD = 1. - Tanh((8.*RD)**3.)
                    If(FD .LE. 0.8) Then
                       D2 = D2  
                    Else
                       D2 = CDES*(AREA(I)*DC(I)*DZ(K))**(1./3.) 
                    Endif 
                    DTIDLE = D1 - FD * Dmax1(0.0,D1-D2)
                    STIDLE = 1.0 * TENSOR(I,K) +
     &                       VIS(I,K) / VK ** 2. / DTIDLE ** 2. * FV2
	              STIDLE = Dmax1(STIDLE, 0.0)
                    R = (VIS(I,K)+UMOL) / VK ** 2. / DTIDLE ** 2. / 
     *                   Dmax1(STIDLE, 1.0E-10)   
	              R = Dmin1(R,10.)
                    G = R + CW2 * R * (R ** 5. - 1.)
	              FW = G*((1.+CW3**6.)/(G**6.+CW3**6.)) ** (1./6.)
                 Endif
C---------------  MODE III                  
                 If(I_DES_MODE(I,K) .EQ. 3) Then
C                    D2 = CDES*(AREA(I)*DC(I)*DZ(K))**(1./3.)
c                    D2 = Dmax1(3**2.*D2/D1*D2,D2)
C                    DTIDLE=Dmin1(D1,CDES2*D2)
                    DTIDLE = D1
                    If(DTIDLE .GE. DWINTERFACE) Then
                       FV1 = 1.0
                       FV2 = 0.0
                       FW = 1.0
                       DTIDLE = CDES*(AREA(I)*DC(I)*DZ(K))**(1./3.)
                       STIDLE = 1.0 * TENSOR(I,K) +
     &                    VIS(I,K) / VK ** 2. / DTIDLE ** 2. * FV2
	                 STIDLE = Dmax1(STIDLE, 0.0)
	                
c                       R = (VIS(I,K)) / VK ** 2. / DTIDLE ** 2. / 
c     *                     Dmax1(STIDLE, 1.0E-10)
c	                 R = Dmin1(R,10.)
c                       G = R + CW2 * R * (R ** 5. - 1.)
c	                 FW = G*((1.+CW3**6.)/(G**6.+CW3**6.)) ** (1./6.)
                    Else
                       STIDLE = 1.0 * TENSOR(I,K) +
     &                    VIS(I,K) / VK ** 2. / DTIDLE ** 2. * FV2
	                 STIDLE = Dmax1(STIDLE, 0.0)
                       R = (VIS(I,K)+UMOL) / VK ** 2. / DTIDLE ** 2. / 
     *                     Dmax1(STIDLE, 1.0E-10)
	                 R = Dmin1(R,10.)
                       G = R + CW2 * R * (R ** 5. - 1.)
	                 FW = G*((1.+CW3**6.)/(G**6.+CW3**6.)) ** (1./6.)
	                
                    Endif
                 Endif    
                 Goto 100    
C============================================================================C                 
              Else
C-----           RANS MODEL
		        DTIDLE = D1
		        STIDLE = 1.0 * TENSOR(I,K) + 
     &                    VIS(I,K) / VK ** 2. / DTIDLE ** 2. * FV2
	          STIDLE = Dmax1(STIDLE, 0.0)
	          R = (VIS(I,K)+UMOL)/VK**2./DTIDLE**2./
     *                Dmax1(STIDLE,1.0E-10)
           	    R = Dmin1(R,10.)
                G = R + CW2 * R * (R ** 5. - 1.)
	          FW = G*((1. + CW3 ** 6.)/(G ** 6. + CW3 ** 6.)) ** (1./6.)
              Endif
C=====================  END OF MODE  =========================================C  
100	        Continue
C-----     source term Cb1
	        VISSOURCE(I,K) = DZ(K) * (1.0 * CB1 * STIDLE -
     &		                   CW1 * FW / DTIDLE **2. * VIS(I,K))
C	        VISF(I,K) = VISF(I,K) + CB1 * AREA(I) * DZ(K) * 
C     &			        STIDLE * VIS(I,K)
c-----     source term Cw1
c	        VISF(I,K) = VISF(I,K) - CW1 * AREA(I) * DZ(K) *
c     &			          FW * (VIS(I,K) / DTIDLE) ** 2.
CCC              VISMIN(I,K) = CB1 / CW1/ FW * DTIDLE ** 2. * STIDLE
c-----     source term Cb2
		VISF(I,K) = VISF(I,K) + CB2 * AREA(I) / SIGMA * DZ(K) *
C     &                (GRADX(I,K)**2.+GRADY(I,K)**2.+GRADZ(I,K)**2.)
     &                (GRADX(I,K) ** 2. +GRADY(I,K) ** 2. + 
     &                 (GRADZ(I,K) + GRADZ(I,K+1))**2. / 4.)
            Enddo
	   Endif
	Enddo
!$OMP END DO
!$OMP BARRIER	
C       STOP
C============================================================================c
c                    open boundary treatments                                c
c============================================================================c
!$OMP MASTER
C-----   elevation boundary condition
	   If(NUMEBC .NE. 0) Then
		 Do N = 1, NUMEBC
	        ID = IEBC(N)
	        IS = IEBCINX(N)
	        Do K = 1, KBM
	           VISF(ID,K) = 0.0
	           UNEBC = UR(ID,K) * CELL_CUV(ID,IS,7) + 
     & 			       VR(ID,K) * CELL_CUV(ID,IS,8)
	           Do J = 1, CELL_POLYGEN(ID)
		          If(CFM(CELL_SIDE(ID,J,1)) .EQ. 1.0) Then
					    VISF(ID,K) = VISF(ID,K) - DZ(K) *
     &				     HQ(CELL_SIDE(ID,J,1),K) * CELL_CUV(ID,J,6) *
     &				     (UN(CELL_SIDE(ID,J,1),K) * CELL_CUV(ID,J,7) +
     &                      VN(CELL_SIDE(ID,J,1),K) * CELL_CUV(ID,J,8))
	              Endif
	           Enddo
			   If(UNEBC .GT. 0.0) Then
	              VISF(ID,K) = VISF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &			             VIS(ID,K) * UNEBC
	           Else
                    VISF(ID,K) = VISF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &			             UNEBC * VISE(N)
	           Endif
	         Enddo
	      Enddo
	   Endif
C-----   astrotidal boundary condition
  	   If(NUMAST .NE. 0) Then
	     Do N = 1, NUMAST
	        ID = IABC(N)
	        IS = IABCINX(N)
	        Do K = 1, KBM
	           VISF(ID,K) = 0.0
	           UNAST = UR(ID,K) * CELL_CUV(ID,IS,7) + 
     & 			       VR(ID,K) * CELL_CUV(ID,IS,8)
	           Do J = 1, CELL_POLYGEN(ID)
		          If(CFM(CELL_SIDE(ID,J,1)) .EQ. 1.0) Then
					    VISF(ID,K) = VISF(ID,K) - DZ(K) *
     &				     HQ(CELL_SIDE(ID,J,1),K) * CELL_CUV(ID,J,6) *
     &				     (UN(CELL_SIDE(ID,J,1),K) * CELL_CUV(ID,J,7) +
     &                      VN(CELL_SIDE(ID,J,1),K) * CELL_CUV(ID,J,8))
	              Endif
	           Enddo
			   If(UNAST .GT. 0.0) Then
	              VISF(ID,K) = VISF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &			             VIS(ID,K) * UNAST
	           Else
                    VISF(ID,K) = VISF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &			             UNAST * Sqrt(1.5) * DC(ID) * DZ(K) *
     &                         UNAST * TINT
	           Endif
	         Enddo
	     Enddo
	   Endif
c-----   discharge boundary condition
	   If(NUMQBC .NE. 0) Then
	      Do N = 1, NUMQBC
	         ID = IQBC(N)
	         IS = IQBCINX(N)
	         ISS = CELL_SIDE(ID,IS,1)
	         Do K = 1, KBM
	            VISF(ID,K) = 0.0
 	            UNQBC = UN(ISS,K) * CELL_CUV(ID,IS,7) + 
     & 			        VN(ISS,K) * CELL_CUV(ID,IS,8)
                  Do J = 1, CELL_POLYGEN(ID)
		           If(CFM(CELL_SIDE(ID,J,1)) .EQ. 1.0) Then
				      VISF(ID,K) = VISF(ID,K) - DZ(K) *
     &				      HQ(CELL_SIDE(ID,J,1),K) * CELL_CUV(ID,J,6) *
     &				     (UN(CELL_SIDE(ID,J,1),K) * CELL_CUV(ID,J,7) +
     &                      VN(CELL_SIDE(ID,J,1),K) * CELL_CUV(ID,J,8))
	               Endif
	            Enddo
	            If(UNQBC .GT. 0.0) Then
	              VISF(ID,K) = VISF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &			                 VIS(ID,K) * UNQBC
	            Else
	              VISF(ID,K) = VISF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &                           UNQBC * VISQ(N,K)
                    
c     Sqrt(1.5) * DC(ID) * DZ(K) *
c     &                         Sqrt(UN(IS,K)**2.+VN(IS,K)**2.) * TINT
	            Endif
	         Enddo
	      Enddo
	   Endif
C-----   velocity boundary condition
	   If(NUMVBC .NE. 0) Then
	      Do N = 1, NUMVBC
	         ID = IVBC(N)
	         IS = IVBCINX(N)
	         Do K = 1, KBM
	            VISF(ID,K) = 0.0
 	            UNVBC = UN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,7) + 
     & 			        VN(CELL_SIDE(ID,IS,1),K) * CELL_CUV(ID,IS,8)
                  Do J = 1, CELL_POLYGEN(ID)
		           If(CFM(CELL_SIDE(ID,J,1)) .EQ. 1.0) Then
				      VISF(ID,K) = VISF(ID,K) - DZ(K) *
     &				      HQ(CELL_SIDE(ID,J,1),K) * CELL_CUV(ID,J,6) *
     &				     (UN(CELL_SIDE(ID,J,1),K) * CELL_CUV(ID,J,7) +
     &                      VN(CELL_SIDE(ID,J,1),K) * CELL_CUV(ID,J,8))
	               Endif
	            Enddo
	            If(UNVBC .GT. 0.0) Then
	              VISF(ID,K) = VISF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &			             VIS(ID,K) * UNVBC
	            Else
	              VISF(ID,K) = VISF(ID,K) - DZ(K) * CELL_CUV(ID,IS,6)*
     &                         UNVBC * Sqrt(1.5) * DC(ID) * DZ(K) *
     &                         Sqrt(UN(IS,K)**2.+VN(IS,K)**2.) * TINT
	            Endif
	         Enddo
	      Enddo
	   Endif
!$OMP END MASTER
!$OMP BARRIER
C============================================================================C
C                     Step forward in time                                   C
C============================================================================C
	Do K = 1, KBM
!$OMP DO
	   Do I = 1, IJM
	      If(CCM(I) .EQ. 1.0) Then
	         VISF(I,K) = VIS(I,K) * AREA(I) * DZ(K) + DTI * VISF(I,K)
	      Endif
         Enddo
!$OMP END DO     
	Enddo
!$OMP END PARALLEL
c======================  end subroutine program  ============================C
      Return
      End
