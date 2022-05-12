C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (18/06/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c      reading the mesh information and initial condition preparation for   #
c      a general model                                                      #
c                                                                           #
c############################################################################
      Subroutine SETDOM
      Include './Include/OCERM_INF'
      Dimension COM(80)
      Dimension EDGE_XY(IJE,2),
     &        QPPORE(IJM*KBM,4),XPORE(IJM*KBM),YPORE(IJM*KBM),
     &        ZPORE(IJM*KBM),PPORE(IJM*KBM)
      Dimension DEL1(KB),DEL2(KB),DEL3(KB),DEL4(KB)
      Dimension DELX1(KB),DELX2(KB),DELX3(KB),DELX4(KB),DELX5(KB)
      Dimension ZZZ(-1:KB+1)
      Write (IUPRT,5000)
C---------------------------------------------------------------------------C
C                    Establish depth array                                  C
C---------------------------------------------------------------------------C
      Read (IUGRD,5200) (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
C---------------------------------------------------------------------------C
C                    Initialize sigma levels                                C
C---------------------------------------------------------------------------C
      Read (IUGRD,5200) (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
C
      Read (IUGRD,5300) IKB
      Write (IUPRT,5400) IKB
      If (IKB .NE. KB) Then
         Write (*,5500) IKB, KB
         Stop
      Endif
C
      Do K = 1, IKB
        Read (IUGRD,5600) Z(K)
      Enddo
      Write (IUPRT,5600) (Z(K),K = 1,IKB)
C
      Do K = 1, KBM
         DZ(K) = Z(K) - Z(K+1)
         DZR(K) = 1. / DZ(K)
         ZZ(K) = .5 * (Z(K)+Z(K+1))
      Enddo
      Do K = 1, KBM - 1
        DZZ(K) = ZZ(K) - ZZ(K+1)
	Enddo
      DZZ(KBM) = ZZ(KBM) - Z(KB) 
c      OPEN(10101,FILE='TEMP1.DAT',STATUS='UNKNOWN')
C      WRITE(1010,'(100F14.6)') (DZ(K)*100,K=1,KBM)
C      STOP
C---------------------------------------------------------------------------C
C            Read depths on the grid points                                 C
C---------------------------------------------------------------------------C
      Read (IUGRD,5200)  (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
C
      Read (IUGRD,5700) IJP1
      Write (IUPRT,5800) IJP1
      If (IJP1 .NE. IJP) Then
         Write (*,5900) IJP1, IJP
         Stop
      Endif
	Do K = 1, IJP
	   Read (IUGRD,*) PXY(K,1), PXY(K,2), HP(K)
	Enddo
c--------------------------------------------------------------------------------C
C            Read the information of the meshes                                  C
C--------------------------------------------------------------------------------C
	Read (IUCUV,5200) (COM(I), I= 1, 80)
      Write (IUPRT,5100) (COM(I),I = 1,80)

	Read (IUCUV,*) IJM1
	If (IJM1 .NE. IJM) Then
         Write (*,5900) IJM1, IJM
         Stop
	Endif

	Do K = 1, IJM
	   Read (IUCUV,*) KCELL, CELL_POLYGEN(K), 
     &	    ((CELL_SIDE(K,I,J), J = 1, 4), I = 1, CELL_POLYGEN(K)),
     &        ((CELL_CUV(K,I,J), J= 1,8), I= 1, CELL_POLYGEN(K)),
     &        CXY(K,1), CXY(K,2), AREA(K)
c       IF(CXY(K,1) .GT. -6. .AND. CXY(K,1) .LT. -5.95 .and.
c     &    CXY(K,2) .GT. 0.0 .AND. CXY(K,2) .LT. 0.04) PRINT*,K,CXY(K,2)
	Enddo

C	DO I = 1,60
C        WRITE(10101,'(2I8,500F12.4)') 
C     &   I, 3,((1.+ZZ(K))*0.12,CXY(I,2),K=1, KBM)	
C	ENDDO
C      stop
	Read (IUCUV,5200) (COM(I), I= 1, 80)
      Write (IUPRT,5100) (COM(I),I = 1,80)

	Read (IUCUV,*) IJE1
	If (IJE1 .NE. IJE) Then
         Write (*,5900) IJE1, IJE
         Stop
	Endif
	Do K = 1, IJE
	   Read (IUCUV,*) KEDGE, (INDEX_EDGE(K,1,J), J = 1, 2),
     &                  (IEND_EDGE(K,J), J = 1, 2), CFM(K)
	Enddo

	Do K = 2, KB
	   Do I = 1, IJE
	      INDEX_EDGE(I,K,1) = INDEX_EDGE(I,1,1)
	      INDEX_EDGE(I,K,2) = INDEX_EDGE(I,1,2)
	   Enddo

	Enddo

	Read (IUCUV,5200) (COM(I), I= 1, 80)
      Write (IUPRT,5100) (COM(I),I = 1,80)

	Read (IUCUV,*) IJP1
	If (IJP1 .NE. IJP) Then
         Write (*,5900) IJP1, IJP
         Stop
	Endif
	Do K = 1, IJP
          Read (IUCUV,*) KPOINT, !#repeated variable#, WangJian, 2020-3-7 11:54:13 !    
     &	    POINT_CELL(K), (INDEX_POINT(K,J), J = 1, POINT_CELL(K))
	Enddo
	Do K = 1, IJM
	   If(CELL_POLYGEN(K) .EQ. 4) Then
	      CELL_INDEX(K,1) = CELL_SIDE(K,1,3)
	      CELL_INDEX(K,2) = CELL_SIDE(K,2,3)
	      CELL_INDEX(K,3) = CELL_SIDE(K,3,3)
	      CELL_INDEX(K,4) = CELL_SIDE(K,4,3)
	   Else
	      CELL_INDEX(K,1) = CELL_SIDE(K,1,3)
	      CELL_INDEX(K,2) = CELL_SIDE(K,2,3)
	      CELL_INDEX(K,3) = CELL_SIDE(K,3,3)
	      CELL_INDEX(K,4) = CELL_SIDE(K,3,3)
	   Endif
	Enddo
C=======================================================================C
C     instantaneous velocity output(compute the index of gauges)        c
C             created by Wang Jian (10/2018)                            c      
c=======================================================================c
      IF(IVEL_INS .NE.0 )THEN                                         
         OPEN(615,FILE='./Include/Gauge_XY.dat',status='unknown')       
         READ(615,*) NGAUGE1                                            
         IF (NGAUGE1.NE.NGAUGE) then                                    
             write(*,*)'the number of gauge is not correct'             
             STOP                                                       
         endif                                                          
         Do N=1,NGAUGE                                                  
             READ(615,*) GAUGE_XY(N,1),GAUGE_XY(N,2)                    
         ENDDO                                                          
         CLOSE(615)                                                     
c-----------check positions-------------------------------------------- 
           DO N = 1, NGAUGE                                             
             D_TOLMIN = 1.E10                                           
             Do K = 1, IJM                                              
              IF(SQRT((CXY(K,1)-GAUGE_XY(N,1))**2                       
     &           + (CXY(K,2)-GAUGE_XY(N,2))**2) .LT.D_TOLMIN )  THEN  
                 D_TOLMIN = SQRT((CXY(K,1)-GAUGE_XY(N,1))**2            
     &           + (CXY(K,2)-GAUGE_XY(N,2))**2)                         
                 INDEX_GAUGE(N)= K                                      
              ENDIF                                                     
             ENDDO                                                      
             D_TOL(N)= D_TOLMIN                                         
           Enddo                                                        
C----------Output gauge coordinates-------------------------------------        
         OPEN (611,FILE='./Results/VEL_INS_position.dat')               
         WRITE(611,*)'TITLE = "GAUGE_position"'                         
         WRITE(611,*)'VARIABLES =  "N", "X", "Y", "D_TOL", "D_ACTUAL" ' 
         WRITE(611,*)'ZONE T= "Model_position" I = ',NGAUGE             
         DO I=1,NGAUGE                                                  
         WRITE(611,'(I10,20F9.4)')INDEX_GAUGE(I),CXY(INDEX_GAUGE(I),1), 
     &         CXY(INDEX_GAUGE(I),2),D_TOL(I),                          
     &         SQRT((CXY(INDEX_GAUGE(I),1)-GAUGE_XY(I,1))**2            
     &         +(CXY(INDEX_GAUGE(I),2)-GAUGE_XY(I,2))**2)               
         ENDDO                                                          
         WRITE(611,*)'ZONE T= "Gauge_position" I = ',NGAUGE             
         DO I=1,NGAUGE                                                  
         WRITE(611,'(I10,20F9.4)')I,GAUGE_XY(I,1),GAUGE_XY(I,2),D_TOL(I)
     &         ,SQRT((CXY(INDEX_GAUGE(I),1)-GAUGE_XY(I,1))**2           
     &         +(CXY(INDEX_GAUGE(I),2)-GAUGE_XY(I,2))**2)               
         ENDDO                                                          
         CLOSE(611)                                                     
      ENDIF                                                             
C=======================================================================C 
c=======================================================================c
C=======================================================================C
C             concentration source                                      c
C             created by Wang Jian (5/2019)                             c      
c=======================================================================c 
      IF (WAQ .NE. 'NEGLECT')THEN                                       
        OPEN(616,FILE='./Include/Csource.QBC',status='old')             
        READ (616,*)                                                  
        READ (616,*) NUM_C1                                             
        IF (NUM_C1.NE. NUM_C) STOP                                      
        DO N=1,NUM_C                                                    
          READ(616,*) X_C(N),Y_C(N)                                     
        ENDDO                                                           
        DO K = 1, NUM_C                                                 
          D_TOLMIN = 1.E10                                              
          DO I=1,IJM                                                    
            if( SQRT( (CXY(I,1)-(X_C(K)))**2 + (CXY(I,2)-(Y_C(K)))**2  )
     &             .LT.  D_TOLMIN)THEN                                  
            D_TOLMIN =SQRT( (CXY(I,1)-(X_C(K)))**2 +                    
     &                      (CXY(I,2)-(Y_C(K)))**2  )                   
            INDEX_C(K) = I                                              
            ENDIF                                                       
          ENDDO                                                         
        ENDDO                                                           
                                                                        
        IF (IREADC==1)THEN                                              
          READ (616,*)                                                  
          READ (616,'(F14.6)')TOTAL_C                                   
          Do N = 1, NUM_C                                               
            Read (616,'(100F14.6)',End=617)                             
     &            (CQBCM(N,K), K = 1, KBM)                              
            DO K=1,KBM                                                  
              CQBCM(N,K)=CQBCM(N,K)*TOTAL_C/100.                        
            ENDDO                                                       
          Enddo                                                         
617       continue                                                      
          CLOSE(616)                                                    
        ENDIF                                                           
                                                                        
        WRITE(IUPRT,*)'Concentration Source:'                           
        DO N =1, NUM_C                                                  
          WRITE(IUPRT,'(A20,I10,2F14.6)')'IJM_INDEX=',INDEX_C(N),       
     &                 CXY(INDEX_C(N),1),CXY(INDEX_C(N),2)              
        ENDDO                                                           
        WRITE(IUPRT,'(A10,100I14)')'K',(INDEX_C(N),N=1, NUM_C)          
        DO K = 1, KBM                                                   
          WRITE(IUPRT,'(I10,100F14.6)')K,(CQBCM(N,K),N=1, NUM_C)        
        ENDDO                                                           
        WRITE(IUPRT,*)'End of concentration Source:'                    
        WRITE(IUPRT,*)                                                  
      ENDIF                                                           
      
C--------------------------------------------------------------------------------C
C   define the depths at the centre, middle point of edges and vertexs of cells   c
c---------------------------------------------------------------------------------c
	Do I = 1, IJE
	   N1 = IEND_EDGE(I,1)
	   N2 = IEND_EDGE(I,2) 
	   HS(I) = .5 * (HP(N1) + HP(N2))
	Enddo
	Do I = 1, IJM
	   If(CELL_POLYGEN(I) .EQ. 4) Then
c	     HC1 = (HP(CELL_SIDE(I,1,3)) + HP(CELL_SIDE(I,2,3)) + 
c     &		    HP(CELL_SIDE(I,3,3))) / 3.
c	     HC2 = (HP(CELL_SIDE(I,3,3)) + HP(CELL_SIDE(I,4,3)) + 
c     &		    HP(CELL_SIDE(I,1,3))) / 3.
c	     X1 = PXY(CELL_SIDE(I,1,3),1)
c	     Y1 = PXY(CELL_SIDE(I,1,3),2)
c	     X2 = PXY(CELL_SIDE(I,2,3),1)
c	     Y2 = PXY(CELL_SIDE(I,2,3),2)
c	     X3 = PXY(CELL_SIDE(I,3,3),1)
c	     Y3 = PXY(CELL_SIDE(I,3,3),2)
c	     X4 = PXY(CELL_SIDE(I,4,3),1)
c	     Y4 = PXY(CELL_SIDE(I,4,3),2)
c	     S1 = Abs(X1 * Y2 + X2 * Y3 + X3 * Y1 -
c     &              X1 * Y3 - X2 * Y1 - X3 * Y2) / 2.
c	     S2 = Abs(X1 * Y3 + X3 * Y4 + X4 * Y1 -
c     &              X1 * Y4 - X3 * Y1 - X4 * Y3) / 2.
	     HC(I) = (HP(CELL_SIDE(I,1,3)) + HP(CELL_SIDE(I,2,3)) +
     &		      HP(CELL_SIDE(I,3,3)) + HP(CELL_SIDE(I,4,3))) / 4.
	   Else
	     HC(I) = (HP(CELL_SIDE(I,1,3)) + HP(CELL_SIDE(I,2,3)) + 
     &		      HP(CELL_SIDE(I,3,3))) / 3.
	   Endif
	Enddo 
c----------------------------------------------------------------------C
	Do I = 1, IJM
	   DC(I) = HC(I) + EL(I)
	Enddo
	Do I = 1, IJE
	   If(INDEX_EDGE(I,1,1) .EQ. -999) Then
       	  DS(I) = HS(I) + EL(INDEX_EDGE(I,1,2))
	      ELFM(I) = EL(INDEX_EDGE(I,1,2))
	   Endif
	   If(INDEX_EDGE(I,1,2) .EQ. -999) Then
      	  DS(I) = HS(I) + EL(INDEX_EDGE(I,1,1))
	      ELFM(I) = EL(INDEX_EDGE(I,1,1))
	   Endif
	   If(INDEX_EDGE(I,1,1).NE.-999 .AND. INDEX_EDGE(I,1,2).NE.-999) 
     &   Then
            DS(I) = HS(I) + 0.5 * 
     &              (EL(INDEX_EDGE(I,1,1)) + EL(INDEX_EDGE(I,1,2)))
	      ELFM(I) = 0.5 * 
     &              (EL(INDEX_EDGE(I,1,1)) + EL(INDEX_EDGE(I,1,2)))
	   Endif
	Enddo
	Do I = 1, IJP
	   RTOL = 0.0
	   ELFV(I) = 0.0
	   Do J = 1, POINT_CELL(I)
	      R1 = Sqrt((PXY(I,1) - CXY(INDEX_POINT(I,J),1)) ** 2. +
     &	            (PXY(I,2) - CXY(INDEX_POINT(I,J),2)) ** 2.) 
c	      R1 = SQRT(AREA(INDEX_POINT(I,J)))
c	      RCOE(I,J) = R1
c		    RTOL = RTOL + R1
	      RCOE(I,J) = 1. / R1
  		    RTOL = RTOL + 1. / R1
            ELFV(I) = ELFV(I) + ELF(INDEX_POINT(I,J)) / R1
	   Enddo
	   ELFV(I) = ELFV(I) / RTOL
	   Do J = 1, POINT_CELL(I)
	      RCOE(I,J) = RCOE(I,J) / RTOL
	   Enddo
	Enddo
C---------------------------------------------------------------------------------C
C       coefficients for least squre method in calculating the gradient           c
c---------------------------------------------------------------------------------c
	Do I = 1, IJM
	   R11 = 0.0
	   R12 = 0.0
	   R22 = 0.0
	   Do J = 1, CELL_POLYGEN(I)
	      If(CELL_SIDE(I,J,2) .EQ. -999) Then
	        N1 = CELL_SIDE(I,J,3)
	        N2 = CELL_SIDE(I,J,4)
	        R110 = 4. * ((PXY(N1,1) + PXY(N2,1)) / 2.- CXY(I,1)) ** 2.
	        R120 = 4. * ((PXY(N1,1) + PXY(N2,1)) / 2. - CXY(I,1)) *
     &	 		        ((PXY(N1,2) + PXY(N2,2)) / 2. - CXY(I,2))
	      Else
	        R110 = (CXY(CELL_SIDE(I,J,2),1) - CXY(I,1)) ** 2.
	        R120 = (CXY(CELL_SIDE(I,J,2),1) - CXY(I,1)) * 
     &		       (CXY(CELL_SIDE(I,J,2),2) - CXY(I,2))
	      Endif
	      R11 = R11 + R110
            R12 = R12 + R120
	   Enddo
	   R11 = Sqrt(R11)
	   R12 = R12 / R11
	   Do J = 1, CELL_POLYGEN(I)
            N1 = CELL_SIDE(I,J,3)
	      N2 = CELL_SIDE(I,J,4)
	      If(CELL_SIDE(I,J,2) .EQ. -999) Then
	        R220 = ((PXY(N1,2) + PXY(N2,2)) - 2. * CXY(I,2) -
     &               ((PXY(N1,1) + PXY(N2,1)) - 2. * CXY(I,1)) *
     &               R12 / R11) **2.
	      Else
	        R220 = ((CXY(CELL_SIDE(I,J,2),2) - CXY(I,2)) - 
     &            (CXY(CELL_SIDE(I,J,2),1) - CXY(I,1)) * R12 / R11) **2.
	      Endif
	      R22 = R22 + R220
	   Enddo
	   R22 = Sqrt(R22)
	   Do J = 1, CELL_POLYGEN(I)
            N1 = CELL_SIDE(I,J,3)
	      N2 = CELL_SIDE(I,J,4)
	      If(CELL_SIDE(I,J,2) .EQ. -999) Then
	      WIX(I,J) = ((PXY(N1,1) + PXY(N2,1)) - 2. * CXY(I,1))/R11**2.
     &                 - R12 / R11 / R22 ** 2. * ((PXY(N1,2) + 
     &                 PXY(N2,2)) - 2. * CXY(I,2) - ((PXY(N1,1) + 
     &                 PXY(N2,1)) - 2. * CXY(I,1)) * R12 / R11)
	      WIY(I,J) = ((PXY(N1,2) + PXY(N2,2)) - 2. * CXY(I,2) - 
     &		         ((PXY(N1,1) + PXY(N2,1)) - 2. * CXY(I,1)) * 
     &                 R12 / R11) / R22 ** 2.
	      Else
	      WIX(I,J) = (CXY(CELL_SIDE(I,J,2),1) - CXY(I,1)) / R11 ** 2.-
     &                 R12 / R11 / R22 ** 2. * (CXY(CELL_SIDE(I,J,2),2)-
     &                 CXY(I,2) - (CXY(CELL_SIDE(I,J,2),1) - CXY(I,1)) *
     &                 R12 / R11)
	      WIY(I,J) = (CXY(CELL_SIDE(I,J,2),2) - CXY(I,2) - R12 / R11 *
     &		         (CXY(CELL_SIDE(I,J,2),1) - CXY(I,1))) / R22 ** 2.
	      Endif
	   Enddo
	Enddo
C---------------------------------------------------------------------------------C
C       Discretization coefficients for the difference on cell faces              c
c---------------------------------------------------------------------------------c
	Do I = 1, IJM
	   Do J = 1, CELL_POLYGEN(I)

	      If(CELL_SIDE(I,J,2) .EQ. -999) Then

               XMID = 0.5 * (PXY(CELL_SIDE(I,J,3),1) + 
     &				       PXY(CELL_SIDE(I,J,4),1))
	         YMID = 0.5 * (PXY(CELL_SIDE(I,J,3),2) + 
     &				       PXY(CELL_SIDE(I,J,4),2))

	         DISCOE(I,J,1) = CELL_CUV(I,J,7) * CELL_CUV(I,J,6) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,4) /
     &                         Sqrt((CXY(I,1) - XMID) ** 2. +
     &                              (CXY(I,2) - YMID) ** 2.)

	         DISCOE(I,J,2) = CELL_CUV(I,J,7) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,3) 

	         DISCOE(I,J,3) = CELL_CUV(I,J,7) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,1) 

	         DISCOE(I,J,4) = CELL_CUV(I,J,7) * CELL_CUV(I,J,6) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,2) /
     &                         Sqrt((CXY(I,1) - XMID) ** 2. +
     &                              (CXY(I,2) - YMID) ** 2.)

	         DISCOE(I,J,5) = CELL_CUV(I,J,8) * CELL_CUV(I,J,6) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,4) /
     &                         Sqrt((CXY(I,1) - XMID) ** 2. +
     &                              (CXY(I,2) - YMID) ** 2.)

	         DISCOE(I,J,6) = CELL_CUV(I,J,8) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,3)

	         DISCOE(I,J,7) = CELL_CUV(I,J,8) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,1)

	         DISCOE(I,J,8) = CELL_CUV(I,J,8) * CELL_CUV(I,J,6) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,2) /
     &                         Sqrt((CXY(I,1) - XMID) ** 2. +
     &                              (CXY(I,2) - YMID) ** 2.)

	      Else

	         DISCOE(I,J,1) = CELL_CUV(I,J,7) * CELL_CUV(I,J,6) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,4) /
     &            Sqrt((CXY(I,1) - CXY(CELL_SIDE(I,J,2),1)) ** 2. +
     &                 (CXY(I,2) - CXY(CELL_SIDE(I,J,2),2)) ** 2.)

	         DISCOE(I,J,2) = CELL_CUV(I,J,7) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,3) 

	         DISCOE(I,J,3) = CELL_CUV(I,J,7) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,1) 

	         DISCOE(I,J,4) = CELL_CUV(I,J,7) * CELL_CUV(I,J,6) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,2) /
     &            Sqrt((CXY(I,1) - CXY(CELL_SIDE(I,J,2),1)) ** 2. +
     &                 (CXY(I,2) - CXY(CELL_SIDE(I,J,2),2)) ** 2.)

	         DISCOE(I,J,5) = CELL_CUV(I,J,8) * CELL_CUV(I,J,6) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,4) /
     &            Sqrt((CXY(I,1) - CXY(CELL_SIDE(I,J,2),1)) ** 2. +
     &                 (CXY(I,2) - CXY(CELL_SIDE(I,J,2),2)) ** 2.)

	         DISCOE(I,J,6) = CELL_CUV(I,J,8) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,3)

	         DISCOE(I,J,7) = CELL_CUV(I,J,8) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,1)

	         DISCOE(I,J,8) = CELL_CUV(I,J,8) * CELL_CUV(I,J,6) / 
     &			             CELL_CUV(I,J,5) * CELL_CUV(I,J,2) /
     &            Sqrt((CXY(I,1) - CXY(CELL_SIDE(I,J,2),1)) ** 2. +
     &                 (CXY(I,2) - CXY(CELL_SIDE(I,J,2),2)) ** 2.)

	      Endif
	   Enddo
	Enddo
	Do K = 1, KBM - 1
	   QZ_INTER(K,1) = (ZZ(K) - Z(K+1)) / (Z(K) - Z(K+1)) * 
     &                   (ZZ(K) - Z(K+2)) / (Z(K) - Z(K+2))
         QZ_INTER(K,2) = (ZZ(K) - Z(K)) / (Z(K+1) - Z(K)) * 
     &                   (ZZ(K) - Z(K+2)) / (Z(K+1) - Z(K+2))	
         QZ_INTER(K,3) = (ZZ(K) - Z(K)) / (Z(K+2) - Z(K)) * 
     &                   (ZZ(K) - Z(K+1)) / (Z(K+2) - Z(K+1)) 
	Enddo
	QZ_INTER(KBM,1) = (ZZ(KBM) - Z(KBM)) / (Z(KBM-1) - Z(KBM)) * 
     &                  (ZZ(KBM) - Z(KB)) / (Z(KBM-1) - Z(KB))
      QZ_INTER(KBM,2) = (ZZ(KBM) - Z(KBM-1)) / (Z(KBM) - Z(KBM-1)) * 
     &                  (ZZ(KBM) - Z(KB)) / (Z(KBM) - Z(KB))	
      QZ_INTER(KBM,3) = (ZZ(KBM) - Z(KBM-1)) / (Z(KB) - Z(KBM-1)) * 
     &                  (ZZ(KBM) - Z(KBM)) / (Z(KB) - Z(KBM)) 
     
	Do K = 2, KBM
	   QZ_INTER_D(K,1) = (ZZ(K) - Z(K)) / (Z(K-1) - Z(K)) * 
     &                   (ZZ(K) - Z(K+1)) / (Z(K-1) - Z(K+1))
         QZ_INTER_D(K,2) = (ZZ(K) - Z(K-1)) / (Z(K) - Z(K-1)) * 
     &                   (ZZ(K) - Z(K+1)) / (Z(K) - Z(K+1))	
         QZ_INTER_D(K,3) = (ZZ(K) - Z(K-1)) / (Z(K+1) - Z(K-1)) * 
     &                   (ZZ(K) - Z(K)) / (Z(K+1) - Z(K))
	Enddo
	QZ_INTER_D(1,1) = QZ_INTER(1,1) 
      QZ_INTER_D(1,2) = QZ_INTER(1,2) 
      QZ_INTER_D(1,3) = QZ_INTER(1,3)
      
     	Do K = 2, KBM - 1
	   W_INTER(K,1) = (Z(K) - ZZ(K)) / (ZZ(K-1) - ZZ(K)) * 
     &                  (Z(K) - ZZ(K+1)) / (ZZ(K-1) - ZZ(K+1))
         W_INTER(K,2) = (Z(K) - ZZ(K-1)) / (ZZ(K) - ZZ(K-1)) * 
     &                  (Z(K) - ZZ(K+1)) / (ZZ(K) - ZZ(K+1))	
         W_INTER(K,3) = (Z(K) - ZZ(K-1)) / (ZZ(K+1) - ZZ(K-1)) * 
     &                  (Z(K) - ZZ(K)) / (ZZ(K+1) - ZZ(K))
	Enddo
	
	W_INTER(KBM,1) = (Z(KBM) - ZZ(KBM)) / (ZZ(KBM-1) - ZZ(KBM)) * 
     &                 (Z(KBM) - Z(KB)) / (ZZ(KBM-1) - Z(KB))
      W_INTER(KBM,2) = (Z(KBM) - ZZ(KBM-1)) / (ZZ(KBM) - ZZ(KBM-1)) * 
     &                 (Z(KBM) - Z(KB)) / (ZZ(KBM) - Z(KB))	
      W_INTER(KBM,3) = (Z(KBM) - ZZ(KBM-1)) / (Z(KB) - ZZ(KBM-1)) * 
     &                 (Z(KBM) - ZZ(KBM)) / (Z(KB) - ZZ(KBM))   
c---------------------------------------------------------------------------------c
c              Corilis force                                                      c
c---------------------------------------------------------------------------------c	
      COR = 0.0
	If(CORLIS .EQ. 'INCLUDE') Then
	   Do I = 1, IJM
	      Call XY2JW(CXY(I,1),CXY(I,2),ALFD,BETD)
	      COR(I) = BETD
            COR(I)=2.*7.292E-5*Sin(COR(I)*2.*3.14159/360.)
	   Enddo
	Endif
C---------------------------------------------------------------------------------C
C              Define the wet and dry cells                                       c
c---------------------------------------------------------------------------------c
	Do I = 1, IJM
	   If(DC(I) .LE. 0.0) Then
	      CCM(I) = 0.0
	   Else
	      CCM(I) = 1.0
	   Endif
	Enddo
	Do I = 1, IJE
	   If(DS(I) .GT. 0.0 .AND. CFM(I) .EQ. 0.0) CFM(I) = 1.0
	Enddo
C	Do I = 1, IJE
C	   If(DS(I) .LE. 0.0) Then
C	      CFM(I) = 0.0
C	   Else
C	      CFM(I) = 1.0
C	   Endif
C	   If(INDEX_EDGE(I,1,1) .EQ. -999 .OR.
C     &      INDEX_EDGE(I,1,2) .EQ. -999) CFM(I) = -1.0
C	Enddo
	Do I = 1, IJM
	   If(CCM(I) .EQ. 0.0) Then
	      Do J = 1, CELL_POLYGEN(I)
	         If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) 
     &			CFM(CELL_SIDE(I,J,1)) = 0.0
	      Enddo
	   Endif
	Enddo
C---------------------------------------------------------------------------------C
C        Distances to the solid walls                                             c
c---------------------------------------------------------------------------------c
	NUM_SOLID_WALL = 0
	Do K = 1, IJE
	   If(CFM(K) .EQ. -2.0) Then
	      NUM_SOLID_WALL = NUM_SOLID_WALL + 1
	      EDGE_XY(NUM_SOLID_WALL,1) = 
     &	            (PXY(IEND_EDGE(K,1),1) + PXY(IEND_EDGE(K,2),1)) / 2.
	      EDGE_XY(NUM_SOLID_WALL,2) = 
     &	            (PXY(IEND_EDGE(K,1),2) + PXY(IEND_EDGE(K,2),2)) / 2.
	   Endif
	   If(CFM(K) .EQ. -1.0 .OR. CFM(K) .EQ. 0.0) Then
	      If(ISLIP .EQ. 0) Then
	      NUM_SOLID_WALL = NUM_SOLID_WALL + 1
	      EDGE_XY(NUM_SOLID_WALL,1) = 
     &	            (PXY(IEND_EDGE(K,1),1) + PXY(IEND_EDGE(K,2),1)) / 2.
	      EDGE_XY(NUM_SOLID_WALL,2) = 
     &	            (PXY(IEND_EDGE(K,1),2) + PXY(IEND_EDGE(K,2),2)) / 2.
            Endif
	   Endif
	Enddo
	Do I = 1, IJM
	   DISMIN = 1.E10
	   Do J = 1, NUM_SOLID_WALL
	      DIS = Sqrt((CXY(I,1) - EDGE_XY(J,1)) ** 2. + 
     &	               (CXY(I,2) - EDGE_XY(J,2)) **2.)
            If(DIS .LT. DISMIN) DISMIN = DIS
	   Enddo
	   D2D(I) = DISMIN
      Enddo  
c----------------------------------------------------------------------------c
c        Initialing the DEM simulation                                       c
c----------------------------------------------------------------------------c
      PI_CON = 2. * ACos(0.)

      If(DEM .NE. 'NEGLECT') THEN
        
        CALL DEMSETTING    
        CALL SEARCH
        Call IPPKPPINITIAL
        
      ENDIF
      
c----------------------------------------------------------------------------c
c        Initialing the PORE module                                          c
c----------------------------------------------------------------------------c
       Do I = 1, IJM
         Do K = 1, KBM                           
             PORE(I,K) = 1.0           
         Enddo
      Enddo      
      
      If(POREMODULE  .EQ. 'INCLUDE') Then
        If(DEM .EQ. 'NEGLECT') Then
         Open(1002,file='./Include/PORE_MODULE.DAT')
         Read(1002,*)
         Read(1002,*)
         

         Do I = 1, IJM * KBM
             Read(1002,*) QPPORE(I,1:4)
         Enddo
          DO I = 1, IJM * KBM
             XPORE(I) = QPPORE(I,1)
             YPORE(I) = QPPORE(I,2)
             ZPORE(I) = QPPORE(I,3)
             PPORE(I) = QPPORE(I,4)
          Enddo
          J = 0
          Do K = 1, KBM
             Do I = 1, IJM
                J = J + 1
                XPPORE(I,K) = XPORE(J)
                YPPORE(I,K) = YPORE(J)
                ZPPORE(I,K) = ZPORE(J)
                PORE(I,K) = PPORE(J)
                APU(I,K) = 4 * (1-PORE(I,K)) /(PI * Dia)
             Enddo
          Enddo
          Close(1002)


c-----Save the information  of the poremodule
      !Call IPKPPOREMODULE    !Caliculate the gridding number that pore is belong	
        endif
	endif      
C---------------------------------------------------------------------------------c
c=================================================================================c
c               ZONES INDEX FOR ZDES                                              C
C---------------------------------------------------------------------------------C
      If(DES .EQ. 'SAZDES ') Then
         Open(1001,file='./Include/ZONE_ZDES.DAT',STATUS='UNKNOWN') 
         Do I = 1, IJM
            Read(1001,*) (I_DES_MODE(I,K), K = 1, KBM) 
         Enddo 
         Close(1001)
         Call INIVOR
      Endif 
c=================================================================================c
c                       WENO COEFFICIENTS                                         c
c=================================================================================c
C----  POSITIVE INTERPLATION
      Do K = 1, KBM
         ZZZ(K) = ZZ(K)
      Enddo
      ZZZ(0) = -ZZ(1)
      ZZZ(KB) = -1.0  + (-1.0 - ZZ(KBM))
      ZZZ(KB+1) = -1.0 + 3. * (-1.0 - ZZ(KBM)) 
      
      Do K = 2, KB -1 

         ALF_PLUX(K,1,1) = (Z(K) - ZZZ(K+1)) * (Z(K) - ZZZ(K)) / 
     &                     ((ZZZ(K+2) - ZZZ(K+1)) * (ZZZ(K+2) - ZZZ(K)))
         ALF_PLUX(K,1,2) = (Z(K) - ZZZ(K+2)) * (Z(K) - ZZZ(K)) / 
     &                     ((ZZZ(K+1) - ZZZ(K+2)) * (ZZZ(K+1) - ZZZ(K)))
         ALF_PLUX(K,1,3) = (Z(K) - ZZZ(K+2)) * (Z(K) - ZZZ(K+1)) / 
     &                     ((ZZZ(K) - ZZZ(K+2)) * (ZZZ(K) - ZZZ(K+1)))

         ALF_PLUX(K,2,1) = (Z(K) - ZZZ(K)) * (Z(K) - ZZZ(K-1)) / 
     &                     ((ZZZ(K+1) - ZZZ(K)) * (ZZZ(K+1) - ZZZ(K-1)))
         ALF_PLUX(K,2,2) = (Z(K) - ZZZ(K+1)) * (Z(K) - ZZZ(K-1)) / 
     &                     ((ZZZ(K) - ZZZ(K+1)) * (ZZZ(K) - ZZZ(K-1)))
         ALF_PLUX(K,2,3) = (Z(K) - ZZZ(K+1)) * (Z(K) - ZZZ(K)) / 
     &                     ((ZZZ(K-1) - ZZZ(K+1)) * (ZZZ(K-1) - ZZZ(K)))

         ALF_PLUX(K,3,1) = (Z(K) - ZZZ(K-1)) * (Z(K) - ZZZ(K-2)) / 
     &                     ((ZZZ(K) - ZZZ(K-1)) * (ZZZ(K) - ZZZ(K-2)))
         ALF_PLUX(K,3,2) = (Z(K) - ZZZ(K)) * (Z(K) - ZZZ(K-2)) / 
     &                     ((ZZZ(K-1) - ZZZ(K)) * (ZZZ(K-1) - ZZZ(K-2)))
         ALF_PLUX(K,3,3) = (Z(K) - ZZZ(K)) * (Z(K) - ZZZ(K-1)) / 
     &                     ((ZZZ(K-2) - ZZZ(K)) * (ZZZ(K-2) - ZZZ(K-1)))
       
      Enddo

      C_PLUX(1) = 1. / 10.
      C_PLUX(2) = 6.0 / 10.
      C_PLUX(3) = 3.0 / 10.

C---   NEGATIVE INTERPLATION
      Do K = 1, KBM
         ZZZ(K) = ZZ(K)
      Enddo
      ZZZ(0) = -ZZ(1)
      ZZZ(-1) = -3. * ZZ(1)
      ZZZ(KB) = -1.0  + (-1.0 - ZZ(KBM))
      
      Do K = 2, KB - 1
         ALF_MINU(K,1,1) = (Z(K) - ZZZ(K-2)) * (Z(K) - ZZZ(K-1)) / 
     &                   ((ZZZ(K-3) - ZZZ(K-2)) * (ZZZ(K-3) - ZZZ(K-1)))
         ALF_MINU(K,1,2) = (Z(K) - ZZZ(K-3)) * (Z(K) - ZZZ(K-1)) / 
     &                   ((ZZZ(K-2) - ZZZ(K-3)) * (ZZZ(K-2) - ZZZ(K-1)))
         ALF_MINU(K,1,3) = (Z(K) - ZZZ(K-3)) * (Z(K) - ZZZ(K-2)) / 
     &                   ((ZZZ(K-1) - ZZZ(K-3)) * (ZZZ(K-1) - ZZZ(K-2)))
         
         ALF_MINU(K,2,1) = (Z(K) - ZZZ(K-1)) * (Z(K) - ZZZ(K)) / 
     &                     ((ZZZ(K-2) - ZZZ(K-1)) * (ZZZ(K-2) - ZZZ(K)))
         ALF_MINU(K,2,2) = (Z(K) - ZZZ(K-2)) * (Z(K) - ZZZ(K)) / 
     &                     ((ZZZ(K-1) - ZZZ(K-2)) * (ZZZ(K-1) - ZZZ(K)))
         ALF_MINU(K,2,3) = (Z(K) - ZZZ(K-2)) * (Z(K) - ZZZ(K-1)) / 
     &                     ((ZZZ(K) - ZZZ(K-2)) * (ZZZ(K) - ZZZ(K-1)))
          
         ALF_MINU(K,3,1) = (Z(K) - ZZZ(K)) * (Z(K) - ZZZ(K+1)) / 
     &                     ((ZZZ(K-1) - ZZZ(K)) * (ZZZ(K-1) - ZZZ(K+1)))
         ALF_MINU(K,3,2) = (Z(K) - ZZZ(K-1)) * (Z(K) - ZZZ(K+1)) / 
     &                     ((ZZZ(K) - ZZZ(K-1)) * (ZZZ(K) - ZZZ(K+1)))
         ALF_MINU(K,3,3) = (Z(K) - ZZZ(K-1)) * (Z(K) - ZZZ(K)) / 
     &                     ((ZZZ(K+1) - ZZZ(K-1)) * (ZZZ(K+1) - ZZZ(K)))

      Enddo
      
      C_MINU(1) = 1. / 10.
      C_MINU(2) = 6.0 / 10.
      C_MINU(3) = 3.0 / 10.
      
      If(IWENOSCHEME .NE. 0) Then
C----------------------------------------------------------------------------C
C       Znoal WENO SCHEME                                                    C
C----------------------------------------------------------------------------C
        Open(1001,file='./Include/ZONE_WENO.DAT',STATUS='UNKNOWN')
        Read(1001,'(100I10)') (IWENO(I),I = 1, IJM)
        Close(1001)
c         ENDDO

         Do I = 1, IJM
            L_DEGREE(I) = IWENO(I) * (IWENO(I) + 1) / 2 - 1
         Enddo                                        
         Call WENO
      Endif   
C----------------------------------------------------------------------------C
C         Immersed boundary method                                           C
C----------------------------------------------------------------------------C

      If (MODULEIBM .NE. 'NEGLECT')   Call IBMINIT
        
      Return
5000  Format (/'...... model starting up from inital conditions .......'
     *    )
5100  Format (/1X,80A1/)
5200  Format (80A1)
5300  Format (I5)
5400  Format (' kb = ',I5,/)
5500  Format (//' number of sigma levels in model_grid',I5,' (ikb)'/
     *    '           does not equal'/
     *    ' number of sigma levels in GEFM_INC   ',I5,' (kb)'/
     *    ' please correct this problem and try again'//)
5600  Format (8F10.6)
5700  Format (2I10)
5800  Format (' ijm = ',I5)
5900  Format (//'     model_grid i-index',I5,' (ijm)',/
     *    '        does not equal'/'     GEFM_INC    i-index',I5,
     *    ' (im)'/' please correct this problem and try again'//)
6000  Format (//'     model_grid j-index',I5,' (ijy)',/
     *    '        does not equal'/'     GEFM_INC    j-index',I5,
     *    ' (jm)'/' please correct this problem and try again'//)
       End

	Subroutine XY2JW(Y,X,ALFD,BETD)
      Double precision alfd,betd,x,y
	Double precision a,b,alf0,pi,bet,alf
	Double precision e,eq,es,eqs,betf,rf
	Double precision sb,sb2,tb,cb,tb2,alf00
	Double precision xn,xm,xa,xa2,xa4,d,d2
	Double precision f1,f2,e4,e6,f3,f4
	pi=2.*acos(0.)
	alf=alfd*pi/180.
      bet=betd*pi/180.
c--------local sys E122
	alf00=123.
	alf0=alf00*pi/180.
c bj54
	a=6378245.
	b=6356863.0188
c wps84
c     a=6378137.
c	b=6356752.3142
c
	eq=1.-(b*b/(a*a))
	e4=eq*eq
	e6=e4*eq
      e =Sqrt(eq)
	esq=(a*a/(b*b))-1.
c
      e1=(1.-b/a)/(1+b/a)
	f1=(1.-eq/4.-3.*e4/64.-5.*e6/256.)
	f2=3.*eq/8.+3.*e4/32.+45.*e6/1024.
	f3=15.*e4/256+45.*e6/1024.
	f4=35.*e6/3072.
	phi=x/a/f1
	betf=phi+(3.*e1/2.-27*e1*e1*e1/32.)*Sin(2.*phi)
     . +(21*e1*e1/16.-55*e1*e1*e1*e1/32.)*Sin(4.*phi)
     . +151*e1*e1*e1/96.*Sin(6.*phi)
	sb=Sin(betf)
	cb=Cos(betf)
	tb=Tan(betf)
	tb2=tb*tb
	tb4=tb2*tb2
	sb2=sb*sb
	c=esq*cb*cb
	xn=a/Sqrt(1.-eq*sb2)
	xm=a*(1.-eq)/(1.-eq*sb2)/Sqrt(1.-eq*sb2)
	d=(y-500000.)/xn
	d2=d*d
      bet=betf-xn*tb/xm*d2*(.5-(5+3*tb2+c-9*tb2*c)*d2/24
     .        +(61+90*tb2+45*tb4)*d2*d2/720.)
	alf=alf0+1./cb*d*(1.-(1+2*tb2+c)*d2+
     .    (5+28*tb2+6*c+8*tb2*c+24*tb4)*d2*d2/120.)
      betd=bet*180/pi
	alfd=alf*180/pi
	End

