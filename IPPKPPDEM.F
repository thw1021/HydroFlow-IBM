c-----------------------------------------------------------------------------------c
c                             SUBROUTINE PROGRAM                                    #
C                             DESINGED BY ChenJun                                   #
C                                                 SHANGHAI JIAO TONG UNIVERSITY     #
C                                                 SHANGHAI, CHINA                   #
c-----------------------------------------------------------------------------------#
c 子程序功能：
c 1------颗粒流场网格寻址及颗粒存储      
c-----------------------------------------------------------------------------------#

      Subroutine IPPKPPDEM
      Include './Include/OCERM_INF'
      Parameter(NSCHEME=2)
      Double Precision X1DEM,X2DEM,Y1DEM,Y2DEM,DIJDEM,XPGDEM
      Double Precision R1,RTOL,ELFPDEM,DCPDEM,HPDEM
      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,IC,X1DEM,X2DEM,Y1DEM,Y2DEM,
!$OMP&         DIJDEM,XPGDEM,R1,RTOL,ELFPDEM,DCPDEM,HPDEM) 
      IF(PSMOD .EQ. 'MOD1')THEN
	  Do K = 1, KBM
!$OMP DO	
	    Do I = 1, IJM
            NCP(I,K) = 0
          Enddo
!$OMP END DO 
        Enddo
      ENDIF
      
!$OMP DO     
      Do I = 1 , I_PARTICLE         

       IF(NSYM_CONT(I) .EQ. 1)THEN          
          
c----颗粒流场网格定位-水平向  
        DO J = 1, (NL-1)*NL*2+1
            
          IC = INE(IPP(I),J) 
          
          Do K = 1, Cell_POLYGEN( IC )  
                
            X1DEM = PXY(CELL_SIDE(IC,K,3),1)
            Y1DEM = PXY(CELL_SIDE(IC,K,3),2)
            X2DEM = PXY(CELL_SIDE(IC,K,4),1)
            Y2DEM = PXY(CELL_SIDE(IC,K,4),2)          
                    
            XPGDEM = ( XPDEM(I) - X1DEM ) * ( Y2DEM - Y1DEM ) -
     &               ( YPDEM(I) - Y1DEM ) * ( X2DEM - X1DEM )
                         
            If(XPGDEM .GT. 0.) Then
              Exit
            Endif  
            
          Enddo 
        
          If(K .EQ. ( Cell_POLYGEN(IC) + 1 ) ) Then     
            IPP(I) = IC
            Exit
          ENDIF
        
        ENDDO
        
        IF(J .EQ. (NL-1)*NL*2+2)THEN
          Print*, 'Error1 of IPP in IPPKPPDEM.F!'
          Print*, 'I = ',I
          Print*, 'XPDEM = ',XPDEM(I)
          Print*, 'YPDEM = ',YPDEM(I)
          Print*, 'ZPDEM = ',ZPDEM(I)           
          Pause 
        ENDIF                      
      
c----颗粒流场网格定位-垂向 
        IF(NSCHEME .EQ. 1)THEN
        !插值计算颗粒水平投影点水位、静水深和总水深
	  RTOL = 0.0
	  ELFPDEM = 0.0
        HPDEM = 0.0
c        DCPDEM = 0.0
        
	  Do J = 1, CELL_POLYGEN(IPP(I))              
c	    R1 = Sqrt((XPDEM(I) - PXY( CELL_INDEX(IPP(I),J) , 1)) ** 2. +
c     &	          (YPDEM(I) - PXY( CELL_INDEX(IPP(I),J) , 2)) ** 2.)   !距离作为权重
	    R1 = (XPDEM(I) - PXY( CELL_INDEX(IPP(I),J) , 1)) ** 2. +
     &	     (YPDEM(I) - PXY( CELL_INDEX(IPP(I),J) , 2)) ** 2.   !距离平方作为权重  
          
          RTOL = RTOL + 1. / R1
          
          ELFPDEM = ELFPDEM + ELFV( CELL_INDEX(IPP(I),J) ) / R1
          
c          HPDEM = HPDEM + HP(CELL_INDEX(IPP(I),J)) / R1 
          
          DCPDEM = DCPDEM + ( HP(CELL_INDEX(IPP(I),J)) + 
     &             ELFV(CELL_INDEX(IPP(I),J)) ) / R1 
	  Enddo
	  ELFPDEM = ELFPDEM / RTOL
        HPDEM = HPDEM / RTOL
        DCPDEM = DCPDEM / RTOL    
      
        !SIGMA_P(I) = ( ZPDEM(I) - ELF(IPP(I)) ) / DC(IPP(I))
        SIGMA_P(I) = ( ZPDEM(I) - ELFPDEM ) / DCPDEM
        !SIGMA_P(I) = ( ZPDEM(I) - ELF(IPP(I)) ) / (ELF(IPP(I))+HPDEM) 
        
        ELSEIF(NSCHEME .EQ. 2)THEN

         HPDEM = ( DEMCOSB(IPP(I),1) * (XPDEM(I) - 
     &             PXY( CELL_INDEX(IPP(I),1) , 1)) +
     &             DEMCOSB(IPP(I),2) * (YPDEM(I) - 
     &             PXY( CELL_INDEX(IPP(I),1) , 2)) ) /
     &             DEMCOSB(IPP(I),3) + HP(CELL_INDEX(IPP(I),1))    
         
         ELFPDEM = - ( DEMCOSA(IPP(I),1) * (XPDEM(I) - 
     &               PXY( CELL_INDEX(IPP(I),1) , 1)) +
     &               DEMCOSA(IPP(I),2) * (YPDEM(I) - 
     &               PXY( CELL_INDEX(IPP(I),1) , 2)) ) /
     &               DEMCOSA(IPP(I),3) + ELFV(CELL_INDEX(IPP(I),1))  
         
         SIGMA_P(I) = ( ZPDEM(I) - ELFPDEM ) / ( HPDEM + ELFPDEM ) 
         
        ELSEIF(NSCHEME .EQ. 3)THEN       

         HPDEM = 0.5 * (
     &           ( DEMCOSB(IPP(I),1) * (XPDEM(I) - 
     &             PXY( CELL_INDEX(IPP(I),1) , 1)) +
     &             DEMCOSB(IPP(I),2) * (YPDEM(I) - 
     &             PXY( CELL_INDEX(IPP(I),1) , 2)) ) /
     &             DEMCOSB(IPP(I),3) + HP(CELL_INDEX(IPP(I),1))  
     &           + 
     &           ( DEMCOSB(IPP(I),4) * (XPDEM(I) - 
     &             PXY( CELL_INDEX(IPP(I),1) , 1)) +
     &             DEMCOSB(IPP(I),5) * (YPDEM(I) - 
     &             PXY( CELL_INDEX(IPP(I),1) , 2)) ) /
     &             DEMCOSB(IPP(I),6) + HP(CELL_INDEX(IPP(I),1)) 
     &           )
         
         ELFPDEM = 0.5 * ( ELFV(CELL_INDEX(IPP(I),1)) -
     &             ( DEMCOSA(IPP(I),1) * (XPDEM(I) - 
     &               PXY( CELL_INDEX(IPP(I),1) , 1)) +
     &               DEMCOSA(IPP(I),2) * (YPDEM(I) - 
     &               PXY( CELL_INDEX(IPP(I),1) , 2)) ) /
     &               DEMCOSA(IPP(I),3)    
     &            +  ELFV(CELL_INDEX(IPP(I),1)) -  
     &             ( DEMCOSA(IPP(I),4) * (XPDEM(I) - 
     &               PXY( CELL_INDEX(IPP(I),1) , 1)) +
     &               DEMCOSA(IPP(I),5) * (YPDEM(I) - 
     &               PXY( CELL_INDEX(IPP(I),1) , 2)) ) /
     &               DEMCOSA(IPP(I),6)  
     &             )  
         
         SIGMA_P(I) = ( ZPDEM(I) - ELFPDEM ) / ( HPDEM + ELFPDEM )
         
        ENDIF 
        
        IF(SIGMA_P(I) .LT. -1.0) THEN
		Write (IUPRT,*),XPDEM(I),YPDEM(I),ZPDEM(I)
 	    RTOL = 0.0
	    ELFPDEM = 0.0
          DCPDEM = 0.0
          HPDEM = 0.0
		Write (IUPRT,*)
	    Do J = 1, CELL_POLYGEN(IPP(I)) 
 
!	      R1 = Sqrt((XPDEM(I) - PXY( CELL_INDEX(IPP(I),J) , 1)) ** 2. +

!     &	            (YPDEM(I) - PXY( CELL_INDEX(IPP(I),J) , 2)) ** 2.) 
	      R1 = (XPDEM(I) - PXY( CELL_INDEX(IPP(I),J) , 1)) ** 2. +
     &	       (YPDEM(I) - PXY( CELL_INDEX(IPP(I),J) , 2)) ** 2. 
            RTOL = RTOL + 1. / R1
            ELFPDEM = ELFPDEM + ELFV( CELL_INDEX(IPP(I),J) ) / R1
            DCPDEM = DCPDEM + ( HP(CELL_INDEX(IPP(I),J)) + 
     &             ELFV(CELL_INDEX(IPP(I),J)) ) / R1 
		  HPDEM = HPDEM + HP(CELL_INDEX(IPP(I),J)) / R1 
            Write (IUPRT,*),PXY( CELL_INDEX(IPP(I),J) , 1),
     &	                  PXY( CELL_INDEX(IPP(I),J) , 2), 
     &	                  ELFV( CELL_INDEX(IPP(I),J) ),
     &	                  HP(CELL_INDEX(IPP(I),J)) 
	    Enddo
	    ELFPDEM = ELFPDEM / RTOL
          DCPDEM = DCPDEM / RTOL 
          HPDEM = HPDEM / RTOL
          Write (IUPRT,*)
		Write (IUPRT,*),ELFPDEM,DCPDEM,HPDEM 
          PRINT*,'Error2 of KPP in IPPKPPDEM.F!'
          PRINT*,'SIGMA_P(I) = ',SIGMA_P(I)          
          PAUSE
        ELSEIF(SIGMA_P(I) .GT. 0.0)THEN            
          KPP(I) = 0          
        ELSEIF(KPP(I) .GT. 1 .AND. KPP(I) .LT. KBM)THEN            
          IF(SIGMA_P(I) .LE. Z(KPP(I)) .AND. 
     &         SIGMA_P(I) .GE. Z(KPP(I)+1))THEN
            KPP(I) = KPP(I) 
          ELSEIF(SIGMA_P(I) .LE. Z(KPP(I)-1) .AND. 
     &         SIGMA_P(I) .GE. Z(KPP(I)))THEN 
            KPP(I) = KPP(I)-1
          ELSEIF(SIGMA_P(I) .LE. Z(KPP(I)+1) .AND. 
     &         SIGMA_P(I) .GE. Z(KPP(I)+2))THEN 
            KPP(I) = KPP(I)+1 
          ELSE
            !PRINT*,'Error1 of KPP in IPPKPPDEM.F!'
            !PRINT*,'SIGMA_P(I) = ',SIGMA_P(I)
            !print*,I,ZPDEM(I),ELF(IPP(I)),DC(IPP(I))
            !PAUSE   
            Do K = 1 , KBM
             If ( SIGMA_P(I) .LE. Z(K) .AND. 
     &        SIGMA_P(I) .GE. Z(K+1) ) Then
              KPP(I) = K           
              exit 
             Endif
            Enddo 
            IF(K .EQ. KBM+1)THEN
              PAUSE
              PRINT*,'Error3 of KPP in IPPKPPDEM.F!'
              PRINT*,'SIGMA_P(I) = ',SIGMA_P(I)
            ENDIF         
          ENDIF
        ELSEIF(KPP(I) .EQ. 1)THEN                    
          IF(SIGMA_P(I) .LE. Z(1) .AND. 
     &         SIGMA_P(I) .GE. Z(2))THEN
            KPP(I) = 1 
          ELSEIF(SIGMA_P(I) .LE. Z(2) .AND. 
     &         SIGMA_P(I) .GE. Z(3))THEN 
            KPP(I) = 2
          ELSEIF(SIGMA_P(I) .GT. Z(1))THEN
            KPP(I) = 0
          ELSE  
            !PRINT*,'Error2 of KPP in IPPKPPDEM.F!'
            !PRINT*,'SIGMA_P(I) = ',SIGMA_P(I)
            !PAUSE 
            Do K = 1 , KBM
             If ( SIGMA_P(I) .LE. Z(K) .AND. 
     &        SIGMA_P(I) .GE. Z(K+1) ) Then
              KPP(I) = K           
              exit 
             Endif
            Enddo 
            IF(K .EQ. KBM+1)THEN
              PAUSE
              PRINT*,'Error4 of KPP in IPPKPPDEM.F!'
              PRINT*,'SIGMA_P(I) = ',SIGMA_P(I)
            ENDIF              
          ENDIF
        ELSEIF(KPP(I) .EQ. KBM)THEN           
          IF(SIGMA_P(I) .LE. Z(KBM) .AND. 
     &         SIGMA_P(I) .GE. Z(KBM+1))THEN
            KPP(I) = KBM 
          ELSEIF(SIGMA_P(I) .LE. Z(KBM-1) .AND. 
     &         SIGMA_P(I) .GE. Z(KBM))THEN 
            KPP(I) = KBM-1
          ELSE
            !PRINT*,'Error3 of KPP in IPPKPPDEM.F!'
            !PRINT*,'SIGMA_P(I) = ',SIGMA_P(I)
            !PAUSE      
            Do K = 1 , KBM
             If ( SIGMA_P(I) .LE. Z(K) .AND. 
     &        SIGMA_P(I) .GE. Z(K+1) ) Then
              KPP(I) = K           
              exit 
             Endif
            Enddo 
            IF(K .EQ. KBM+1)THEN
              PAUSE
              PRINT*,'Error5 of KPP in IPPKPPDEM.F!'
              PRINT*,'SIGMA_P(I) = ',SIGMA_P(I)
            ENDIF       
          ENDIF       
        ELSEIF(KPP(I) .EQ. 0)THEN           
          IF(SIGMA_P(I) .LE. Z(1) .AND. 
     &         SIGMA_P(I) .GE. Z(2))THEN
            KPP(I) = 1
          ELSEIF(SIGMA_P(I) .GT. 0.0)THEN  
            KPP(I) = 0 
          ELSE
            !PRINT*,'Error4 of KPP in IPKPDEM.F!'
            !PRINT*,'SIGMA_P(I) = ',SIGMA_P(I)
            !PAUSE
            Do K = 1 , KBM
             If ( SIGMA_P(I) .LE. Z(K) .AND. 
     &        SIGMA_P(I) .GE. Z(K+1) ) Then
              KPP(I) = K           
              exit 
             Endif
            Enddo 
            IF(K .EQ. KBM+1)THEN
              PAUSE
              PRINT*,'Error6 of KPP in IPKPDEM.F!'
              PRINT*,'SIGMA_P(I) = ',SIGMA_P(I)
            ENDIF            
          ENDIF          
        ENDIF      

C- 流场网格搜索颗粒存储
        IF(PSMOD .EQ. 'MOD1')THEN  
          NCP(IPP(I),KPP(I)) = NCP(IPP(I),KPP(I)) + 1 
          IF(NCP(IPP(I),KPP(I)) .GT. ICPMAX)THEN
            PRINT*,'Error7 of ICPMAX in IPPKPPDEM.F!'
            PRINT*,'NCP=',NCP(IPP(I),KPP(I))
            PAUSE
          ENDIF
          ICP(IPP(I),KPP(I),NCP(IPP(I),KPP(I))) = I
        ENDIF        
       
       ENDIF 
        
      Enddo
      
!$OMP END DO
!$OMP END PARALLEL                       
      Return 
      End
