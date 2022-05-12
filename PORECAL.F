c-----------------------------------------------------------------------------------c
C                               CHEN JUN, 2019.09.25                                #
C                               SHANGHAI JIAO TONG UNIVERSITY                       #
C                               SHANGHAI, CHINA                                     #
c-----------------------------------------------------------------------------------#
c                     根据格心孔隙率计算水平和垂向单元面上的孔隙率                     #
c-----------------------------------------------------------------------------------c
      Subroutine PORECAL
      Include './Include/OCERM_INF'
      INTEGER I, J, K, IL, IR
      !Double Precision VMESH, VPATICLE
      !Print*,'Calculating PORE'
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I, J, K, IL, IR) 
      Do K = 1 , KBM
!$OMP DO        
        Do J = 1 , IJE
            IL = INDEX_EDGE(J,K,1)
	      IR = INDEX_EDGE(J,K,2)
            IF(IL .EQ. -999)THEN
                PORE_HF(J,K) = PORE(IR,K)
            ELSEIF(IR .EQ. -999)THEN
                PORE_HF(J,K) = PORE(IL,K)
            ELSE                
                PORE_HF(J,K)= 0.5 * (PORE(IL,K) + PORE(IR,K))
            ENDIF
            
            !PRINT*,J,K,PORE_HF(J,K)
         Enddo
!$OMP END DO                    
      Enddo
      
      Do K = 1 , KB
!$OMP DO        
        Do I = 1 , IJM  
            IF(K .EQ. 1)THEN
               PORE_VF(I,K) = PORE(I,K)
            ELSEIF(K .EQ. KB)THEN
               PORE_VF(I,K) = PORE(I,K-1) 
            ELSE    
               PORE_VF(I,K) = 0.5 * (PORE(I,K-1) + PORE(I,K))          
            ENDIF
            !PRINT*,I,K,PORE_VF(I,K)
         Enddo
!$OMP END DO                    
      Enddo 
      
!$OMP END PARALLEL 
      !PRINT*,MINVAL(PORE_HF),MINVAL(PORE_VF)
      Return
      End 