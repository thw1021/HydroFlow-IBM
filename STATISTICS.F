C############################################################################
C                                                                           #
C                     VARIABLE STATISTICS                                   #
C                   BY CHEN JUN(21/04/2019)                                 #
C                                                                           #
C############################################################################
      
      Subroutine STATISTICS(DAYI,HOURSTAR)
      Include './Include/OCERM_INF'
      Character*10 FN

	If(   NSTEP .GE. (NSTART+ISTART) 
     &	 .AND. ISKILL .NE. 0  .AND. Mod(NSTEP,ISKILL) .EQ. 0) Then
         
           TMIDDLE = TIME - DTI * DAYI / SKILLI / 2.       
           IAVENUM = (NSTEP-(ISTART+NSTART))/ISKILL   
C          WRITE(*,*) IAVENUM    
C          T = TMIDDLE * 24. + HOURSTAR) * 3600.           
      IF(IAVEHOT .NE. 1 .AND. IAVENUM.EQ.0) THEN
          Do K = 1, KB
              Do I = 1, IJP
                  UTRV(IJP,KB) = 0.0
                  AVEURV(IJP,KB) = 0.0
                  AVEVRV(IJP,KB) = 0.0
                  AVEWRV(IJP,KB) = 0.0
                  AVEUTRV(IJP,KB) = 0.0
                  AVEUURV(IJP,KB) = 0.0
                  AVEUVRV(IJP,KB) = 0.0
                  AVEUWRV(IJP,KB) = 0.0
                  AVEVVRV(IJP,KB) = 0.0
                  AVEVWRV(IJP,KB) = 0.0
                  AVEWWRV(IJP,KB) = 0.0
                  RYNUU(IJP,KB) = 0.0
                  RYNUV(IJP,KB) = 0.0
                  RYNUW(IJP,KB) = 0.0
                  RYNVV(IJP,KB) = 0.0
                  RYNVW(IJP,KB) = 0.0
                  RYNWW(IJP,KB) = 0.0
              ENDDO
          ENDDO
      ENDIF                            
      
      IF(IAVEVEL .EQ. 1) THEN
          Do K = 1, KB
              Do I = 1, IJP
                  UTRV(I,K) = SQRT(URV(I,K)*URV(I,K)+VRV(I,K)*VRV(I,K)+
     &                WRV(I,K)*WRV(I,K))
                  AVEURV(I,K) = ( AVEURV(I,K) * IAVENUM + URV(I,K))
     &                / ( IAVENUM +1 )
                  AVEVRV(I,K) = ( AVEVRV(I,K) * IAVENUM + VRV(I,K))
     &                /( IAVENUM +1 )
                  AVEWRV(I,K) = ( AVEWRV(I,K) * IAVENUM + WRV(I,K))
     &                /( IAVENUM +1 )
                  AVEUTRV(I,K) = ( AVEUTRV(I,K) * IAVENUM + UTRV(I,K))
     &                /( IAVENUM +1 )
              ENDDO
          ENDDO
      END IF
      
      IF(IAVERYN .EQ. 1) THEN
          Do K = 1, KB
              Do I = 1, IJP
                  AVEUURV(I,K) = (AVEUURV(I,K)* IAVENUM +
     &                URV(I,K)*URV(I,K))/(IAVENUM+1)
                  AVEUVRV(I,K) = (AVEUVRV(I,K)* IAVENUM +
     &                URV(I,K)*VRV(I,K))/(IAVENUM+1)
                  AVEUWRV(I,K) = (AVEUWRV(I,K)* IAVENUM +
     &                URV(I,K)*WRV(I,K))/(IAVENUM+1)             
                  AVEVVRV(I,K) = (AVEVVRV(I,K)* IAVENUM +
     &                VRV(I,K)*VRV(I,K))/(IAVENUM+1)
                  AVEVWRV(I,K) = (AVEVWRV(I,K)* IAVENUM +
     &                VRV(I,K)*WRV(I,K))/(IAVENUM+1)
                  AVEWWRV(I,K) = (AVEWWRV(I,K)* IAVENUM +
     &                WRV(I,K)*WRV(I,K))/(IAVENUM+1)
                  
                  RYNUU(I,K) = AVEUURV(I,K)-AVEURV(I,K)*AVEURV(I,K)
                  RYNUV(I,K) = AVEUVRV(I,K)-AVEURV(I,K)*AVEVRV(I,K)
                  RYNUW(I,K) = AVEUWRV(I,K)-AVEURV(I,K)*AVEWRV(I,K)
                  RYNVV(I,K) = AVEVVRV(I,K)-AVEVRV(I,K)*AVEURV(I,K)
                  RYNVW(I,K) = AVEVWRV(I,K)-AVEVRV(I,K)*AVEWRV(I,K)
                  RYNWW(I,K) = AVEWWRV(I,K)-AVEWRV(I,K)*AVEWRV(I,K)
              ENDDO
          ENDDO
      END IF
      
      ENDIF
      
      If(  NSTEP .GE. (NSTART+ISTART).AND. 
     &     Mod(NSTEP,ISAVE) .EQ. 0.0   ) Then
          REWIND(IAVEDAT)
          READ(IAVEDAT,*)
          READ(IAVEDAT,*)
          Write(FN,'(F10.4)')  TIME * 24. + HOURSTAR
          If(KBM .GT. 1) Then
		    IJK = IJP * KB
	        IJCELL = IJM * KBM
	         WRITE(IAVEDAT,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20)')
     &        'ZONE T=', '"',FN,'"',' I=',IJK,' K=',IJCELL,',F=FEPOINT',
     &        ', ET=BRICK'
	        Do K = 1, KB
	           Do I = 1, IJP
		          ZP = HP(I)+ELFV(I)
                    If(ZP .LE. 0.0) Then
				     ZP = 0.001
			  	     ZZZ = Z(K) * ZP + (ZP - HP(I))
			  	  Else
				     ZZZ = Z(K)*(HP(I)+ELFV(I))+ELFV(I)
				  Endif	       
	              Write(IAVEDAT,'(20E16.5)') 
     &	 	         PXY(I,1), PXY(I,2), ZZZ, AVEURV(I,K), AVEVRV(I,K),
     &               AVEWRV(I,K), AVEUTRV(I,K), RYNUU(I,K), RYNUV(I,K),
     &               RYNUW(I,K), RYNVV(I,K), RYNVW(I,K),RYNWW(I,K)
	           Enddo
              Enddo
              Do K = 1, KBM 
	            Do I = 1, IJM
	               Write(IAVEDAT,'(8I8)') 
     &                 IJP * (K-1) + CELL_INDEX(I,1), 
     &                 IJP * (K-1) + CELL_INDEX(I,2), 
     &                 IJP * (K-1) + CELL_INDEX(I,3), 
     &                 IJP * (K-1) + CELL_INDEX(I,4), 
     &                 IJP * (K) + CELL_INDEX(I,1), 
     &                 IJP * (K) + CELL_INDEX(I,2), 
     &                 IJP * (K) + CELL_INDEX(I,3), 
     &                 IJP * (K) + CELL_INDEX(I,4)
	            Enddo   
	         Enddo
	      Else
C-----      2D
	        WRITE(IAVEDAT,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20,A20)')
     &        'ZONE T=', '"',FN,'"',' I=',IJP,' K=',IJM,',F=FEPOINT',
     &        ', ET=QUADRILATERAL',', D=(FECONNECT)'
	        Do I = 1, IJP
	           Write(IAVEDAT,'(20E16.4)') 
     &	 	         PXY(I,1), PXY(I,2), AVEURV(I,K), AVEVRV(I,K)
              Enddo
              Do I = 1, IJM
	            Write(IAVEDAT,'(8I8)') CELL_INDEX(I,1), CELL_INDEX(I,2), 
     &                               CELL_INDEX(I,3), CELL_INDEX(I,4) 
              Enddo
          EndIF
      ENDIF
      
      RETURN
      END
      
      