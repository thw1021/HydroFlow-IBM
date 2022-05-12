    
!!############################################################################!!                              
!!#    Program      : ARCHIVE_BIN                                             !!  
!!#    Version      : 2.0 (2019-12-26 21:18:55)                               !!    
!!#    Aothorized by: Wang Jian                                               !!
!!#    Organization : SHANGHAI JIAO TONG UNIVERSITY                           !!           
!!#    Address      : SHANGHAI, CHINA                                         !!   
!!#    Contact      : ***@****.com                                            !!
!!#--------------------------------------------------------------------------#!!
!!#  Description: This subroutine is used to save data as binary file         !!                                   
!!#                                                                           !!
!!#    input files : ' '                                                      !!
!!#    output files: 'VEL_TRAM.BIN','C_TRAM.BIN'                              !!                                            
!!############################################################################!!                          
      
	Subroutine ARCHIVE_BIN
      Include './Include/OCERM_INF' 
      DOUBLE PRECISION TIMECOUNT

C-----   SAVING THE COMPUTIONAL RESULTS OF THE FLUID FIELD  ----------------C
	  If(Mod(NSTEP,ISAVE) .EQ. 0.0 .OR. NSTEP .EQ. ISTART) Then
    	    If(NSTEP .EQ. ISTART) Then
	      TIMECOUNT = TIME * 24. + HOURSTAR
	      If(KBM .GT. 1) Then
			 IJK = IJP * KB
	         IJCELL = IJM * KBM
	         WRITE(ITF_BIN) TIMECOUNT, IJK, IJCELL
               IF (WAQ .NE. 'NEGLECT')THEN
               WRITE(ITC_BIN) TIMECOUNT, IJK, IJCELL
               ENDIF
               
	         Do K = 1, KB
	            Do I = 1, IJP
	               ZP = HP(I)+ELFV(I)
	               If(ZP .LE. 0.0) Then
				      ZP = 0.001
					  ZZZ = Z(K) * ZP + (ZP - HP(I))
				   Else
				      ZZZ = Z(K)*(HP(I)+ELFV(I))+ELFV(I)
				   Endif	       
	               Write(ITF_BIN) 
     &	 	         PXY(I,1), PXY(I,2), ZZZ, ELFV(I), PNV(I,K),
     &               URV(I,K), VRV(I,K), WRV(I,K), VORTICITYV(I,K),
     &               VORXV(I,K), VORYV(I,K), VORZV(I,K),OMIGAV(I,K), 
     &               KMV(I,K),TKERV(I,K),TDISSRV(I,K)
                     
                     IF (WAQ .NE. 'NEGLECT')THEN
                     Write(ITC_BIN) 
     &	 	         PXY(I,1), PXY(I,2), ZZZ, CRV(I,K)
                     ENDIF                     
	            Enddo
	         Enddo
	         Do K = 1, KBM 
	            Do I = 1, IJM
	               Write(ITF_BIN) 
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
              IF(WAQ .NE. 'NEGLECT')THEN
	         Do K = 1, KBM 
	            Do I = 1, IJM
	               Write(ITC_BIN) 
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
              ENDIF 
C---------------------                                          
	      Else
C-----      2D
	         WRITE(ITF_BIN) TIMECOUNT, IJP, IJM
	         Do I = 1, IJP
	            Write(ITF_BIN) 
     &	 	       PXY(I,1), PXY(I,2), ELFV(I), PNV(I,1),
     &               URV(I,1), VRV(I,1), KMV(I,1)
	         Enddo
	         Do I = 1, IJM
	            Write(ITF_BIN) CELL_INDEX(I,1), CELL_INDEX(I,2), 
     &                           CELL_INDEX(I,3), CELL_INDEX(I,4) 
	         Enddo
	      Endif
c-------------zone>=2--------------------------------------------------
	    Else
	      TIMECOUNT = TIME * 24. + HOURSTAR
	      If(KBM .GT. 1) Then
		    IJK = IJP * KB
	        IJCELL = IJM * KBM
	        WRITE(ITF_BIN) TIMECOUNT, IJK, IJCELL
              IF (WAQ .NE. 'NEGLECT')THEN
              WRITE(ITC_BIN) TIMECOUNT, IJK, IJCELL
              ENDIF    
	        Do K = 1, KB
	           Do I = 1, IJP
		          ZP = HP(I)+ELFV(I)
                    If(ZP .LE. 0.0) Then
				     ZP = 0.001
			  	     ZZZ = Z(K) * ZP + (ZP - HP(I))
			  	  Else
				     ZZZ = Z(K)*(HP(I)+ELFV(I))+ELFV(I)
				  Endif	       
	              Write(ITF_BIN) 
     &	 	         PXY(I,1), PXY(I,2), ZZZ, ELFV(I), PNV(I,K),
     &               URV(I,K), VRV(I,K), WRV(I,K), VORTICITYV(I,K),
     &               VORXV(I,K), VORYV(I,K), VORZV(I,K),OMIGAV(I,K), 
     &               KMV(I,K),TKERV(I,K),TDISSRV(I,K)
                    
                     IF (WAQ .NE. 'NEGLECT')THEN
                     Write(ITC_BIN) 
     &	 	         PXY(I,1), PXY(I,2), ZZZ, CRV(I,K)
                     ENDIF                                          
	           Enddo
	        Enddo
	      Else
C-----      2D
	        WRITE(ITF_BIN)TIMECOUNT, IJP, IJM
	        Do I = 1, IJP
	           Write(ITF_BIN) 
     &	 	         PXY(I,1), PXY(I,2), ELFV(I), PNV(I,1),
     &                 URV(I,1), VRV(I,1), KMV(I,1)
	        Enddo
	      Endif
	    Endif
      Endif
           
      
	Return
	End
