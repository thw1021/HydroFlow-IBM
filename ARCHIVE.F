C############################################################################
C                                                                           #
C               SAVING THE COMPUTATIONAL RESULTS                            #
C                 VERSION 1.0 (10/08/2010)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C############################################################################
	Subroutine ARCHIVE(DAYI,HOURSTAR)
      Include './Include/OCERM_INF'

	Dimension UINS(200,100),VINS(200,100),WINS(200,100),PINS(200,100),
     &            VORXINS(200,100), VORYINS(200,100), VORZINS(200,100),
     &            TXINS(200), ELFINS(200), VISINS(200,100) 	 
      Dimension K_GAUGE(NGAUGE),R_GAUGE(NGAUGE)
      Character*10 FN

C============================================================================C
C         Saving results of fixed observation points                                      
C============================================================================C

      Do N = 1, EPTS
          I = INXE(N)
c         ESAVE(N) = ESAVE(N) + EL(I) * SKILLI 
          ESAVE(N) = EL(I)
          Do K = 1, KBM
              USAVE(N,K) = USAVE(N,K) + UR(I,K) * SKILLI
	        VSAVE(N,K) = VSAVE(N,K) + VR(I,K) * SKILLI
	        WSAVE(N,K) = WSAVE(N,K) + WR(I,K) * SKILLI
          Enddo
      Enddo

      If(ISKILL .NE. 0 .AND. Mod(NSTEP,ISKILL) .EQ. 0 
     &	 .AND. NSTEP .GT. ISTART) Then
      TMIDDLE = TIME - DTI * DAYI / SKILLI / 2.
           
      Write(IUTSR,5600) (TMIDDLE * 24. + HOURSTAR) * 3600.,
C     &	   (ESAVE(N), (USAVE(N,K), K = 1, KBM), WSAVE(N,K), 
C     &       N = 1, EPTS)
     &        (ESAVE(N), N = 1, EPTS) 

C----------------- output at constant height ------------------- 

      If (IVEL_INS == 1) Then
          Do N = 1, NGAUGE
              ZZ_GAUGE = (-0.5*HC(1) - ELF(INDEX_GAUGE(N))) /
     &                   (HC(INDEX_GAUGE(N)) + ELF(INDEX_GAUGE(N)) ) 
              Do K = 1,KBM-1
              If(ZZ(K).GE.ZZ_GAUGE .AND. ZZ(K+1).LE.ZZ_GAUGE) THEN
                  K_GAUGE(N)= K
	            R_GAUGE(N) =(ZZ(K)-ZZ_GAUGE)/(ZZ(K)-ZZ(K+1))
              Endif
              Enddo
          Enddo

C----------------- output at middle layer ------------------- 

      Else if (IVEL_INS == 2) Then          
          Do K = 1,KBM-1
              If(ZZ(K).GE.-0.50 .AND. ZZ(K+1).LE.-0.50) Then
                  K_GAUGE(1) = K
	            R_GAUGE(1) = (ZZ(K)+0.50)/(ZZ(K)-ZZ(K+1))
              Endif
          Enddo	 
          Do N = 2, NGAUGE
              K_GAUGE(N) = K_GAUGE(1)
	        R_GAUGE(N) = R_GAUGE(1)
          Enddo
      Endif      
                  
C----------------- output at gauge points ------------------- 
                  
      Write(21,5600) (TMIDDLE * 24. + HOURSTAR) * 3600.,
     & (USAVE(N,K_GAUGE(1))+(USAVE(N,K_GAUGE(1)+1)-USAVE(N,K_GAUGE(1)))*
     &  R_GAUGE(1), N = 1, EPTS) 

      
      If (IVEL_INS .NE. 0) Then              
           write(IUINS,5600)(TMIDDLE * 24. + HOURSTAR) * 3600.,
     & (UR(INDEX_GAUGE(N),K_GAUGE(N)) +(UR(INDEX_GAUGE(N),K_GAUGE(N)+1)
     & -UR(INDEX_GAUGE(N),K_GAUGE(N)))* R_GAUGE(N), N = 1, NGAUGE)
           write(IVINS,5600)(TMIDDLE * 24. + HOURSTAR) * 3600.,
     & (VR(INDEX_GAUGE(N),K_GAUGE(N)) +(VR(INDEX_GAUGE(N),K_GAUGE(N)+1)
     & -VR(INDEX_GAUGE(N),K_GAUGE(N)))* R_GAUGE(N), N = 1, NGAUGE)
           write(IWINS,5600)(TMIDDLE * 24. + HOURSTAR) * 3600.,
     & (WR(INDEX_GAUGE(N),K_GAUGE(N)) +(WR(INDEX_GAUGE(N),K_GAUGE(N)+1)
     & -WR(INDEX_GAUGE(N),K_GAUGE(N)))* R_GAUGE(N), N = 1, NGAUGE)
      Endif
		 
      Do N = 1, EPTS
          ESAVE(N) = 0.0
          Do K = 1, KBM
              USAVE(N,K) = 0.0
              VSAVE(N,K) = 0.0
              WSAVE(N,K) = 0.0
          Enddo
      Enddo
      Endif
      
C============================================================================C
C         Saving results of fluid field                                      
C============================================================================C
      
	  If(Mod(NSTEP,ISAVE) .EQ. 0.0 .OR. NSTEP .EQ. ISTART) Then
    	    If(NSTEP .EQ. ISTART) Then
	      Write(FN,'(F10.3)')  TIME * 24. + HOURSTAR
	      If(KBM .GT. 1) Then
			 IJK = IJP * KB
	         IJCELL = IJM * KBM
	         WRITE(ITF,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20)')
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
	               Write(ITF,'(20E16.5)') 
     &	 	         PXY(I,1), PXY(I,2), ZZZ, ELFV(I), PNV(I,K),
     &               URV(I,K), VRV(I,K), WRV(I,K), VORTICITYV(I,K),
     &               VORXV(I,K), VORYV(I,K), VORZV(I,K),OMIGAV(I,K), 
     &               QCRIV(I,K),KMV(I,K),TKEV(I,K),TDISSV(I,K)
	            Enddo
	         Enddo
	         Do K = 1, KBM 
	            Do I = 1, IJM
	               Write(ITF,'(8I8)') 
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
	         WRITE(ITF,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20)')
     &        'ZONE T=', '"',FN,'"',' I=',IJP,' K=',IJM,',F=FEPOINT',
     &        ', ET=QUADRILATERAL'
	         Do I = 1, IJP
	            Write(ITF,'(20E16.5)') 
     &	 	       PXY(I,1), PXY(I,2), ELFV(I), PNV(I,1),
     &               URV(I,1), VRV(I,1), KMV(I,1)
	         Enddo
	         Do I = 1, IJM
	            Write(ITF,'(8I8)') CELL_INDEX(I,1), CELL_INDEX(I,2), 
     &                               CELL_INDEX(I,3), CELL_INDEX(I,4) 
	         Enddo
	      Endif
	    Else
	      Write(FN,'(F10.4)')  TIME * 24. + HOURSTAR
	      If(KBM .GT. 1) Then
		    IJK = IJP * KB
	        IJCELL = IJM * KBM
	        WRITE(ITF,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20,A20)')
     &        'ZONE T=', '"',FN,'"',' I=',IJK,' K=',IJCELL,',F=FEPOINT',
     &        ', ET=BRICK',', D=(FECONNECT)'
	        Do K = 1, KB
	           Do I = 1, IJP
		          ZP = HP(I)+ELFV(I)
                    If(ZP .LE. 0.0) Then
				     ZP = 0.001
			  	     ZZZ = Z(K) * ZP + (ZP - HP(I))
			  	  Else
				     ZZZ = Z(K)*(HP(I)+ELFV(I))+ELFV(I)
				  Endif	       
	              Write(ITF,'(20E16.5)') 
     &	 	         PXY(I,1), PXY(I,2), ZZZ, ELFV(I), PNV(I,K),
     &               URV(I,K), VRV(I,K), WRV(I,K), VORTICITYV(I,K),
     &               VORXV(I,K), VORYV(I,K), VORZV(I,K),OMIGAV(I,K), 
     &               QCRIV(I,K),KMV(I,K),TKEV(I,K),TDISSV(I,K)
	           Enddo
	        Enddo
	      Else
C-----      2D
	        WRITE(ITF,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20,A20)')
     &        'ZONE T=', '"',FN,'"',' I=',IJP,' K=',IJM,',F=FEPOINT',
     &        ', ET=QUADRILATERAL',', D=(FECONNECT)'
	        Do I = 1, IJP
	           Write(ITF,'(20E16.4)') 
     &	 	         PXY(I,1), PXY(I,2), ELFV(I), PNV(I,1),
     &                 URV(I,1), VRV(I,1), KMV(I,1)
	        Enddo
	      Endif
	    Endif
      Endif

C============================================================================C
C         Temporal saving for hot start                                       
C============================================================================C

      If(Mod(NSTEP,IHOTSTART).EQ.0.0) Then
	     Write(FN,'(F10.4)') NSTEP*DTI/3600.
	     Open(IUSTA,Form='unformatted',File='./Restart/'//FN//'.sta',
     &	      Status='unknown')
	     Write (IUSTA) NSTEP, DZR, Z,ZZ,DZ, DZZ, HC, HS, HP, DC, DS, 
     &         EL, ELFV, ELFM, COR, WTSURF, WTBOT, AAM, KM, KH, DT,
     &         KQ, U, V, W, UR, VR, WR, UV, VV, WV, UN, VN, CBC, VIS, 
     &         T, SED, RMEAN, Z01, PN, QZ, CCM, CFM, INDEX_EDGE,
     &         THOUR, TBX, TBY, TKE, TDISS, STRESS,
     &         VIS, TENSOR, VORTICITY,C,QZ,GRADP
C     *        NNDEM,TTTDEM, XPDEM,YPDEM,ZPDEM, 
C     *        UPDEM,VPDEM,WPDEM,OMXDEM,OMYDEM,OMZDEM,
C     *        PORE,GRADP, 
C     *        DELTAT_PP1X,DELTAT_PP1Y,DELTAT_PP1Z,N_P2,N_PP2 
	       Close(IUSTA)

           !#hot restart for averaging#, WangJian, 2020-3-6 19:36:24 !    
           IF(IAVESAVE==1)THEN
           Open(IAVESTA,Form='unformatted',File='./Restart/AVE.STA',
     &	      Status='unknown')
	     Write (IAVESTA) NSTEP-1, AVEURV, AVEVRV,AVEWRV, AVEUTRV,
     &                  RYNUU, RYNUV, RYNUW, RYNVV, RYNVW,RYNWW  
	     Close(IAVESTA)
           ENDIF
	  Endif

        If(DEM .NE. 'NEGLECT')THEN
          IF( (IDEMRES .EQ. 0 .AND. MOD(NSTEP,IHOTSTART).EQ.0) .OR.
     *        (IDEMRES .NE. 0 .AND. MOD(NNDEM,IDEMRES).EQ.0) )THEN
            WRITE(FN,'(F10.4)') NNDEM*TDEM/3600.
            Open(IDEMSTA,Form='unformatted',File='./Restart/'//FN//
     *	       'DEM.sta',Status='unknown')
	      Write (IDEMSTA) NNDEM,TTTDEM,XPDEM,YPDEM,ZPDEM,
     *            UPDEM,VPDEM,WPDEM,OMXDEM,OMYDEM,OMZDEM,
     *            PORE,DELTAT_PP1X,DELTAT_PP1Y,DELTAT_PP1Z,
     *            N_P2,N_PP2,UPSUMDEM,VPSUMDEM,WPSUMDEM,
     *            UFPDEM_PRE,VFPDEM_PRE,WFPDEM_PRE,
     *            UFDEM_PRE,VFDEM_PRE,WFDEM_PRE,
     *            N_SYNFA,DEL_HIS,NSYM_CONT  
            Close(IDEMSTA)            
          ENDIF 
      ENDIF
      
	Return
 5600 Format (F15.5,1000E15.5)
	End
