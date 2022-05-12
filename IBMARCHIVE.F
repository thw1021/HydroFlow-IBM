C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (12/05/2022)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #     
C                    Save result data of IBM module                         #                                           
C                                                                           #
C############################################################################
      Subroutine IBMARCHIVE
      
      Include './Include/OCERM_INF'
      
      Double precision Cp1,Cp2,Cp3,FD,FL,CD,CL
      Double precision UREF,ELFREF,PNREF
      Integer SAVE_ID(IJM,KBM)
      Integer IJK,IJCELL
      Character*10 FN
      
      
C===========================================================================C      
C         Save results of time averaged
C===========================================================================C
      
      Do I = 1,IJP
          SAE(I) = SAE(I) + ELFV(I)
          Do K = 1,KB
              SAU(I,K) = SAU(I,K) + URV(I,K)
              SAV(I,K) = SAV(I,K) + VRV(I,K)
              SAW(I,K) = SAW(I,K) + WRV(I,K)
              SAP(I,K) = SAP(I,K) + PNV(I,K)
          Enddo
      Enddo

      If (NSTEP.EQ.IEND) Then
          Do I = 1,IJP
              ELFV(I) = SAE(I) / (IEND - ISTART)
              Do K = 1,KB
	            URV(I,K) = SAU(I,K) / (IEND - ISTART)
 	            VRV(I,K) = SAV(I,K) / (IEND - ISTART)
	            WRV(I,K) = SAW(I,K) / (IEND - ISTART)
                  PNV(I,K) = SAP(I,K) / (IEND - ISTART)
              Enddo
          Enddo
          
          Write(FN,'(F10.4)')  TIME * 24. + HOURSTAR
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
     &             PXY(I,1), PXY(I,2), ZZZ, ELFV(I), PNV(I,K),
     &             URV(I,K), VRV(I,K), WRV(I,K), VORTICITYV(I,K),
     &             VORXV(I,K), VORYV(I,K), VORZV(I,K),OMIGAV(I,K), 
     &             QCRIV(I,K),KMV(I,K),TKEV(I,K),TDISSV(I,K)
          Enddo
          Enddo    
      Endif
      
C===========================================================================C      
C         Calculate and save Cp with THETA 
C===========================================================================C
      
      Do I = 1,IJM
          SAE2(I) = SAE2(I) + ELF(I)
          Do K = 1,KBM
              SAP2(I,K) = SAP2(I,K) + PN(I,K)
          Enddo
      Enddo
              
      If (NSTEP.EQ.IEND) Then
          Do I = 1,IJM
              SAE2(I) = SAE2(I) / (IEND - ISTART)
              Do K = 1,KBM
	            SAP2(I,K) = SAP2(I,K) / (IEND - ISTART)
              Enddo
          Enddo
          
          Cp1 = 0
          Cp2 = 0
          Cp3 = 0
          Do K = 1,KBM
          Do I = 1,IJM
          If (IDC(I,K).GT.1.9) Then
              ELFREF = ELF(INFCELL)
              UREF = U(INFCELL,K) / (HC(INFCELL) + ELFREF)
              PNREF = PN(INFCELL,K)
              Cp1 = (ELF(I) - ELFREF) * 2 * GRAV / (UREF**2)
              Cp2 = (PN(I,K) - PNREF) / (500 * UREF**2)
              Cp3 = Cp1 + Cp2
              Write(IBMCP,'(7F16.8)') CXY(I,1), CXY(I,2), ZZ(k),
     &                                Cp1, Cp2, Cp3, ELFREF
          Endif
          Enddo
          Enddo
      Endif

C===========================================================================C      
C         Calculate and save CD,CL with time 
C===========================================================================C

      If ((ISKILL.NE.0).AND.(Mod(NSTEP,ISKILL).EQ.0) 
     &	 .AND.(NSTEP.GT.ISTART)) Then
          TMIDDLE = TIME - DTI * DAYI / SKILLI / 2.
          TSAVE = (TMIDDLE * 24. + HOURSTAR) * 3600.
          FD = 0
          FL = 0
          Do I = 1,IJM
              If (IDC(I,K).GT.1.9) Then
                  FD = FD - BFX(I,K) * AREA(I)/(HC(I)+ELF(I))
                  FL = FL - BFY(I,K) * AREA(I)/(HC(I)+ELF(I))
              Endif
          Enddo
          CD = 20*FD/(UINF**2)
          CL = 20*FL/(UINF**2)
          
          Write(IBMCDCL,'(4F16.10)') TSAVE,CD,CL,UINF
      Endif     
          
C===========================================================================C      
C         Save results of ghost cells
C===========================================================================C
          
      Do K = 1,KBM
      Do I = 1,IJM
          SAVE_ID(I,K) = 0
      Enddo
      Enddo
      
      Do K = 1,KBM
      Do I = 1,IJM
          If (IDC(I,K).GT.0.5) Then
          SAVE_ID(I,K) = 1
          Do J = 1,CELL_POLYGEN(I)
              If (IDC(CELL_SIDE(I,J,2),K).LT.0.5) SAVE_ID(I,K) = 1
          Enddo
          Endif
      Enddo
      Enddo
      
      Do K = 1,KBM
      Do I = 1,IJM
          If (SAVE_ID(I,K).EQ.1) Then
          Do J = 1,CELL_POLYGEN(I)
              If (IDC(CELL_SIDE(I,J,2),K).LT.0.5) SAVE_ID(I,K) = 1
          Enddo
          Endif
      Enddo
      Enddo
      
      Do K = 1,KBM
          Do I = 1,IBM
              If (SAVE_ID(I,K).LT.0.5) Then
                  URV(CELL_SIDE(I,J,3),K+1) = 0
                  VRV(CELL_SIDE(I,J,3),K+1) = 0
                  WRV(CELL_SIDE(I,J,3),K+1) = 0
                  ELFV(CELL_SIDE(I,J,3)) = 0
                  PNV(CELL_SIDE(I,J,3),K+1) = 0  
                  TKEV(CELL_SIDE(I,J,3),K+1) = 0
                  TDISSV(CELL_SIDE(I,J,3),K+1) = 0
              Endif
          Enddo
      Enddo

      If((Mod(NSTEP,ISAVE).EQ.0).OR.(NSTEP.EQ.ISTART)) Then
          If(NSTEP.EQ.ISTART) Then
              Write(FN,'(F10.3)')  TIME * 24. + HOURSTAR
              IJK = IJP * KB
              IJCELL = IJM * KBM
              Write(IBMDATA,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20)')
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
                  Write(IBMDATA,'(10E16.8)') 
     &                PXY(I,1), PXY(I,2), ZZZ, ELFV(I), PNV(I,K),
     &                URV(I,K), VRV(I,K), WRV(I,K),TKEV(I,K),TDISSV(I,K)
              Enddo
              Enddo
              
              Do K = 1, KBM 
              Do I = 1, IJM
                  Write(IBMDATA,'(8I8)') 
     &                IJP * (K-1) + CELL_INDEX(I,1), 
     &                IJP * (K-1) + CELL_INDEX(I,2), 
     &                IJP * (K-1) + CELL_INDEX(I,3), 
     &                IJP * (K-1) + CELL_INDEX(I,4), 
     &                IJP * (K) + CELL_INDEX(I,1), 
     &                IJP * (K) + CELL_INDEX(I,2), 
     &                IJP * (K) + CELL_INDEX(I,3), 
     &                IJP * (K) + CELL_INDEX(I,4)
              Enddo   
              Enddo
	    Else
	        Write(FN,'(F10.4)')  TIME * 24. + HOURSTAR
		    IJK = IJP * KB
	        IJCELL = IJM * KBM
	        Write(IBMDATA,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20,A20)')
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
                  Write(IBMDATA,'(10E16.8)') 
     &                PXY(I,1), PXY(I,2), ZZZ, ELFV(I), PNV(I,K),
     &                URV(I,K), VRV(I,K), WRV(I,K),TKEV(I,K),TDISSV(I,K)
              Enddo
	        Enddo
	    Endif
      Endif
      
      End Subroutine