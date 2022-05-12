C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (12/05/2022)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #     
C                     Visualization of data results                         #                                           
C                                                                           #
C############################################################################
      Subroutine IBMVISUAL(FIGMOD,FRAMENUM)
      
      Include './Include/OCERM_INF'
      
      Dimension ZPHC(IJM,KBM),ZPHP(IJP,KB)
      Integer FIGMOD,FRAMENUM
      
      
      Do K = 1,KBM
      Do I = 1,IJM
          ZPHC(I,K) = ZZ(K)*(HC(I)+ELF(I))+ELF(I)
      Enddo
      Enddo 
      
      Do K = 1,KB
      Do I = 1,IJP
          ZPHP(I,K) = Z(K)*(HP(I)+ELFV(I))+ELFV(I)
      Enddo
      Enddo    
      
C===========================================================================C     
C         Visualization of ghost cells  
C===========================================================================C  
      
      If ((FIGMOD.EQ.1).AND.(IBMSCHEME.EQ.0)) Then
          Do nn = 1,NUMIB
          Do J = 1,IBN(nn)   
              Write(IBMB,1000) BXY(nn,J,1),BXY(nn,J,2),BXY(nn,J,3)
          Enddo
          Enddo

          Do K = 1,KBM
          Do I = 1,IJM
              If ((IDC(I,K).LT.0.1)) Then
                  Write(IBMF,1000) CXY(I,1),CXY(I,2),ZPHC(I,K)
              Else if ((IDC(I,K).EQ.1)) Then
                  Write(IBMS,1000) CXY(I,1),CXY(I,2),ZPHC(I,K)
              Else If (IDC(I,K).GT.1.9) Then
                  Write(IBMGC,1000) CXY(I,1),CXY(I,2),ZPHC(I,K)
              Endif
          Enddo
          Enddo
          
          If (IBGCSCHEME.GT.0) Then
          near_points = 28
          Do N = 1,NUMGH
              Write(IBMGCI,1000) PGC(N,1),PGC(N,2),PGC(N,3)
              Do J = 1,near_points
              If (IDC(NGC(N,J,1),NGC(N,J,2)).LT.0.9) Then
                  DK = (PGC(N,1) - CXY(NGC(N,J,1),1)) ** 2
     &                +(PGC(N,2) - CXY(NGC(N,J,1),2)) ** 2
                  If (DK.LE.RITP**2) Then
                  Write(IBMP,1000) CXY(NGC(N,J,1),1),CXY(NGC(N,J,1),2),
     &                             ZPHC(NGC(N,J,1),NGC(N,J,2))
                  Endif
              Endif
              Enddo
          Enddo
          Endif
          
      Endif
      
C===========================================================================C      
C         Visualization of ghost points      
C===========================================================================C  

      If ((FIGMOD.EQ.2).AND.(IBMSCHEME.EQ.0)) Then        
          Do nn = 1,NUMIB
          Do J = 1,IBN(nn)   
              Write(IBMB,1000) BXY(nn,J,1),BXY(nn,J,2),BXY(nn,J,3)
          Enddo
          Enddo
      
          Do K = 1,KB
          Do I = 1,IJP
              If (IDP(I,K).LT.0.5) Then
                  Write(IBMF,1000) PXY(I,1),PXY(I,2),ZPHP(I,K)
              Else if (IDP(I,K).GT.0.5) Then
                  Write(IBMS,1000) PXY(I,1),PXY(I,2),ZPHP(I,K)
              Endif
          Enddo
          Enddo
      Endif

C===========================================================================C
C         Visualization of IBM mesh
C===========================================================================C

      If ((FIGMOD.EQ.3).AND.(IBMSCHEME.EQ.0)) Then
          Do nn = 1,NUMIB
          Do J = 1,IBN(nn)   
              Write(IBMB,1000) BXY(nn,J,1),BXY(nn,J,2),BXY(nn,J,3)
          Enddo
          Enddo

          Do I = 1,IJM
              If((CXY(I,1).GT.-0.1).AND.(CXY(I,1).LT.0.1).AND.
     &           (CXY(I,2).GT.-0.1).AND.(CXY(I,2).LT.0.1)) Then
          Do J = 1,CELL_POLYGEN(I)
              Write(IBMF,1000) PXY(CELL_SIDE(I,J,3),1),
     &                         PXY(CELL_SIDE(I,J,3),2),ZPHC(I,1)
          Enddo
          Write(IBMF,1000) PXY(CELL_SIDE(I,1,3),1),
     &                     PXY(CELL_SIDE(I,1,3),2),ZPHC(I,1)
          Endif
          Enddo
          
          K = 5
          Do I = 1,IJM
              If ((IDC(I,K).EQ.1)) Then
              Do J = 1,CELL_POLYGEN(I)
                  Write(IBMS,1000) PXY(CELL_SIDE(I,J,3),1),
     &                             PXY(CELL_SIDE(I,J,3),2),ZPHC(I,K)
              Enddo
              Else If (IDC(I,K).GT.1.9) Then
              Do J = 1,CELL_POLYGEN(I)
                  Write(IBMGC,1000) PXY(CELL_SIDE(I,J,3),1),
     &                              PXY(CELL_SIDE(I,J,3),2),ZPHC(I,K)
              Enddo
              Endif
          Enddo
      Endif
          
C===========================================================================C
C         Visualization of moving boundary
C===========================================================================C
      If ((SAVESTEP.GT.1).AND.(IBMSCHEME.GE.1)) Then 
      If ((Mod(NSTEP,FRAMENUM) .EQ. 0.0 ).OR.(NSTEP.EQ.ISTART)) Then
          Do nn = 1,NUMIB
          Do J = 1,IBN(nn)   
             ! Write(IBMB,1000) BXY(nn,J,1),BXY(nn,J,2),BXY(nn,J,3)
          Enddo
          Enddo

          Do K = 1,KBM
          Do I = 1,IJM
              If ((IDC(I,K).LT.0.1)) Then
                  Write(IBMF,1000) CXY(I,1),CXY(I,2),ZPHC(I,K)
              Else if ((IDC(I,K).EQ.1)) Then
                  Write(IBMS,1000) CXY(I,1),CXY(I,2),ZPHC(I,K)
              Else If (IDC(I,K).GT.1.9) Then
                  Write(IBMGC,1000) CXY(I,1),CXY(I,2),ZPHC(I,K)
              Endif
          Enddo
          Enddo
          
          Write(IBMB,1100) 0,0,1
          Write(IBMS,1100) 0,0,1
          Write(IBMGC,1100) 0,0,1
      Endif
      Endif
      
1000  Format (3F16.8) 
1100  Format (3I16)
      
      
      End Subroutine