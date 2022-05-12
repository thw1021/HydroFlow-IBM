C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (12/05/2022)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #     
C               Ghost Cell Immersed Boundary Method moudle                  #                                           
C                                                                           #
C############################################################################
      Subroutine IBMALGC
      
      Include './Include/OCERM_INF'
      
      Real(kind=8) BFU(IJM,KBM),BFV(IJM,KBM),BFW(IJM,KBM)
      Real(kind=8) BFP(IJM,KBM),BFE(IJM)
      Double precision UTEMP,VTEMP,WTEMP,PTEMP,ETEMP,YTEMP,ETAK,WM
      Double precision RM,DK,UCOS,USIN,UFR,UPT,UGT,UGN,UPL,YPL,REST
      Integer near_points,NUM

      
C===========================================================================C      
C         Initialize the arrays      
C===========================================================================C
          
      Do I = 1,IJM
          Do K = 1,KBM
              BFU(I,K) = 0
              BFV(I,K) = 0
              BFW(I,K) = 0
              BFP(I,K) = 0
              BFX(I,K) = 0
              BFY(I,K) = 0
              BFZ(I,K) = 0
          Enddo
      Enddo
              
C===========================================================================C      
C         Direct assignment method £¨Static boundary£©
C===========================================================================C
              
      If (IBGCSCHEME.EQ.0) Then 
      If (NSTEP.EQ.ISTART) Write(*,*) 'GCMODE is direct assignment.'
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      
      Do N = 1,NUMGH
          BFU(IGC(N,1),IGC(N,2)) = 0
          BFV(IGC(N,1),IGC(N,2)) = 0
          BFW(IGC(N,1),IGC(N,2)) = 0
      Enddo
      
      Do N = 1,NUMSL
          BFU(ISC(N,1),ISC(N,2)) = 0
          BFV(ISC(N,1),ISC(N,2)) = 0
          BFW(ISC(N,1),ISC(N,2)) = 0
      Enddo
      
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      Endif

C===========================================================================C
C         Linear weight interpolation
C===========================================================================C 
      
      If (IBGCSCHEME.EQ.1) Then
      If (NSTEP.EQ.ISTART) Then
          Write(*,*) 'GCMODE is linear weight interpolation.'
      Endif
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      
      near_points = 28

      Do N = 1,NUMGH
          ETAK = 0
              
c          DK = (PGC(N,1) - BXY(1,PGC(N,4),1)) ** 2
c     &        +(PGC(N,2) - BXY(1,PGC(N,4),2)) ** 2
          Do J = 1,near_points
          If (IDC(NGC(N,J,1),NGC(N,J,2)).LT.0.9) Then
              DK = (PGC(N,1) - CXY(NGC(N,J,1),1)) ** 2
     &            +(PGC(N,2) - CXY(NGC(N,J,1),2)) ** 2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      ETAK = ETAK + 1/DK
                  Else
                      ETAK = ETAK + 0
                  Endif
              Endif
          Endif
          Enddo
          
          WM = 0
          UTEMP = 0
          VTEMP = 0
          WTEMP = 0
          PTEMP = 0
          ETEMP = 0
          
          Do J = 1,near_points
          If (IDC(NGC(N,J,1),NGC(N,J,2)).LT.0.9) Then
              DK = (PGC(N,1) - CXY(NGC(N,J,1),1)) ** 2
     &            +(PGC(N,2) - CXY(NGC(N,J,1),2)) ** 2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      WM = (1/ETAK) * (1/DK)
                  Else
                      WM = 0
                  Endif
                  UTEMP = UTEMP + WM * U(NGC(N,J,1),NGC(N,J,2))
                  VTEMP = VTEMP + WM * V(NGC(N,J,1),NGC(N,J,2))
                  WTEMP = WTEMP + WM * W(NGC(N,J,1),NGC(N,J,2))
                  ETEMP = ETEMP + WM * ELF(NGC(N,J,1))
                  PTEMP = PTEMP + WM * PN(NGC(N,J,1),NGC(N,J,2))
              Endif
          Endif
          Enddo
              
          If (IBGCMODE.EQ.0) Then
              BFU(IGC(N,1),IGC(N,2)) = -UTEMP
              BFV(IGC(N,1),IGC(N,2)) = -VTEMP
              BFW(IGC(N,1),IGC(N,2)) = -WTEMP
          Else if (IBGCMODE.EQ.1) Then
              BFU(IGC(N,1),IGC(N,2)) = 0.5 * UTEMP
              BFV(IGC(N,1),IGC(N,2)) = 0.5 * VTEMP
              BFW(IGC(N,1),IGC(N,2)) = 0.5 * WTEMP
          Endif
          BFE(IGC(N,1)) = ETEMP
          BFP(IGC(N,1),IGC(N,2)) = PTEMP
      Enddo
      
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Endif

C===========================================================================C
C         Wall function interpolation
C===========================================================================C 

      If (IBGCSCHEME.EQ.2) Then
      If (NSTEP.EQ.ISTART) Then
          Write(*,*) 'GCMODE is wall function interpolation.'
      Endif
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      
      near_points = 28

      Do N = 1,NUMGH
          ETAK = 0
          
          Do J = 1,near_points
              If (IDC(NGC(N,J,1),NGC(N,J,2)).LT.0.9) Then
                  DK = (PGC(N,1) - CXY(NGC(N,J,1),1)) ** 2
     &                +(PGC(N,2) - CXY(NGC(N,J,1),2)) ** 2
                  If (DK.LE.RITP**2) Then
                      If (DK.GT.1.0E-8) Then
                          ETAK = ETAK + 1/DK
                      Else
                          ETAK = ETAK + 0
                      Endif
                  Endif
              Endif
          Enddo
          
          WM = 0
          UTEMP = 0
          VTEMP = 0
          WTEMP = 0
          UGN = 0
          UGT = 0
          
          Do J = 1,near_points
          If (IDC(NGC(N,J,1),NGC(N,J,2)).LT.0.9) Then
              DK = (PGC(N,1) - CXY(NGC(N,J,1),1)) ** 2
     &            +(PGC(N,2) - CXY(NGC(N,J,1),2)) ** 2
              If (DK.LE.RITP**2) Then 
                  If (DK.GT.1.0E-10) Then
                      WM = (1/ETAK) * (1/DK)
                  Else
                      WM = 0
                  Endif
                  UTEMP = UTEMP + WM * U(NGC(N,J,1),NGC(N,J,2))
                  VTEMP = VTEMP + WM * V(NGC(N,J,1),NGC(N,J,2))
                  WTEMP = WTEMP + WM * QZ(NGC(N,J,1),NGC(N,J,2))
              If (NSTEP.EQ.ISTART) Then
              Write(IBMP,1000) CXY(NGC(N,J,1),1),CXY(NGC(N,J,1),2),
     &                         ZZ(NGC(N,J,2))                  
              Endif
              Endif
          Endif
          Enddo
          
          
          UPT = UTEMP * PGC(N,8) - VTEMP * PGC(N,7)
          
          NUM = 0
	    REST = 100
	    YPL = 20
	  
	    Do while (REST.GT.1E-10) 
              YTEMP = YPL
		    YPL = (YPL + 0.4*PGC(N,6)*Abs(UPT)/UMOL)/(1+log(9.8*YPL))
		    REST = Abs(YPL - YTEMP) 
		    NUM = NUM + 1
		    If (NUM.GT.1000) REST = 0
          End do
          
	    UFR = YPL * UMOL / PGC(N,6)
	    UPL = 2.5 * log(9.8*0.5*PGC(N,6)*UFR / UMOL)
          
          If (UPT.GE.0) Then 
	        UGT = UPL * UFR
          Else
              UGT = -1 * UPL * UFR
          Endif
          
	    UGN = 0
	    TKE(IGC(N,1),IGC(N,2)) = UFR**2/0.3
	    TDISS(IGC(N,1),IGC(N,2)) = Sqrt((0.00008/(PGC(N,6)**2))**2+
     &            (Sqrt(TKE(IGC(N,1),IGC(N,2)))/(0.219*PGC(N,6)))**2)

c          UGT = 0.5 * (UTEMP * PGC(N,8) - VTEMP * PGC(N,7)) 
c          UGN = 0.5 * (UTEMP * PGC(N,7) + VTEMP * PGC(N,8))    

          BFU(IGC(N,1),IGC(N,2)) = UGN * PGC(N,7) + UGT * PGC(N,8)
          BFV(IGC(N,1),IGC(N,2)) = UGN * PGC(N,8) - UGT * PGC(N,7)

          If (NSTEP.EQ.20000) Then
              If (UGN.LT.0) Then
                  UCOS = -1 * PGC(N,7)
                  USIN = -1 * PGC(N,8)
              Else
                  UCOS = PGC(N,7)
                  USIN = PGC(N,8)
              Endif
              Write(IBMN,1100) CXY(IGC(N,1),1),CXY(IGC(N,1),2),UCOS,USIN
                      
              If (UGT.LT.0) Then
                  UCOS = PGC(N,7)
                  USIN = -1 * PGC(N,8)
              Else
                  UCOS = -1 * PGC(N,7)
                  USIN = PGC(N,8)
              Endif
              Write(IBMT,1100) CXY(IGC(N,1),1),CXY(IGC(N,1),2),USIN,UCOS
          Endif
      Enddo
      
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Endif
      
C===========================================================================C      
C         Assign value to solid cells and ghost cells
C===========================================================================C
            
      Do K = 1, KBM
      Do I = 1, IJM
          If (IDC(I,K).GT.0.9) Then
              BFX(I,K) = (BFU(I,K) - U(I,K)) / DTI
              BFY(I,K) = (BFV(I,K) - V(I,K)) / DTI
              BFZ(I,K) = (BFW(I,K) - W(I,K)) / DTI           
          Endif
      Enddo
      Enddo
      
      Do K = 1, KBM
      Do I = 1, IJM
          If (IDC(I,K).GT.0.9) Then
              U(I,K) = BFU(I,K)
              V(I,K) = BFV(I,K)
              QZ(I,K) = 0
              W(I,K) = 0
              ELF(I) = BFE(I)
              PN(I,K) = BFP(I,K)
              If(IDC(I,K).EQ.1) Then
                  TKE(I,K) = 0
                  TDISS(I,K) = 100000
                  ELF(I) = 0
              Endif  
          Endif
      Enddo
      Enddo

1000  Format (3F8.4)
1100  Format (4F10.6)
      
      
      End Subroutine