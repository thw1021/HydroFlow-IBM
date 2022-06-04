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
      Double precision UTEMP,VTEMP,WTEMP,PTEMP,ETEMP,YTEMP
      Double precision UTEMP2,VTEMP2,WTEMP2,PTEMP2,ETEMP2,YTEMP2
      Double precision ETAK,WM,DK,ZALPHA
      Double precision UCOS,USIN,UFR,UPT,UGT,UGN,UPL,YPL,REST
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
      
      Do K = 1,KBM
      Do I = 1,IJM
          ZRC(I,K) = ZZ(K)*(HC(I)+ELF(I))+ELF(I)
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
      
      If (IBMDOMAIN.EQ.2) Then
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      
      near_points = 28

      Do N = 1,NUMGH
      If (IDC(IGC(N,1),IGC(N,2)).EQ.2) Then  
          ETAK = 0
          Do J = 1,near_points
          I = NBR(IGC(N,1),J)
          K = PGC(N,4)
          If (IDC(I,K).LT.0.9) Then
              DK = (PGC(N,1)-CXY(I,1))**2 + (PGC(N,2)-CXY(I,2))**2
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
          I = NBR(IGC(N,1),J)
          K = PGC(N,4)
          If (IDC(I,K).LT.0.9) Then
              DK = (PGC(N,1)-CXY(I,1))**2 + (PGC(N,2)-CXY(I,2))**2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      WM = (1/ETAK) * (1/DK)
                  Else
                      WM = 0
                  Endif
                  UTEMP = UTEMP + WM * U(I,K)
                  VTEMP = VTEMP + WM * V(I,K)
                  WTEMP = WTEMP + WM * W(I,K)
                  ETEMP = ETEMP + WM * ELF(I)
                  PTEMP = PTEMP + WM * PN(I,K)
              Endif
          Endif
          Enddo
          
          BFU(IGC(N,1),IGC(N,2)) = 0.5 * UTEMP
          BFV(IGC(N,1),IGC(N,2)) = 0.5 * VTEMP
          BFW(IGC(N,1),IGC(N,2)) = 0.5 * WTEMP
          BFE(IGC(N,1)) = ETEMP
          BFP(IGC(N,1),IGC(N,2)) = PTEMP
      Else if (IDC(IGC(N,1),IGC(N,2)).GE.3) Then
          ZALPHA = 1 - Abs(PGC(N,3)-ZRC(IGC(N,1),IGC(N,2))) / PGC(N,7)
          
          BFU(IGC(N,1),IGC(N,2)) = ZALPHA * U(I,K)
          BFV(IGC(N,1),IGC(N,2)) = ZALPHA * V(I,K)
          BFW(IGC(N,1),IGC(N,2)) = ZALPHA * W(I,K)
          BFE(IGC(N,1)) = ELF(I)
          BFP(IGC(N,1),IGC(N,2)) = PN(I,K)
      Endif
      Enddo
      
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Endif
      
      If (IBMDOMAIN.EQ.3) Then
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      near_points = 28

      Do N = 1,NUMGH
! Step1: K layer distance weight interpolation
          ETAK = 0
          Do J = 1,near_points
          I = NBR(IGC(N,1),J)
          K = PGC(N,4)
          If (IDC(I,K).LT.0.9) Then
              DK = (PGC(N,1)-CXY(I,1))**2 + (PGC(N,2)-CXY(I,2))**2
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
          I = NBR(IGC(N,1),J)
          K = PGC(N,4)
          If (IDC(I,K).LT.0.9) Then
              DK = (PGC(N,1)-CXY(I,1))**2 + (PGC(N,2)-CXY(I,2))**2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      WM = (1/ETAK) * (1/DK)
                  Else
                      WM = 0
                  Endif
                  UTEMP = UTEMP + WM * U(I,K)
                  VTEMP = VTEMP + WM * V(I,K)
                  WTEMP = WTEMP + WM * W(I,K)
                  ETEMP = ETEMP + WM * ELF(I)
                  PTEMP = PTEMP + WM * PN(I,K)
              Endif
          Endif
          Enddo
          
! Step2: K-1 layer distance weight interpolation          
          ETAK = 0
          Do J = 1,near_points
          I = NBR(IGC(N,1),J)
          K = PGC(N,4) - 1
          If (IDC(I,K).LT.0.9) Then
              DK = (PGC(N,1)-CXY(I,1))**2 + (PGC(N,2)-CXY(I,2))**2
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
          UTEMP2 = 0
          VTEMP2 = 0
          WTEMP2 = 0
          PTEMP2 = 0
          ETEMP2 = 0
          
          Do J = 1,near_points
          I = NBR(IGC(N,1),J)
          K = PGC(N,4) - 1
          If (IDC(I,K).LT.0.9) Then
              DK = (PGC(N,1)-CXY(I,1))**2 + (PGC(N,2)-CXY(I,2))**2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      WM = (1/ETAK) * (1/DK)
                  Else
                      WM = 0
                  Endif
                  UTEMP2 = UTEMP2 + WM * U(I,K)
                  VTEMP2 = VTEMP2 + WM * V(I,K)
                  WTEMP2 = WTEMP2 + WM * W(I,K)
                  PTEMP2 = PTEMP2 + WM * PN(I,K)
              Endif
          Endif
          Enddo
          
! Step3: K & K-1 layer linear interpolation  
          I = IGC(N,1)
          K = PGC(N,4)
          ZALPHA = Abs(ZRC(I,K-1)-PGC(N,3)) / Abs(ZRC(I,K-1)-ZRC(I,K)) 
                 
          BFU(IGC(N,1),IGC(N,2)) = ZALPHA * UTEMP + (1-ZALPHA) * UTEMP2
          BFV(IGC(N,1),IGC(N,2)) = ZALPHA * VTEMP + (1-ZALPHA) * VTEMP2
          BFW(IGC(N,1),IGC(N,2)) = ZALPHA * WTEMP + (1-ZALPHA) * WTEMP2
          BFE(IGC(N,1)) = ETEMP
          BFP(IGC(N,1),IGC(N,2)) = ZALPHA * PTEMP + (1-ZALPHA) * PTEMP2
      Enddo
      
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Endif

      Endif

C===========================================================================C
C         Wall function interpolation
C===========================================================================C 

      If (IBGCSCHEME.EQ.2) Then
      If (NSTEP.EQ.ISTART) Then
          Write(*,*) 'GCMODE is wall function interpolation.'
      Endif
      
      If (IBMDOMAIN.EQ.2) Then
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      
      near_points = 28

      Do N = 1,NUMGH
      If (IDC(IGC(N,1),IGC(N,2)).EQ.2) Then 
          ETAK = 0
          Do J = 1,near_points
          I = NBR(IGC(N,1),J)
          K = PGC(N,4)
          If (IDC(I,K).LT.0.9) Then
              DK = (PGC(N,1)-CXY(I,1))**2 + (PGC(N,2)-CXY(I,2))**2
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
          I = NBR(IGC(N,1),J)
          K = PGC(N,4)
          If (IDC(I,K).LT.0.9) Then
              DK = (PGC(N,1)-CXY(I,1))**2 + (PGC(N,2)-CXY(I,2))**2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      WM = (1/ETAK) * (1/DK)
                  Else
                      WM = 0
                  Endif
                  UTEMP = UTEMP + WM * U(I,K)
                  VTEMP = VTEMP + WM * V(I,K)
                  WTEMP = WTEMP + WM * W(I,K)
                  ETEMP = ETEMP + WM * ELF(I)
                  PTEMP = PTEMP + WM * PN(I,K)
              Endif
          Endif
          Enddo
      Else if (IDC(IGC(N,1),IGC(N,2)).GE.3) Then
          ZALPHA = 1 - Abs(PGC(N,3)-ZRC(IGC(N,1),IGC(N,2))) / PGC(N,7)
          
          UTEMP = ZALPHA * U(I,K)
          VTEMP = ZALPHA * V(I,K)
          WTEMP = ZALPHA * W(I,K)
          ETEMP = ELF(I)
          PTEMP = PN(I,K)
      Endif
          
      UPT = UTEMP * PGC(N,9) - VTEMP * PGC(N,8)
          
      NUM = 0
      REST = 100
      YPL = 20
      Do while (REST.GT.1E-10) 
          YTEMP = YPL
          YPL = (YPL + 0.4*PGC(N,7)*Abs(UPT)/UMOL)/(1+log(9.8*YPL))
          REST = Abs(YPL - YTEMP) 
          NUM = NUM + 1
          If (NUM.GT.1000) REST = 0
      End do
          
      UFR = YPL * UMOL / PGC(N,7)
      UPL = 2.5 * log(9.8*0.5*PGC(N,7)*UFR / UMOL)
          
      If (UPT.GE.0) Then 
          UGT = UPL * UFR
      Else
          UGT = -1 * UPL * UFR
      Endif
          
      UGN = 0
      TKE(IGC(N,1),IGC(N,2)) = UFR**2/0.3
      TDISS(IGC(N,1),IGC(N,2)) = Sqrt((0.00008/(PGC(N,7)**2))**2+
     &            (Sqrt(TKE(IGC(N,1),IGC(N,2)))/(0.219*PGC(N,7)))**2)

c     UGT = 0.5 * (UTEMP * PGC(N,8) - VTEMP * PGC(N,7)) 
c     UGN = 0.5 * (UTEMP * PGC(N,7) + VTEMP * PGC(N,8))    

      BFU(IGC(N,1),IGC(N,2)) = UGN * PGC(N,8) + UGT * PGC(N,9)
      BFV(IGC(N,1),IGC(N,2)) = UGN * PGC(N,9) - UGT * PGC(N,8)
      BFW(IGC(N,1),IGC(N,2)) = 0.5 * WTEMP
      BFE(IGC(N,1)) = ETEMP
      BFP(IGC(N,1),IGC(N,2)) = PTEMP

      If (NSTEP.EQ.IEND) Then
          If (UGN.LT.0) Then
              UCOS = -1 * PGC(N,8)
              USIN = -1 * PGC(N,9)
          Else
              UCOS = PGC(N,8)
              USIN = PGC(N,9)
          Endif
          Write(IBMN,1100) CXY(IGC(N,1),1),CXY(IGC(N,1),2),UCOS,USIN
                      
          If (UGT.LT.0) Then
              UCOS = PGC(N,8)
              USIN = -1 * PGC(N,9)
          Else
              UCOS = -1 * PGC(N,8)
              USIN = PGC(N,9)
          Endif
          Write(IBMT,1100) CXY(IGC(N,1),1),CXY(IGC(N,1),2),USIN,UCOS
      Endif
      Enddo
      
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Endif
      
      If (IBMDOMAIN.EQ.3) Then
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      near_points = 28

      Do N = 1,NUMGH
! Step1: K layer distance weight interpolation
          ETAK = 0
          Do J = 1,near_points
          I = NBR(IGC(N,1),J)
          K = PGC(N,4)
          If (IDC(I,K).LT.0.9) Then
              DK = (PGC(N,1)-CXY(I,1))**2 + (PGC(N,2)-CXY(I,2))**2
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
          I = NBR(IGC(N,1),J)
          K = PGC(N,4)
          If (IDC(I,K).LT.0.9) Then
              DK = (PGC(N,1)-CXY(I,1))**2 + (PGC(N,2)-CXY(I,2))**2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      WM = (1/ETAK) * (1/DK)
                  Else
                      WM = 0
                  Endif
                  UTEMP = UTEMP + WM * U(I,K)
                  VTEMP = VTEMP + WM * V(I,K)
                  WTEMP = WTEMP + WM * W(I,K)
                  ETEMP = ETEMP + WM * ELF(I)
                  PTEMP = PTEMP + WM * PN(I,K)
              Endif
          Endif
          Enddo
          
! Step2: K-1 layer distance weight interpolation          
          ETAK = 0
          Do J = 1,near_points
          I = NBR(IGC(N,1),J)
          K = PGC(N,4) - 1
          If (IDC(I,K).LT.0.9) Then
              DK = (PGC(N,1)-CXY(I,1))**2 + (PGC(N,2)-CXY(I,2))**2
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
          UTEMP2 = 0
          VTEMP2 = 0
          WTEMP2 = 0
          PTEMP2 = 0
          ETEMP2 = 0
          
          Do J = 1,near_points
          I = NBR(IGC(N,1),J)
          K = PGC(N,4) - 1
          If (IDC(I,K).LT.0.9) Then
              DK = (PGC(N,1)-CXY(I,1))**2 + (PGC(N,2)-CXY(I,2))**2
              If (DK.LE.RITP**2) Then
                  If (DK.GT.1.0E-10) Then
                      WM = (1/ETAK) * (1/DK)
                  Else
                      WM = 0
                  Endif
                  UTEMP2 = UTEMP2 + WM * U(I,K)
                  VTEMP2 = VTEMP2 + WM * V(I,K)
                  WTEMP2 = WTEMP2 + WM * W(I,K)
                  PTEMP2 = PTEMP2 + WM * PN(I,K)
              Endif
          Endif
          Enddo
          
! Step3: K & K-1 layer linear interpolation  
          I = IGC(N,1)
          K = PGC(N,4)
          ZALPHA = Abs(ZRC(I,K-1)-PGC(N,3)) / Abs(ZRC(I,K-1)-ZRC(I,K)) 
                 
          UTEMP = ZALPHA * UTEMP + (1-ZALPHA) * UTEMP2
          VTEMP = ZALPHA * VTEMP + (1-ZALPHA) * VTEMP2
          WTEMP = ZALPHA * WTEMP + (1-ZALPHA) * WTEMP2
          
          UPT = UTEMP * PGC(N,9) - VTEMP * PGC(N,8)
          
          NUM = 0
          REST = 100
          YPL = 20
          Do while (REST.GT.1E-10) 
              YTEMP = YPL
              YPL = (YPL + 0.4*PGC(N,7)*Abs(UPT)/UMOL)/(1+log(9.8*YPL))
              REST = Abs(YPL - YTEMP) 
              NUM = NUM + 1
              If (NUM.GT.1000) REST = 0
          End do
          
          UFR = YPL * UMOL / PGC(N,7)
          UPL = 2.5 * log(9.8*0.5*PGC(N,7)*UFR / UMOL)
          
          If (UPT.GE.0) Then 
              UGT = UPL * UFR
          Else
              UGT = -1 * UPL * UFR
          Endif
          
          UGN = 0
          TKE(IGC(N,1),IGC(N,2)) = UFR**2/0.3
          TDISS(IGC(N,1),IGC(N,2)) = Sqrt((0.00008/(PGC(N,7)**2))**2+
     &            (Sqrt(TKE(IGC(N,1),IGC(N,2)))/(0.219*PGC(N,7)))**2)

c         UGT = 0.5 * (UTEMP * PGC(N,8) - VTEMP * PGC(N,7)) 
c         UGN = 0.5 * (UTEMP * PGC(N,7) + VTEMP * PGC(N,8))    

          BFU(IGC(N,1),IGC(N,2)) = UGN * PGC(N,8) + UGT * PGC(N,9)
          BFV(IGC(N,1),IGC(N,2)) = UGN * PGC(N,9) - UGT * PGC(N,8)
          BFW(IGC(N,1),IGC(N,2)) = 0.5 * WTEMP
          BFE(IGC(N,1)) = ETEMP
          BFP(IGC(N,1),IGC(N,2)) = ZALPHA * PTEMP + (1-ZALPHA) * PTEMP2
      Enddo
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      Endif
      
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
              If(IDC(I,K).EQ.1) Then
                  TKE(I,K) = 0
                  TDISS(I,K) = 100000
              Endif  
          Endif
      Enddo
      Enddo
      
      If (IBGCMODE.EQ.1) Then
          Do K = 1, KBM
          Do I = 1, IJM
              If (IDC(I,K).GT.0.9) Then
                  ELF(I) = BFE(I)
                  PN(I,K) = BFP(I,K)  
              Endif
          Enddo
          Enddo
      Endif

1000  Format (3F8.4)
1100  Format (4F10.6)
      
      
      End Subroutine