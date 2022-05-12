C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (12/05/2022)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #     
C            Direct Force Immersed Boundary Method moudle                   #                                           
C                                                                           #
C############################################################################
      Subroutine IBMALDF
      
      Include './Include/OCERM_INF'
      
      Real(kind=8) BDF(10000,KBM,3),FX(IJM,KBM),FY(IJM,KBM),FZ(IJM,KBM)
      Double precision DELTAH,DHX,DHY,DHZ,BDU,BDV,BDW
      Double precision XN,YN,ZN,RN,d,XMIN,XMAX,YMIN,YMAX
      
      
      If (NSTEP.EQ.ISTART) Then
          Write(*,*) 'DFMODE is exetuting.'
      Endif
      
      Do K = 1,KBM      
          Do I = 1,IJM
              FX(I,K) = 0
              FY(I,K) = 0
              FZ(I,K) = 0
          Enddo
      Enddo 
      
C===========================================================================C      
C         Delta-Force source term caculation
C===========================================================================C       
  
      Do 200 K = 1,KBM
C---------------------------------------------------------------------------C
      Do J = 1,IBN(1),IBSTEP  
      BDU = 0
      BDV = 0
      XMIN = BXY(1,J,1) - 1.5*d
      XMAX = BXY(1,J,1) + 1.5*d
      YMIN = BXY(1,J,2) - 1.5*d
      YMAX = BXY(1,J,2) + 1.5*d
      
      If(NSTEP.EQ.ISTART) Then
          Write(IBMDIB,1300) BXY(1,J,1),BXY(1,J,2),BXY(1,J,3)
      Endif
           
! Caculate BDU,BDV,BDW 
      Do I = 1,IJM
          DELTAH = 0
          DHX = 0
          DHY = 0
          DHZ = 0
          
          If ((CXY(I,1).GT.XMIN).AND.(CXY(I,1).LT.XMAX).AND. 
     &        (CXY(I,2).GT.YMIN).AND.(CXY(I,2).LT.YMAX)) Then
              XN = Abs((CXY(I,1) - BXY(1,J,1))/d)
              YN = Abs((CXY(I,2) - BXY(1,J,2))/d)
              ZN = Abs((ZZ(K) - BXY(1,J,3))/((d-HP(1))/HP(1))) 
              RN = Sqrt(XN**2 + YN**2)
              
              If (XN.LE.0.5) Then
                  DHX = (2 + 2*Sqrt(1-3*XN**2))/(6*d)
              Else if ((XN.GT.0.5).AND.(XN.LE.1.5)) Then
                  DHX = (5 - 3*XN - Sqrt(1-3*(1-XN)**2))/(6*d)
              Else
                  DHX = 0
              Endif
              
              If (YN.LE.0.5) Then
                  DHY = (2 + 2*Sqrt(1-3*YN**2))/(6*d)
              Else if ((YN.GT.0.5).AND.(YN.LE.1.5)) Then
                  DHY = (5 - 3*YN - Sqrt(1-3*(1-YN)**2))/(6*d)
              Else
                  DHY = 0
              Endif
              
              If (IBMDOMAIN.EQ.2) Then
                  DHZ = 1
              Else if (IBMDOMAIN.EQ.3) Then
                  If (ZN.LE.0.5) Then
                      DHZ = (2 + 2*Sqrt(1-3*ZN**2))/(6*d)
                  Else if ((ZN.GT.0.5).AND.(ZN.LE.1.5)) Then
                      DHZ = (5 - 3*ZN - Sqrt(1-3*(1-ZN)**2))/(6*d)
                  Else
                      DHZ = 0
                  Endif  
              Endif
              
              DELTAH = DHX * DHY * DHZ
              
              If (NSTEP.EQ.ISTART) Then
              If (RN.LE.1.5)
     &            Write(IBMDF,1300) CXY(I,1),CXY(I,2),ZZ(K)
              Endif
          Endif 
          BDU = BDU + (d*d) * U(I,K) * DELTAH
          BDV = BDV + (d*d) * V(I,K) * DELTAH
          BDW = BDW + (d*d) * W(I,K) * DELTAH
      Enddo
      
! Caculate BDF      
      BDF(J,K,1) = (0 - BDU) / DTI
      BDF(J,K,2) = (0 - BDV) / DTI
      BDF(J,K,3) = (0 - BDW) / DTI
      Enddo
C---------------------------------------------------------------------------C
200   Continue

      Do K = 1,KBM
      Do I = 1,IJM
          Do J = 1,IBN(1),IBSTEP
              DELTAH = 0
              DHX = 0
              DHY = 0
              DHZ = 0
              
              XN = Abs((CXY(I,1) - BXY(1,J,1))/d)
              YN = Abs((CXY(I,2) - BXY(1,J,2))/d)
              ZN = Abs((ZZ(K) - BXY(1,J,3))/d)
              RN = Sqrt(XN**2 + YN**2)
              
              If (XN.LE.0.5) Then
                  DHX = (2 + 2*Sqrt(1-3*XN**2))/(6*d)
              Else if ((XN.GT.0.5).AND.(XN.LE.1.5)) Then
                  DHX = (5 - 3*XN - Sqrt(1-3*(1-XN)**2))/(6*d)
              Else
                  DHX = 0
              Endif
              
              If (YN.LE.0.5) Then
                  DHY = (2 + 2*Sqrt(1-3*YN**2))/(6*d)
              Else if ((YN.GT.0.5).AND.(YN.LE.1.5)) Then
                  DHY = (5 - 3*YN - Sqrt(1-3*(1-YN)**2))/(6*d)
              Else
                  DHY = 0
              Endif
              
              If (IBMDOMAIN.EQ.2) Then
                  DHZ = 1
              Else if (IBMDOMAIN.EQ.3) Then
                  If (ZN.LE.0.5) Then
                      DHZ = (2 + 2*Sqrt(1-3*ZN**2))/(6*d)
                  Else if ((ZN.GT.0.5).AND.(ZN.LE.1.5)) Then
                      DHZ = (5 - 3*ZN - Sqrt(1-3*(1-ZN)**2))/(6*d)
                  Else
                      DHZ = 0
                  Endif  
              Endif
              
              DELTAH = DHX * DHY * DHZ
              
              FX(I,K) = FX(I,K) + RIFL * BDF(J,K,1) * DELTAH
              FY(I,K) = FY(I,K) + RIFL * BDF(J,K,2) * DELTAH
              FZ(I,K) = FZ(I,K) + RIFL * BDF(J,K,2) * DELTAH
          Enddo
      Enddo
      Enddo
      
C===========================================================================C         
C         Caculate final source term BFX,BFY,BFZ
C===========================================================================C 
      
      Do K = 1,KBM
          Do I = 1,IJM
              BFX(I,K) = FX(I,K)
              BFY(I,K) = FY(I,K)
              BFZ(I,K) = FZ(I,K)
          Enddo
      Enddo

1000  Format (100A1)
1100  Format (/1X,80A1/)
1200  Format (F8.4,I8,F10.7) 
1300  Format (3F8.4) 
      
      
      End Subroutine