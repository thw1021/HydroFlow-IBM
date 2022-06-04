C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (12/05/2022)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #
C        Detecte differnet kind of cells and interpolation points.          #                                           
C                                                                           #
C############################################################################
      Subroutine IBMALIDC
      
      Include './Include/OCERM_INF'
      
      Real(kind=8) BXY2(10000,KB,3),BXY3(10000,3,3)
      Real(kind=8) PLANE(3,3),LINE(2,2),CRP(3)
      Double precision XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN,FL,FR,GK,ZIB,Z2D
      Double precision LA,LB,LC,LDIS,LDMIN,XIM,YIM,ZIM,ZNE
      Double precision PA,PB,PC,PD,PDIS,PDMIN,PNDMIN,KDMIN,TLINE
      Double precision IPX,IPY,IPZ,IPK,IPDIS,IPCOS,IPSIN
      Integer nn,nc,NNB,NTC,S1,N_C,N_SL,N_GH,NBMIN,NPMIN,KMIN
      Logical is_coincide
      

C===========================================================================C      
C         Initialize the arrays      
C===========================================================================C

      NUMSL = 0
      NUMGH = 0
      Z2D = -1
      
      Do K = 1,KBM
      Do I = 1,IJM
          IDC(I,K) = 0
          ZRC(I,K) = ZZ(K)*(HC(I)+ELF(I))+ELF(I)
      Enddo
      Enddo 
      
      Do N = 1,IJM
          IGC(N,1) = 0
          IGC(N,2) = 0
          Do J = 1,9
              PGC(N,J) = 0
          Enddo
      Enddo  
      
      Do N = 1,IJM*2
          ISC(N,1) = 0
          ISC(N,2) = 0
      Enddo
      
      Do I = 1,3
          CRP(I) = 0
      Enddo
      
C===========================================================================C      
C         2D Immersed boundary detection      
C===========================================================================C   
      
      If (IBMDOMAIN.EQ.2) Then 
          
      Do 200 nn = 1,NUMIB
        
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      Do I = 1,10000
      Do K = 1,KBM
          BXY2(I,K,1) = 0
          BXY2(I,K,2) = 0
          BXY2(I,K,3) = 0
      Enddo
      Enddo
          
      Do K = 1,KBM
      Do J = 1,IBN(nn)
          BXY(nn,J,3) = ZZ(K)
          BXY2(J,K,1) = BXY(nn,J,1)
          BXY2(J,K,2) = BXY(nn,J,2)
          BXY2(J,K,3) = BXY(nn,J,3)
      Enddo
      Enddo
      
! Get XMAX,XMIN,YMAX,YMIN 
      XMAX = BXY(nn,1,1)
      XMIN = BXY(nn,1,1)
      YMAX = BXY(nn,1,2)
      YMIN = BXY(nn,1,2)
     
      Do I = 1,IBN(nn)
          If (BXY(nn,I,1).GE.XMAX) Then
              XMAX = BXY(nn,I,1)
          Endif
          If (BXY(nn,I,1).LE.XMIN) Then
              XMIN = BXY(nn,I,1)
          Endif
          If (BXY(nn,I,2).GE.YMAX) Then
              YMAX = BXY(nn,I,2)
          Endif
          If (BXY(nn,I,2).LE.YMIN) Then
              YMIN = BXY(nn,I,2)
          Endif
      Enddo 

! Solid cells detection 
      Do 201 K = 1,KBM    
      Do I = 1,IJM    
          
      If (ZRC(I,K).LE.Z2D) Then
          
C***************************************************************************C
C             Ray-casting Method Begins                                     C
C***************************************************************************C            
      If ((CXY(I,1).LT.XMAX).AND.(CXY(I,1).GT.XMIN).AND. 
     &    (CXY(I,2).LT.YMAX).AND.(CXY(I,2).GT.YMIN)) Then
      S1 = 0   
      Do J = 1,IBN(nn) 
          If (J.LT.IBN(nn)) Then
              CRP(1)=(CXY(I,1)-BXY2(J,K,1))*(BXY2(J+1,K,2)-CXY(I,2))
     &              -(BXY2(J+1,K,1)-CXY(I,1))*(CXY(I,2)-BXY2(J,K,2))
                           
              If (((BXY2(J,K,1).LE.CXY(I,1)).AND.
     &            (BXY2(J+1,K,1).GT.CXY(I,1)).AND.(CRP(1).GE.0))
     &            .OR.
     &            ((BXY2(J,K,1).GE.CXY(I,1)).AND.
     &            (BXY2(J+1,K,1).LT.CXY(I,1)).AND.(CRP(1).LE.0))) 
     &        Then               
                  S1 = S1 + 1
                  
                  If (CRP(1) .EQ. 0) Then       ! Case 1
                      S1 = 1
                      Exit
                  Endif
                  
                  If(J.EQ.1)Then       ! Case 2
                      If(((BXY2(J,K,1).EQ.CXY(I,1)).AND.
     &                    (BXY2(2,K,1).GE.CXY(I,1)).AND.
     &                    (BXY2(IBN(nn),K,1).GE.CXY(I,1))) 
     &                    .OR. 
     &                    ((BXY2(J,K,1).EQ.CXY(I,1)).AND.
     &                    (BXY2(2,K,1).LE.CXY(I,1)).AND.
     &                    (BXY2(IBN(nn),K,1).LE.CXY(I,1)))) 
     &                Then
                          S1 = S1 - 1
                      Endif
                  Else 
                      If(((BXY2(J,K,1).EQ.CXY(I,1)).AND.
     &                    (BXY2(J+1,K,1).GT.CXY(I,1)).AND.
     &                    (BXY2(J-1,K,1).GT.CXY(I,1))) 
     &                        .OR. 
     &                    ((BXY2(J,K,1).EQ.CXY(I,1)).AND.
     &                    (BXY2(J+1,K,1).LT.CXY(I,1)).AND.
     &                    (BXY2(J-1,K,1).LT.CXY(I,1)))) 
     &                Then
                          S1 = S1 - 1
                  Endif
                  Endif
                  
              Endif

              If ((BXY2(J,K,1).EQ.CXY(I,1)).AND.       
     &           (BXY2(J+1,K,1).EQ.CXY(I,1))) Then      ! Case 3
                  If (((CXY(I,2).LE.BXY2(J,K,2)).AND.
     &                (CXY(I,2).GE.BXY2(J+1,K,2)))
     &                .OR.
     &                ((CXY(I,2).GE.BXY2(J,K,2)).AND.
     &                (CXY(I,2).LE.BXY2(J+1,K,2))))
     &            Then
                      S1 = 1
                      Exit
                  Endif
              Endif
          
          Else if(J.EQ.IBN(nn)) Then
              CRP(1) = (CXY(I,1)-BXY2(J,K,1))*(BXY2(1,K,2)-CXY(I,2))
     &                -(BXY2(1,K,1)-CXY(I,1))*(CXY(I,2)-BXY2(J,K,2))
              
              If (((BXY2(J,K,1).LE.CXY(I,1)).AND.
     &            (BXY2(1,K,1).GT.CXY(I,1)).AND.(CRP(1).GE.0))
     &            .OR.
     &            ((BXY2(J,K,1).GE.CXY(I,1)).AND.
     &            (BXY2(1,K,1).LT.CXY(I,1)).AND.(CRP(1).LE.0))) 
     &        Then                             
                  S1 = S1 + 1
                  
                  If (CRP(1) .EQ. 0) Then       ! Case 1
                      S1 = 1
                      Exit
                  Endif
                  
                  If(((BXY2(J,K,1).EQ.CXY(I,1)).AND.       ! Case 2
     &               (BXY2(1,K,1).GT.CXY(I,1)).AND.
     &               (BXY2(J-1,K,1).GT.CXY(I,1))) 
     &                .OR. 
     &               ((BXY2(J,K,1).EQ.CXY(I,1)).AND.
     &               (BXY2(1,K,1).LT.CXY(I,1)).AND.
     &               (BXY2(J-1,K,1).LT.CXY(I,1)))) 
     &            Then
                     S1 = S1 - 1
                  Endif

              Endif
                  
              If ((BXY2(J,K,1).EQ.CXY(I,1)).AND.       
     &           (BXY2(1,K,1).EQ.CXY(I,1))) Then      ! Case 3
                  If (((CXY(I,2).LE.BXY2(J,K,2)).AND.
     &                (CXY(I,2).GE.BXY2(1,K,2)))
     &                .OR.
     &                ((CXY(I,2).GE.BXY2(J,K,2)).AND.
     &                (CXY(I,2).LE.BXY2(1,K,2))))
     &            Then
                      S1 = 1
                      Exit
                  Endif
              Endif
      
          Endif
      Enddo
    
      If (mod(S1,2) .GT. 0) Then
          IDC(I,K) = 1     
      Endif
      Endif     
C***************************************************************************C
C             End of Ray-casting method                                     C
C***************************************************************************C 
      
      Endif
      Enddo 
    
201   Continue  
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      
200   Continue
                               
      Endif
      
C===========================================================================c     
C         3D Immersed boundary detection     
C===========================================================================c

      If (IBMDOMAIN.EQ.3) Then
          
      Do 300 nn = 1,NUMIB
          
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C  
      Do I = 1,10000
      Do K = 1,3
          BXY3(I,K,1) = 0
          BXY3(I,K,1) = 0
          BXY3(I,K,1) = 0
      Enddo
      Enddo
      
      NTC = 0     ! NTC is used to count triangle element
      Do J = 1,IBN(nn)/3
          NTC = NTC+1
              
          BXY3(NTC,1,1) = BXY(nn,3*J-2,1)
          BXY3(NTC,1,2) = BXY(nn,3*J-2,2)
          BXY3(NTC,1,3) = BXY(nn,3*J-2,3)

          BXY3(NTC,2,1) = BXY(nn,3*J-1,1)
          BXY3(NTC,2,2) = BXY(nn,3*J-1,2)
          BXY3(NTC,2,3) = BXY(nn,3*J-1,3)

          BXY3(NTC,3,1) = BXY(nn,3*J,1)
          BXY3(NTC,3,2) = BXY(nn,3*J,2)
          BXY3(NTC,3,3) = BXY(nn,3*J,3)
      Enddo

! Get XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN
      XMAX = BXY(1,1,1)
      XMIN = BXY(1,1,1)
      YMAX = BXY(1,1,2)
      YMIN = BXY(1,1,2)
      ZMAX = BXY(1,1,3)
      ZMIN = BXY(1,1,3)

      Do I = 1,IBN(nn)
          If (BXY(nn,I,1).GE.XMAX) Then
              XMAX = BXY(nn,I,1)
          Endif
          If (BXY(nn,I,1).LE.XMIN) Then
              XMIN = BXY(nn,I,1)
          Endif
          If (BXY(nn,I,2).GE.YMAX) Then
              YMAX = BXY(nn,I,2)
          Endif
          If (BXY(nn,I,2).LE.YMIN) Then
              YMIN = BXY(nn,I,2)
          Endif
          If (BXY(nn,I,3).GE.ZMAX) Then
              ZMAX = BXY(nn,I,3)
          Endif
          If (BXY(nn,I,3).LE.ZMIN) Then
              ZMIN = BXY(nn,I,3)
          Endif
      Enddo 
      
! Solid cell detection 
      Do 301 K = 1,KBM   
      Do I = 1,IJM   
      
C***************************************************************************C
C             Ray-casting Method Begins                                     C
C***************************************************************************C            
      If ((CXY(I,1).LE.XMAX).AND.(CXY(I,1).GE.XMIN).AND. 
     &    (CXY(I,2).LE.YMAX).AND.(CXY(I,2).GE.YMIN).AND.
     &    (ZRC(I,K).LE.ZMAX).AND.(ZRC(I,K).GE.ZMIN)) Then 
      
      S1 = 0
      Do 302 nc = 1,NTC
      GK = (BXY3(nc,2,1)-BXY3(nc,1,1))*(BXY3(nc,3,2)-BXY3(nc,1,2))
     &    -(BXY3(nc,2,2)-BXY3(nc,1,2))*(BXY3(nc,3,1)-BXY3(nc,1,1))
      If (GK.NE.0) Then
          is_coincide = (((CXY(I,1).EQ.BXY3(nc,1,1)).AND.
     &                (CXY(I,2).EQ.BXY3(nc,1,2))).OR.
     &                ((CXY(I,1).EQ.BXY3(nc,2,1)).AND.
     &                (CXY(I,2).EQ.BXY3(nc,2,2))).OR.
     &                ((CXY(I,1).EQ.BXY3(nc,3,1)).AND.
     &                (CXY(I,2).EQ.BXY3(nc,3,2))))
          If (is_coincide) Exit

          CRP(1) = (BXY3(nc,1,1)-CXY(I,1))*(BXY3(nc,2,2)-CXY(I,2))
     &            -(BXY3(nc,1,2)-CXY(I,2))*(BXY3(nc,2,1)-CXY(I,1))
          CRP(2) = (BXY3(nc,2,1)-CXY(I,1))*(BXY3(nc,3,2)-CXY(I,2))
     &            -(BXY3(nc,2,2)-CXY(I,2))*(BXY3(nc,3,1)-CXY(I,1))
          CRP(3) = (BXY3(nc,3,1)-CXY(I,1))*(BXY3(nc,1,2)-CXY(I,2))
     &            -(BXY3(nc,3,2)-CXY(I,2))*(BXY3(nc,1,1)-CXY(I,1))

          If (((CRP(1).GE.0).AND.(CRP(2).GE.0).AND.(CRP(3).GE.0)).OR.
     &        ((CRP(1).LE.0).AND.(CRP(2).LE.0).AND.(CRP(3).LE.0))) Then
          If ((CRP(1).EQ.0).OR.(CRP(2).EQ.0).OR.(CRP(3).EQ.0)) Exit

          FL=((BXY3(nc,2,1)-BXY3(nc,1,1))*(BXY3(nc,3,2)-BXY3(nc,1,2))-
     &        (BXY3(nc,2,2)-BXY3(nc,1,2))*(BXY3(nc,3,1)-BXY3(nc,1,1))) 
          FR =(((CXY(I,2)-BXY3(nc,1,2))*(BXY3(nc,2,1)-BXY3(nc,1,1))*
     &        (BXY3(nc,3,3)-BXY3(nc,1,3)))+
     &        ((CXY(I,1)-BXY3(nc,1,1))*(BXY3(nc,2,3)-BXY3(nc,1,3))*
     &        (BXY3(nc,3,2)-BXY3(nc,1,2)))-
     &        ((CXY(I,1)-BXY3(nc,1,1))*(BXY3(nc,2,2)-BXY3(nc,1,2))*
     &        (BXY3(nc,3,3)-BXY3(nc,1,3)))-
     &        ((CXY(I,2)-BXY3(nc,1,2))*(BXY3(nc,2,3)-BXY3(nc,1,3))*
     &        (BXY3(nc,3,1)-BXY3(nc,1,1))))
          ZIB = FR/FL+BXY3(nc,1,3)
          
          If (ZRC(I,K).LT.ZIB) Then
              S1 = S1 + 1
          Else If (ZRC(I,K).EQ.ZIB) Then
              S1 = 1
              Exit
          Endif
          Endif
      Endif 
                   
302   Continue
      
      If (mod(S1,2) .GT. 0) Then
          IDC(I,K) = 1
      Endif
      Endif
C***************************************************************************C
C             End of Ray-casting method                                     C
C***************************************************************************C
         
      Enddo
      
301   Continue  
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C  

300   Continue

      Endif     
          
C===========================================================================C
C         Search ghost cells and image points
C===========================================================================C 
      
      Do K = 1,KBM
      Do I = 1,IJM
          If (IDC(I,K).EQ.1) Then  
              Do J = 1,CELL_POLYGEN(I)
              If ((CELL_SIDE(I,J,2).GT.0).AND.
     &           (IDC(CELL_SIDE(I,J,2),K).LT.0.9)) Then
                  IDC(CELL_SIDE(I,J,2),K) = 2
              Endif
              If (K.EQ.1) Then
                  If (IDC(I,2).LT.0.9) IDC(I,2) = 3
              Else if (K.EQ.KBM) Then
                  If (IDC(I,KBM-1).LT.0.9) IDC(I,KBM-1) = 4
              Else
                  If (IDC(I,K-1).LT.0.9) IDC(I,K-1) = 4
                  If (IDC(I,K+1).LT.0.9) IDC(I,K+1) = 3
              Endif
              Enddo
          
              Do J = 1,CELL_POLYGEN(I)
              Do N = 1,POINT_CELL(CELL_SIDE(I,J,3))
              If ((INDEX_POINT(CELL_SIDE(I,J,3),N).GT.0).AND.
     &           (IDC(INDEX_POINT(CELL_SIDE(I,J,3),N),K).LT.0.9)) Then
                  IDC(INDEX_POINT(CELL_SIDE(I,J,3),N),K) = 2
              Endif
              Enddo
              Enddo
          Endif
      Enddo
      Enddo

C===========================================================================C
C         IBM variable assignment
C===========================================================================C 

      N_SL = 0
      N_GH = 0
      NNB = 30
     
      Do K = 1,KBM
      Do I = 1,IJM
          If (IDC(I,K).EQ.1) Then
              N_SL = N_SL + 1
              ISC(N_SL,1) = I
              ISC(N_SL,2) = K
          Endif
          If (IDC(I,K).GT.1.9) Then
              N_GH = N_GH + 1
              IGC(N_GH,1) = I
              IGC(N_GH,2) = K 
          Endif
      Enddo
      Enddo
      
      NUMSL = N_SL
      NUMGH = N_GH
      
      ! search image point
      If (IBMDOMAIN.EQ.2) Then
      Do N = 1,NUMGH
          
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
      LDMIN = 10000
      NPMIN = 1
      NBMIN = 1
      
      Do I = 1,2
      Do J = 1,2
          LINE(I,J) = 1
      Enddo
      Enddo
	               
      Do nn = 1,NUMIB
      Do J = 1,IBN(nn)
          If (J.LT.IBN(nn)) Then
          LDIS =(CXY(IGC(N,1),1)-0.5*(BXY(nn,J,1)+BXY(nn,J+1,1)))**2 
     &         +(CXY(IGC(N,1),2)-0.5*(BXY(nn,J,2)+BXY(nn,J+1,2)))**2
          Else
          LDIS =(CXY(IGC(N,1),1)-0.5*(BXY(nn,IBN(nn),1)+BXY(nn,1,1)))**2
     &         +(CXY(IGC(N,1),2)-0.5*(BXY(nn,IBN(nn),2)+BXY(nn,1,2)))**2
          Endif
		  
          If (LDIS.LE.LDMIN) Then 
              LDMIN = LDIS
              NPMIN = J
              NBMIN = nn
              If (J.LT.IBN(nn)) Then
			    LINE(1,1) = BXY(nn,J,1)
			    LINE(1,2) = BXY(nn,J,2)
			    LINE(2,1) = BXY(nn,J+1,1)
			    LINE(2,2) = BXY(nn,J+1,2)
			Else
			    LINE(1,1) = BXY(nn,IBN(nn),1)
			    LINE(1,2) = BXY(nn,IBN(nn),2)
			    LINE(2,1) = BXY(nn,1,1)
			    LINE(2,2) = BXY(nn,1,2)
			Endif
          Endif
      Enddo
      Enddo
      
      If (IDC(IGC(N,1),IGC(N,2)).EQ.2) Then
          LA = LINE(2,2)-LINE(1,2)
          LB = LINE(1,1)-LINE(2,1)
          LC = LINE(1,2)*LINE(2,1)-LINE(2,2)*LINE(1,1)
          XIM = (-2*LA*LC+(LB**2-LA**2)*CXY(IGC(N,1),1)
     &         -2*LA*LB*CXY(IGC(N,1),2))/(LA**2+LB**2)
          YIM = (-2*LB*LC+(LA**2-LB**2)*CXY(IGC(N,1),2)
     &         -2*LA*LB*CXY(IGC(N,1),1))/(LA**2+LB**2)
          ZIM = ZRC(IGC(N,1),IGC(N,2))
          
          IPX = 1.5*CXY(IGC(N,1),1) - 0.5*XIM
          IPY = 1.5*CXY(IGC(N,1),2) - 0.5*YIM
          IPZ = ZIM
          IPK = IGC(N,2)
          IPDIS = Abs(LA*IPX+LB*IPY+LC) / Sqrt(LA**2+LB**2)
          IPCOS = LA / Sqrt(LA**2+LB**2)
          IPSIN = LB / Sqrt(LA**2+LB**2)
      Else if (IDC(IGC(N,1),IGC(N,2)).GE.3) Then
          XIM = CXY(IGC(N,1),1)
          YIM = CXY(IGC(N,1),2)
          IPX = XIM
          IPY = YIM
          If (IDC(IGC(N,1),IGC(N,2)).EQ.3) Then
              If (IGC(N,2).GT.1) IPK = IGC(N,2) + 1
          Else if (IDC(IGC(N,1),IGC(N,2)).EQ.4) Then
              If (IGC(N,2).LT.KBM) IPK = IGC(N,2) - 1
          Endif  
          IPZ = ZRC(IGC(N,1),IPK)
          IPDIS = Abs(IPZ - Z2D)
          IPCOS = 0
          IPSIN = 1
      Endif
     
      PGC(N,1) = IPX
      PGC(N,2) = IPY
      PGC(N,3) = IPZ
      PGC(N,4) = IPK
      PGC(N,5) = NBMIN
      PGC(N,6) = NPMIN
      PGC(N,7) = IPDIS
      PGC(N,8) = IPCOS
      PGC(N,9) = IPSIN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C 
          
      Enddo
      Endif
      
      If (IBMDOMAIN.EQ.3) Then
      Do N = 1,NUMGH
          
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      PDMIN = 10000
      KDMIN = 10000
      
      Do I = 1,3
      Do J = 1,3
          PLANE(I,J) = 1
      Enddo
      Enddo
          
      Do nn = 1,NUMIB    
          Do I = 1,10000
          Do K = 1,3
              BXY3(I,K,1) = 0
              BXY3(I,K,1) = 0
              BXY3(I,K,1) = 0
          Enddo
          Enddo
      
          NTC = 0     ! NTC is used to count triangle element
          Do J = 1,IBN(nn)/3
              NTC = NTC+1
              
              BXY3(NTC,1,1) = BXY(nn,3*J-2,1)
              BXY3(NTC,1,2) = BXY(nn,3*J-2,2)
              BXY3(NTC,1,3) = BXY(nn,3*J-2,3)

              BXY3(NTC,2,1) = BXY(nn,3*J-1,1)
              BXY3(NTC,2,2) = BXY(nn,3*J-1,2)
              BXY3(NTC,2,3) = BXY(nn,3*J-1,3)

              BXY3(NTC,3,1) = BXY(nn,3*J,1)
              BXY3(NTC,3,2) = BXY(nn,3*J,2)
              BXY3(NTC,3,3) = BXY(nn,3*J,3)
          Enddo
          
	    Do m = 1,NTC
		PDIS = Sqrt((CXY(IGC(N,1),1)-((BXY3(m,1,1))
     &                +(BXY3(m,2,1))+(BXY3(m,3,1)))/3)**2+
     &                 (CXY(IGC(N,1),2)-((BXY3(m,1,2))
     &                +(BXY3(m,2,2))+(BXY3(m,3,2)))/3)**2+
     &                 (ZRC(IGC(N,1),IGC(N,2))-((BXY3(m,1,3))
     &                +(BXY3(m,2,3))+(BXY3(m,3,3)))/3)**2)
     
		If (PDIS.LE.PDMIN) Then
			PDMIN = PDIS
			NPMIN = m
              NBMIN = nn
              PLANE(1,1) = BXY3(m,1,1)
              PLANE(1,2) = BXY3(m,1,2)
              PLANE(1,3) = BXY3(m,1,3)
              PLANE(2,1) = BXY3(m,2,1)
              PLANE(2,2) = BXY3(m,2,2)
              PLANE(2,3) = BXY3(m,2,3)
              PLANE(3,1) = BXY3(m,3,1)
              PLANE(3,2) = BXY3(m,3,2)
              PLANE(3,3) = BXY3(m,3,3)
		Endif
          Enddo 
      Enddo
      
      PA = (PLANE(2,2)-PLANE(1,2)) * (PLANE(3,3)-PLANE(1,3))
     &    -(PLANE(2,3)-PLANE(1,3)) * (PLANE(3,2)-PLANE(1,2))
      PB = (PLANE(2,3)-PLANE(1,3)) * (PLANE(3,1)-PLANE(1,1))
     &    -(PLANE(2,1)-PLANE(1,1)) * (PLANE(3,3)-PLANE(1,3))
      PC = (PLANE(2,1)-PLANE(1,1)) * (PLANE(3,2)-PLANE(1,2))
     &    -(PLANE(2,2)-PLANE(1,2)) * (PLANE(3,1)-PLANE(1,1))
      PD = -1 * (PA*PLANE(1,1)+PB*PLANE(1,2)+PC*PLANE(1,3))
      
      TLINE = -1 * (PA*CXY(IGC(N,1),1) + PB*CXY(IGC(N,1),2)+
     &         PC*ZRC(IGC(N,1),IGC(N,2)) + PD) /
     &        (PA**2 + PB**2 + PC**2)
     
      XIM = 2*PA*TLINE + CXY(IGC(N,1),1)
      YIM = 2*PB*TLINE + CXY(IGC(N,1),2)
      ZIM = 2*PC*TLINE + ZRC(IGC(N,1),IGC(N,2))
      
      IPX = 1.5*CXY(IGC(N,1),1) - 0.5*XIM
      IPY = 1.5*CXY(IGC(N,1),2) - 0.5*YIM
      IPZ = 1.5*ZRC(IGC(N,1),IGC(N,2)) - 0.5*ZIM
      IPDIS = Abs(PA*XIM+PB*YIM+PC*ZIM+PD) / Sqrt(PA**2+PB**2+PC**2)
      IPCOS = PA / Sqrt(PA**2+PB**2+PC**2)
      IPSIN = PB / Sqrt(PA**2+PB**2+PC**2)
      
      ZNE = 1
      Do While (ZRC(IGC(N,1),ZNE).GE.IPZ)
          ZNE = ZNE + 1
      Enddo
      IPK = ZNE

      PGC(N,1) = IPX
      PGC(N,2) = IPY
      PGC(N,3) = IPZ
      PGC(N,4) = IPK
      PGC(N,5) = NBMIN
      PGC(N,6) = NPMIN
      PGC(N,7) = IPDIS
      PGC(N,8) = IPCOS
      PGC(N,9) = IPSIN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
      
      Enddo
      Endif
          
      End Subroutine
