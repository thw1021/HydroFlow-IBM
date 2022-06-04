C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (12/05/2022)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #
C                Detecte differnet kind of points.                          #                                           
C                                                                           #
C############################################################################      
      Subroutine IBMALIDP
      
      Include './Include/OCERM_INF'
      
      Real(kind=8) BXY2(10000,KB,3),BXY3(10000,3,3),CRP(3)
      Double precision XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN,FL,FR,ZIB,Z2D,GK
      Integer nn,nc,NTC,S1,S2,S3
      

C===========================================================================C      
C         Initialize the arrays      
C===========================================================================C

      Z2D = 100
      
      Do K = 1,KB
      Do I = 1,IJP
          IDP(I,K) = 0
          ZRP(I,K) = Z(K)*(HP(I)+ELFV(I))+ELFV(I)
      Enddo
      Enddo

      Do I = 1,3
          CRP(I) = 0
      Enddo
      
C===========================================================================C      
C      2D Immersed boundary detection     
C===========================================================================C   
      
      If (IBMDOMAIN.EQ.2) Then
         
      Do 200 nn = 1,NUMIB
          
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C   
      Do I = 1,10000
      Do K = 1,KB
          BXY2(I,K,1) = 0
          BXY2(I,K,2) = 0
          BXY2(I,K,3) = 0
      Enddo
      Enddo
      
      Do K = 1,KB
      Do J = 1,IBN(nn)
          BXY(nn,J,3) = Z(K)
          BXY2(J,K,1) = BXY(nn,J,1)
          BXY2(J,K,2) = BXY(nn,J,2)
          BXY2(J,K,3) = BXY(nn,J,3)
      Enddo
      Enddo   
      
! Get XMAX,XMIN,YMAX,YMIN 
      XMAX = BXY(1,1,1)
      XMIN = BXY(1,1,1)
      YMAX = BXY(1,1,2)
      YMIN = BXY(1,1,2)

      Do i = 1,IBN(nn)
          If (BXY(nn,i,1).GE.XMAX) Then
              XMAX = BXY(nn,i,1)
          Endif
          If (BXY(nn,i,1).LE.XMIN) Then
              XMIN = BXY(nn,i,1)
          Endif
          If (BXY(nn,i,2).GE.YMAX) Then
              YMAX = BXY(nn,i,2)
          Endif
          If (BXY(nn,i,2).LE.YMIN) Then
              YMIN = BXY(nn,i,2)
          Endif
      Enddo 
      
! Solid cells detection 
      Do 201 K = 1,KB     
      Do I = 1,IJP        
          
      If (ZRP(I,K).LE.Z2D) Then
          
C***************************************************************************C
C             Ray-casting Method Begins                                     C
C***************************************************************************C            
      If ((PXY(I,1).LT.XMAX).AND.(PXY(I,1).GT.XMIN).AND. 
     &    (PXY(I,2).LT.YMAX).AND.(PXY(I,2).GT.YMIN)) Then
     &
      S1 = 0   
      Do J = 1,IBN(nn) 
          If (J.LT.IBN(nn)) Then
              CRP(1)=(PXY(I,1)-BXY2(J,K,1))*(BXY2(J+1,K,2)-PXY(I,2))
     &              -(BXY2(J+1,K,1)-PXY(I,1))*(PXY(I,2)-BXY2(J,K,2))
                           
              If (((BXY2(J,K,1).LE.PXY(I,1)).AND.
     &            (BXY2(J+1,K,1).GT.PXY(I,1)).AND.(CRP(1).GE.0))
     &            .OR.
     &            ((BXY2(J,K,1).GE.PXY(I,1)).AND.
     &            (BXY2(J+1,K,1).LT.PXY(I,1)).AND.(CRP(1).LE.0))) 
     &        Then               
                  S1 = S1 + 1
                  
                  If (CRP(1) .EQ. 0) Then       ! Case 1
                      S1 = 1
                      Exit
                  Endif
                  
                  If(J.EQ.1)Then       ! Case 2
                      If(((BXY2(J,K,1).EQ.PXY(I,1)).AND.
     &                    (BXY2(2,K,1).GE.PXY(I,1)).AND.
     &                    (BXY2(IBN(nn),K,1).GE.PXY(I,1))) 
     &                    .OR. 
     &                    ((BXY2(J,K,1).EQ.PXY(I,1)).AND.
     &                    (BXY2(2,K,1).LE.PXY(I,1)).AND.
     &                    (BXY2(IBN(nn),K,1).LE.PXY(I,1)))) 
     &                Then
                          S1 = S1 - 1
                      Endif
                  Else 
                      If(((BXY2(J,K,1).EQ.PXY(I,1)).AND.
     &                    (BXY2(J+1,K,1).GT.PXY(I,1)).AND.
     &                    (BXY2(J-1,K,1).GT.PXY(I,1))) 
     &                        .OR. 
     &                    ((BXY2(J,K,1).EQ.PXY(I,1)).AND.
     &                    (BXY2(J+1,K,1).LT.PXY(I,1)).AND.
     &                    (BXY2(J-1,K,1).LT.PXY(I,1)))) 
     &                Then
                          S1 = S1 - 1
                  Endif
                  Endif
                  
              Endif

              If ((BXY2(J,K,1).EQ.PXY(I,1)).AND.       
     &           (BXY2(J+1,K,1).EQ.PXY(I,1))) Then      ! Case 3
                  If (((PXY(I,2).LE.BXY2(J,K,2)).AND.
     &                (PXY(I,2).GE.BXY2(J+1,K,2)))
     &                .OR.
     &                ((PXY(I,2).GE.BXY2(J,K,2)).AND.
     &                (PXY(I,2).LE.BXY2(J+1,K,2))))
     &            Then
                      S1 = 1
                      Exit
                  Endif
              Endif
          
          Else if(J.EQ.IBN(nn)) Then
              CRP(1) = (PXY(I,1)-BXY2(J,K,1))*(BXY2(1,K,2)-PXY(I,2))
     &                -(BXY2(1,K,1)-PXY(I,1))*(PXY(I,2)-BXY2(J,K,2))
              
              If (((BXY2(J,K,1).LE.PXY(I,1)).AND.
     &            (BXY2(1,K,1).GT.PXY(I,1)).AND.(CRP(1).GE.0))
     &            .OR.
     &            ((BXY2(J,K,1).GE.PXY(I,1)).AND.
     &            (BXY2(1,K,1).LT.PXY(I,1)).AND.(CRP(1).LE.0))) 
     &        Then                             
                  S1 = S1 + 1
                  
                  If (CRP(1) .EQ. 0) Then       ! Case 1
                      S1 = 1
                      Exit
                  Endif
                  
                  If(((BXY2(J,K,1).EQ.PXY(I,1)).AND.       ! Case 2
     &               (BXY2(1,K,1).GT.PXY(I,1)).AND.
     &               (BXY2(J-1,K,1).GT.PXY(I,1))) 
     &                .OR. 
     &               ((BXY2(J,K,1).EQ.PXY(I,1)).AND.
     &               (BXY2(1,K,1).LT.PXY(I,1)).AND.
     &               (BXY2(J-1,K,1).LT.PXY(I,1)))) 
     &            Then
                     S1 = S1 - 1
                  Endif

              Endif
                  
              If ((BXY2(J,K,1).EQ.PXY(I,1)).AND.       
     &           (BXY2(1,K,1).EQ.PXY(I,1))) Then      ! Case 3
                  If (((PXY(I,2).LE.BXY2(J,K,2)).AND.
     &                (PXY(I,2).GE.BXY2(1,K,2)))
     &                .OR.
     &                ((PXY(I,2).GE.BXY2(J,K,2)).AND.
     &                (PXY(I,2).LE.BXY2(1,K,2))))
     &            Then
                      S1 = 1
                      Exit
                  Endif
              Endif
      
          Endif
      Enddo
    
      If (mod(S1,2) .GT. 0) Then
          IDP(I,K) = 1     
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
      
C===========================================================================C      
C        3D Immersed boundary detection      
C===========================================================================C

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
      Do 301 K = 1,KB    
      Do I = 1,IJP    

C***************************************************************************C
C             Ray-casting Method Begins                                     C
C***************************************************************************C            
      If ((PXY(I,1).LE.XMAX).AND.(PXY(I,1).GE.XMIN).AND. 
     &    (PXY(I,2).LE.YMAX).AND.(PXY(I,2).GE.YMIN).AND.
     &    (ZRP(I,K).LE.ZMAX).AND.(ZRP(I,K).GE.ZMIN)) Then 
     
      S1 = 0
      Do 302 nc = 1,NTC
      GK = (BXY3(nc,2,1)-BXY3(nc,1,1))*(BXY3(nc,3,2)-BXY3(nc,1,2))
     &    -(BXY3(nc,2,2)-BXY3(nc,1,2))*(BXY3(nc,3,1)-BXY3(nc,1,1))
      If (GK.NE.0) Then
          is_coincide = (((PXY(I,1).EQ.BXY3(nc,1,1)).AND.
     &                (PXY(I,2).EQ.BXY3(nc,1,2))).OR.
     &                ((PXY(I,1).EQ.BXY3(nc,2,1)).AND.
     &                (PXY(I,2).EQ.BXY3(nc,2,2))).OR.
     &                ((PXY(I,1).EQ.BXY3(nc,3,1)).AND.
     &                (PXY(I,2).EQ.BXY3(nc,3,2))))
          If (is_coincide) Exit

          CRP(1) = (BXY3(nc,1,1)-PXY(I,1))*(BXY3(nc,2,2)-PXY(I,2))
     &            -(BXY3(nc,1,2)-PXY(I,2))*(BXY3(nc,2,1)-PXY(I,1))
          CRP(2) = (BXY3(nc,2,1)-PXY(I,1))*(BXY3(nc,3,2)-PXY(I,2))
     &            -(BXY3(nc,2,2)-PXY(I,2))*(BXY3(nc,3,1)-PXY(I,1))
          CRP(3) = (BXY3(nc,3,1)-PXY(I,1))*(BXY3(nc,1,2)-PXY(I,2))
     &            -(BXY3(nc,3,2)-PXY(I,2))*(BXY3(nc,1,1)-PXY(I,1))

          If (((CRP(1).GE.0).AND.(CRP(2).GE.0).AND.(CRP(3).GE.0)).OR.
     &        ((CRP(1).LE.0).AND.(CRP(2).LE.0).AND.(CRP(3).LE.0))) Then
          If ((CRP(1).EQ.0).OR.(CRP(2).EQ.0).OR.(CRP(3).EQ.0)) Exit

          FL=((BXY3(nc,2,1)-BXY3(nc,1,1))*(BXY3(nc,3,2)-BXY3(nc,1,2))-
     &        (BXY3(nc,2,2)-BXY3(nc,1,2))*(BXY3(nc,3,1)-BXY3(nc,1,1))) 
          FR =(((PXY(I,2)-BXY3(nc,1,2))*(BXY3(nc,2,1)-BXY3(nc,1,1))*
     &        (BXY3(nc,3,3)-BXY3(nc,1,3)))+
     &        ((PXY(I,1)-BXY3(nc,1,1))*(BXY3(nc,2,3)-BXY3(nc,1,3))*
     &        (BXY3(nc,3,2)-BXY3(nc,1,2)))-
     &        ((PXY(I,1)-BXY3(nc,1,1))*(BXY3(nc,2,2)-BXY3(nc,1,2))*
     &        (BXY3(nc,3,3)-BXY3(nc,1,3)))-
     &        ((PXY(I,2)-BXY3(nc,1,2))*(BXY3(nc,2,3)-BXY3(nc,1,3))*
     &        (BXY3(nc,3,1)-BXY3(nc,1,1))))
          ZIB = FR/FL+BXY3(nc,1,3)
          
          If (ZRP(I,K).LT.ZIB) Then
              S1 = S1 + 1
          Else If (ZRP(I,K).EQ.ZIB) Then
              S1 = 1
              Exit
          Endif
          Endif
      Endif 
                   
302   Continue
      
      If (mod(S1,2) .GT. 0) Then
          IDP(I,K) = 1
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

      
      End Subroutine