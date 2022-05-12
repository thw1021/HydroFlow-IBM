C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (28/07/2019)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c      Caculate the coefficients for WENO schemes                           #
c                                                                           #
c############################################################################
      Subroutine WENO
      Include './Include/OCERM_INF'
      Dimension A_LOCAL(30,20), 
     &           A_TEMP(L_DEGREE,L_DEGREE)      ! Coefficient array
      Double precision A_LOCAL,A_TEMP,
     &           EDGE_FIT(4,20),CELL_ENO(20,30)

      If(IWENO .EQ. 3) NUM_STENCIL = 10        ! 3rd ENO scheme
      If(IWENO .EQ. 4) NUM_STENCIL = 12        ! 4th ENO scheme
      If(IWENO .EQ. 5) NUM_STENCIL = 20        ! 5th ENO scheme
      If(IWENO .EQ. 6) NUM_STENCIL = 28        ! 6th ENO scheme
c----------------------------------------------------------------------------c
c                 Select cells to form the stencil                           c
c----------------------------------------------------------------------------c      
      Do I = 1, IJM
         K = 0
         Do J = 1, CELL_POLYGEN(I)
            If(CELL_SIDE(I,J,2) .NE. -999) Then
              K = K + 1
              I_STENCIL(I,K) = CELL_SIDE(I,J,2)
            Endif   
         Enddo
         NUM_CLASS1 = K                        ! Norber cells with common edges
C         If(K .GE. 25) Goto 100
         Do K1 = 1, NUM_CLASS1
            If(CELL_POLYGEN(I_STENCIL(I,K1)) .EQ. 4) Then
               Do J = 1, CELL_POLYGEN(I_STENCIL(I,K1)) 
                  INOR = CELL_SIDE(I_STENCIL(I,K1),J,2)
                  If(INOR .EQ. I) Then
                     If(J .EQ. 1) INEXT = 3
                     If(J .EQ. 2) INEXT = 4
                     If(J .EQ. 3) INEXT = 1
                     If(J .EQ. 4) INEXT = 2
                  Endif
               Enddo
               I_CELL = CELL_SIDE(I_STENCIL(I,K1),INEXT,2)
               If(I_CELL .NE. -999) Then
                  IDEF = 0
                  Do K2 = 1, K
                     If(I_STENCIL(I,K2) .EQ. I_CELL) IDEF = 1
                  Enddo
                  If(IDEF .EQ. 0) Then
                     K = K + 1
                    I_STENCIL(I,K) = I_CELL
                  Endif  
               Endif    
            Else
            
            Endif
         Enddo
         NUM_CLASS2 = K                         ! Norber cells calssII
         Do K1 = 1, NUM_CLASS2
            Do J = 1, CELL_POLYGEN(I_STENCIL(I,K1))
               INOR = CELL_SIDE(I_STENCIL(I,K1),J,2)
               If(INOR .NE. -999 .AND. INOR .NE. I) Then
                  IDEF = 0
                  Do K2 = 1, K
                     If(I_STENCIL(I,K2) .EQ. INOR) IDEF = 1
                  Enddo
                  If(IDEF .EQ. 0) Then
                     K = K + 1
                     If(K .GT. 30) Goto 100
                     I_STENCIL(I,K) = INOR
                  Endif
               Endif
            Enddo
         Enddo
         NUM_CLASS3 = K                       ! Norber cells calssII  
         Do K1 = NUM_CLASS2 + 1, NUM_CLASS3
            Do J = 1, CELL_POLYGEN(I_STENCIL(I,K1))
               INOR = CELL_SIDE(I_STENCIL(I,K1),J,2)
               If(INOR .NE. -999 .AND. INOR .NE. I) Then
                  IDEF = 0
                  Do K2 = 1, K
                     If(I_STENCIL(I,K2) .EQ. INOR) IDEF = 1
                  Enddo
                  If(IDEF .EQ. 0) Then
                     K = K + 1
                     If(K .GT. 30) Goto 100
                     I_STENCIL(I,K) = INOR
                  Endif
               Endif
            Enddo
         Enddo
         NUM_CLASS4 = K                       ! Norber cells calssII  
100   Continue
      Enddo
C----------------------------------------------------------------------------C
C                Calculate the polynomal coefficients                        c
c----------------------------------------------------------------------------c
      Do I = 1, IJM
         
         SCALING_X = 0.0
         SCALING_Y = 0.0
         Do J = 1, CELL_POLYGEN(I)
            DEL_X = CELL_CUV(I,J,6) * Abs(CELL_CUV(I,J,8))
            DEL_Y = CELL_CUV(I,J,6) * Abs(CELL_CUV(I,J,7))
            If(DEL_X .GT. SCALING_X)  SCALING_X = DEL_X
            If(DEL_Y .GT. SCALING_Y)  SCALING_Y = DEL_Y
         Enddo
         SCALING_X = 1. / SCALING_X
         SCALING_Y = 1. / SCALING_Y
C         SCALING = 1.0 / Sqrt(AREA(I))
         
         Do J = 1, NUM_STENCIL

            A_LOCAL(J,1) = (CXY(I_STENCIL(I,J),1) -CXY(I,1)) * SCALING_X
            
            A_LOCAL(J,2) = (CXY(I_STENCIL(I,J),2) -CXY(I,2)) * SCALING_Y
            
            A_LOCAL(J,3) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 2) *
     &                      SCALING_X ** 2.
     
            A_LOCAL(J,4) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) * 
     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2))) * 
     &                      SCALING_X * SCALING_Y 
     
            A_LOCAL(J,5) = ((CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 2.) * 
     &                      SCALING_Y ** 2.       
     
            A_LOCAL(J,6) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 3.) *
     &                      SCALING_X ** 3.        
                           
            A_LOCAL(J,7) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 2. * 
     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2))) * 
     &                     SCALING_X ** 2. * SCALING_Y 
     
            A_LOCAL(J,8) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) * 
     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 2.) *
     &                      SCALING_X ** 1. * SCALING_Y ** 2. 
     
            A_LOCAL(J,9) = ((CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 3.) *
     &                      SCALING_Y ** 3.        
     
            A_LOCAL(J,10) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 4) *
     &                      SCALING_X ** 4. 
                           
            A_LOCAL(J,11) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 3. * 
     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2))) * 
     &                      SCALING_X ** 3. * SCALING_Y ** 1.0  
     
            A_LOCAL(J,12) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 2. * 
     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 2.) *
     &                      SCALING_X ** 2. * SCALING_Y ** 2.0 
     
            A_LOCAL(J,13) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) * 
     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 3.) * 
     &                      SCALING_X ** 1. * SCALING_Y ** 3.0
     
            A_LOCAL(J,14) = ((CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 4.) *
     &                      SCALING_Y ** 4.        

            A_LOCAL(J,15) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 5.) *
     &                      SCALING_X ** 5.        

            A_LOCAL(J,16) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 4. * 
     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2))) * 
     &                      SCALING_X ** 4. * SCALING_Y ** 1.0  

            A_LOCAL(J,17) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 3. * 
     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 2.) * 
     &                      SCALING_X ** 3. * SCALING_Y ** 2.0  

            A_LOCAL(J,18) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 2. * 
     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 3.) * 
     &                      SCALING_X ** 2. * SCALING_Y ** 3.0 

            A_LOCAL(J,19) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 1. * 
     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 4.) * 
     &                      SCALING_X ** 1. * SCALING_Y ** 4.0  

            A_LOCAL(J,20) = ((CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 5.) *
     &                      SCALING_Y ** 5.       
           
C            A_LOCAL(J,1) = (CXY(I_STENCIL(I,J),1) -CXY(I,1)) * SCALING
            
C            A_LOCAL(J,2) = (CXY(I_STENCIL(I,J),2) -CXY(I,2)) * SCALING
            
C            A_LOCAL(J,3) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 2) *
C     &                      SCALING ** 2.
     
C            A_LOCAL(J,4) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) * 
C     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2))) * 
C     &                      SCALING ** 2. 
     
C            A_LOCAL(J,5) = ((CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 2.) * 
C     &                      SCALING ** 2.       
     
C            A_LOCAL(J,6) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 3.) *
C     &                      SCALING ** 3.        
                           
C            A_LOCAL(J,7) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 2. * 
C     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2))) * 
C     &                     SCALING ** 3. 
     
C            A_LOCAL(J,8) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) * 
C     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 2.) *
C     &                      SCALING ** 3. 
     
C            A_LOCAL(J,9) = ((CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 3.) *
C     &                      SCALING ** 3.        
     
C            A_LOCAL(J,10) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 4) *
C     &                      SCALING ** 4. 
                           
C            A_LOCAL(J,11) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 3. * 
C     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2))) * 
C     &                      SCALING ** 4.  
     
C            A_LOCAL(J,12) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 2. * 
C     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 2.) *
C     &                      SCALING ** 4. 
     
C            A_LOCAL(J,13) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) * 
C     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 3.) * 
C     &                      SCALING ** 4. 
     
C            A_LOCAL(J,14) = ((CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 4.) *
C     &                      SCALING ** 4.        

C            A_LOCAL(J,15) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 5.) *
C     &                      SCALING ** 5.        

C            A_LOCAL(J,16) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 4. * 
C     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2))) * 
C     &                      SCALING ** 5.  

C            A_LOCAL(J,17) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 3. * 
C     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 2.) * 
C     &                      SCALING ** 5.  

C            A_LOCAL(J,18) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 2. * 
C     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 3.) * 
C     &                      SCALING ** 5.  

C            A_LOCAL(J,19) = ((CXY(I_STENCIL(I,J),1) -CXY(I,1)) ** 1. * 
C     &                     (CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 4.) * 
C     &                      SCALING ** 5.  

C            A_LOCAL(J,20) = ((CXY(I_STENCIL(I,J),2) -CXY(I,2)) ** 5.) *
C     &                      SCALING ** 5.        
     
         Enddo

C-----   calculate the Trans(A)*A  householder transformation: A=QR
         Do J = 1, L_DEGREE
            Do K = 1, L_DEGREE
               A_TEMP(J,K) = 0.0
               Do IT = 1, NUM_STENCIL
                 A_TEMP(J,K) = A_TEMP(J,K)+A_LOCAL(IT,J)*A_LOCAL(IT,K)
               Enddo
            Enddo
         Enddo
C-----   calculate the INV(Trans(A)*A)
         Call BRINV(A_TEMP,L_DEGREE,L) 
         Do IH = 1, L_DEGREE
            Do IV = 1, NUM_STENCIL
               CELL_ENO(IH,IV) = 0.0
               Do IT = 1, L_DEGREE
                  CELL_ENO(IH,IV) = CELL_ENO(IH,IV) + 
     &                             A_TEMP(IH,IT) * A_LOCAL(IV,IT)
               Enddo 
            Enddo
         Enddo 
C-----   Recalculate the interplation coefficients         
         Do K = 1, NUM_STENCIL
            CELL_ENO(1,K) = CELL_ENO(1,K) * SCALING_X
            CELL_ENO(2,K) = CELL_ENO(2,K) * SCALING_Y
            CELL_ENO(3,K) = CELL_ENO(3,K) * SCALING_X** 2.
            CELL_ENO(4,K) = CELL_ENO(4,K) * SCALING_X*SCALING_Y
            CELL_ENO(5,K) = CELL_ENO(5,K) * SCALING_Y**2.
            CELL_ENO(6,K) = CELL_ENO(6,K) * SCALING_X**3.
            CELL_ENO(7,K) = CELL_ENO(7,K) * SCALING_X**2.* SCALING_Y
            CELL_ENO(8,K) = CELL_ENO(8,K) * SCALING_X**1.* SCALING_Y**2.
            CELL_ENO(9,K) = CELL_ENO(9,K) * SCALING_Y**3.
            CELL_ENO(10,K) = CELL_ENO(10,K)*SCALING_X**4.
            CELL_ENO(11,K) = CELL_ENO(11,K)*SCALING_X**3.* SCALING_Y
            CELL_ENO(12,K) = CELL_ENO(12,K)*SCALING_X**2.* SCALING_Y**2.
            CELL_ENO(13,K) = CELL_ENO(13,K)*SCALING_X**1.* SCALING_Y**3.
            CELL_ENO(14,K) = CELL_ENO(14,K)*SCALING_Y**4.
            CELL_ENO(15,K) = CELL_ENO(10,K)*SCALING_X**5.
            CELL_ENO(16,K) = CELL_ENO(11,K)*SCALING_X**4.* SCALING_Y
            CELL_ENO(17,K) = CELL_ENO(12,K)*SCALING_X**3.* SCALING_Y**2.
            CELL_ENO(18,K) = CELL_ENO(13,K)*SCALING_X**2.* SCALING_Y**3.
            CELL_ENO(19,K) = CELL_ENO(14,K)*SCALING_X**1.* SCALING_Y**4.
            CELL_ENO(20,K) = CELL_ENO(14,K)*SCALING_Y**5.
C            CELL_ENO(1,K) = CELL_ENO(1,K) * SCALING
C            CELL_ENO(2,K) = CELL_ENO(2,K) * SCALING
C            CELL_ENO(3,K) = CELL_ENO(3,K) * SCALING ** 2.
C            CELL_ENO(4,K) = CELL_ENO(4,K) * SCALING ** 2.
C            CELL_ENO(5,K) = CELL_ENO(5,K) * SCALING ** 2.
C            CELL_ENO(6,K) = CELL_ENO(6,K) * SCALING ** 3.
C            CELL_ENO(7,K) = CELL_ENO(7,K) * SCALING ** 3.
C            CELL_ENO(8,K) = CELL_ENO(8,K) * SCALING ** 3.
C            CELL_ENO(9,K) = CELL_ENO(9,K) * SCALING ** 3.
C            CELL_ENO(10,K) = CELL_ENO(10,K) * SCALING ** 4.
C            CELL_ENO(11,K) = CELL_ENO(11,K) * SCALING ** 4.
C            CELL_ENO(12,K) = CELL_ENO(12,K) * SCALING ** 4.
C            CELL_ENO(13,K) = CELL_ENO(13,K) * SCALING ** 4.
C            CELL_ENO(14,K) = CELL_ENO(14,K) * SCALING ** 4.
         Enddo  
c------  shape function interplation along the cell edge
         X0 = CXY(I,1)
         Y0 = CXY(I,2)
         Do J = 1, CELL_POLYGEN(I)
            X1 = PXY(CELL_SIDE(I,J,3),1)
            Y1 = PXY(CELL_SIDE(I,J,3),2)
            EDGE_FIT(J,1) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),1,0)
            EDGE_FIT(J,2) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),0,1)
            EDGE_FIT(J,3) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),2,0)
            EDGE_FIT(J,4) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),1,1)
            EDGE_FIT(J,5) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),0,2)
            EDGE_FIT(J,6) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),3,0)
            EDGE_FIT(J,7) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),2,1)
            EDGE_FIT(J,8) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),1,2)
            EDGE_FIT(J,9) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),0,3)
            EDGE_FIT(J,10) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),4,0)
            EDGE_FIT(J,11) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),3,1)
            EDGE_FIT(J,12) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),2,2)
            EDGE_FIT(J,13) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),1,3)
            EDGE_FIT(J,14) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),0,4)
         Enddo
C-----   Forming the interplation coefficients        
         Do J = 1, CELL_POLYGEN(I)
            A_ENO(I,J,0) = 0.0 
            Do K = 1, NUM_STENCIL
               A_ENO(I,J,K) = 0.0
               Do I1 = 1, L_DEGREE
                   A_ENO(I,J,K) = A_ENO(I,J,K) + 
     &                            CELL_ENO(I1,K) * EDGE_FIT(J,I1) /
     &                            CELL_CUV(I,J,6)
               Enddo
               A_ENO(I,J,0) = A_ENO(I,J,0) + A_ENO(I,J,K)  
            Enddo
            A_ENO(I,J,0) = 1.0 - A_ENO(I,J,0)
         Enddo 
      Enddo
      PRINT*, (A_ENO(16760,3,K),K=0,10)
c----------------------------------------------------------------------------c
c                Calculate the boundary integration                          c
c----------------------------------------------------------------------------c


c----------------------------------------------------------------------------c          
      Return
      End
      
      Function FIT(X1,Y1,X0,Y0,DS,COS_ALF,SIN_ALF,L1,L2)
     	Dimension T(5),C(5)
      
      Double precision FIT,X1,Y1,X0,Y0,DS,COS_ALF,SIN_ALF
	Double precision A,B,F,G,T,C,S,P,H,AA,BB,W,X,Q
      Integer L1,L2
	Data T/-0.9061798459,-0.5384693101,0.0,
     *          0.5384693101,0.9061798459/
	Data C/0.2369268851,0.4786286705,0.5688888889,
     *         0.4786286705,0.2369268851/
      A = 0.0
      B = DS
      EPS = 1.E-6
	M=1
	S=(B-A)*0.001
	P=0.0
10	H=(B-A)/M
	G=0.0
	Do 30 I=1,M
	  AA=A+(I-1)*H
	  BB=A+I*H
	  W=0.0
	  Do 20 J=1,5
	    X=((BB-AA)*T(J)+(BB+AA))/2.0
	    W=W+F(X1,Y1,X0,Y0,COS_ALF,SIN_ALF,L1,L2,X)*C(J)
20	  Continue
	  G=G+W
30	Continue
	G=G*H/2.0
	Q=Abs(G-P)/(1.0+Abs(G))
	If ((Q.GE.EPS).AND.(Abs(H).GT.Abs(S))) Then
	  P=G
	  M=M+1
	  Goto 10
	Endif
      FIT = G
      Return
      End
      
      Function F(X1,Y1,X0,Y0,COS_ALF,SIN_ALF,L1,L2,X)
      Double precision F,X1,Y1,X0,Y0,COS_ALF,SIN_ALF,X
      Integer L1,L2
      F = (X1-X0-X*SIN_ALF)**L1 * (Y1-Y0+X*COS_ALF) **L2 
      Return
      End
