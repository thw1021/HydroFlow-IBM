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
C     &           INDEX_GHOST(-10000:-1),GHOSTP(-10000:-1,2)
      Double precision A_LOCAL,A_TEMP,
     &           EDGE_FIT(4,30),CELL_ENO(20,30)
      Dimension I1_TEM(2),I2_TEM(2)
      Dimension INDEX_POINT_G(IJP,20)
      Dimension POINT_CELL_G(IJP)
     
      Do I = 1, IJP
         POINT_CELL_G(I) = POINT_CELL(I)
         Do J = 1, POINT_CELL(I)
            INDEX_POINT_G(I,J) = INDEX_POINT(I,J)
         Enddo
      Enddo
      
      If(IWENO .EQ. 3) NUM_STENCIL = 8         ! 3rd ENO scheme
      If(IWENO .EQ. 4) NUM_STENCIL = 15        ! 4th ENO scheme
      If(IWENO .EQ. 5) NUM_STENCIL = 20        ! 5th ENO scheme
      If(IWENO .EQ. 6) NUM_STENCIL = 28        ! 6th ENO scheme

C----------------------------------------------------------------------------C
C           Generate ghost cells for the higher order schemes                c
c----------------------------------------------------------------------------c
       NUM_GHOST = 0
       Do I = 1, IJM
          Do I1 = 1, 2
             I1_TEM(I1) = 0
             I2_TEM(I1) = 0
          Enddo
          KKK = 0
          Do J = 1, CELL_POLYGEN(I)
             If(CELL_SIDE(I,J,2) .EQ. -999) Then
                K1 = CELL_SIDE(I,J,3)
                K2 = CELL_SIDE(I,J,4)
                X_MID = (PXY(K1,1) + PXY(K2,1)) / 2.
                Y_MID = (PXY(K1,2) + PXY(K2,2)) / 2.
                NUM_GHOST = NUM_GHOST - 1
                If(NUM_GHOST .LT. -10000) Then
                   Print*, 'OVER THE UP LIMIT OF THE ARRAY'
                   Stop
                Endif   
                INDEX_GHOST(NUM_GHOST) = I
                GHOSTP(NUM_GHOST,1) = 2. * X_MID - CXY(I,1)
                GHOSTP(NUM_GHOST,2) = 2. * Y_MID - CXY(I,2)
 
                KKK = KKK + 1                  ! New cell in the corner
                I1_TEM(KKK) = NUM_GHOST
                I2_TEM(KKK) = 1
C                If(KKK .EQ. 2) Then
C                   If(K1 .EQ. K3) IPOINT_COM = K1
C                   If(K1 .EQ. K4) IPOINT_COM = K1
C                   If(K2 .EQ. K3) IPOINT_COM = K2
C                   If(K2 .EQ. K4) IPOINT_COM = K2
C                Endif
C                K3 = K1
C                K4 = K2
                
                POINT_CELL_G(CELL_SIDE(I,J,3)) =
     &                      POINT_CELL_G(CELL_SIDE(I,J,3)) + 1
                NEW_VER_CELL = POINT_CELL_G(CELL_SIDE(I,J,3))
                INDEX_POINT_G(CELL_SIDE(I,J,3),NEW_VER_CELL) = NUM_GHOST
                
                POINT_CELL_G(CELL_SIDE(I,J,4)) =
     &                      POINT_CELL_G(CELL_SIDE(I,J,4)) + 1
                NEW_VER_CELL = POINT_CELL_G(CELL_SIDE(I,J,4))
                INDEX_POINT_G(CELL_SIDE(I,J,4),NEW_VER_CELL) = NUM_GHOST
                
             Endif
          Enddo
          If(I2_TEM(1) * I2_TEM(2) .EQ. 1) Then
            X_MID = (GHOSTP(I1_TEM(1),1) + GHOSTP(I1_TEM(2),1)) / 2.
            Y_MID = (GHOSTP(I1_TEM(1),2) + GHOSTP(I1_TEM(2),2)) / 2.
            NUM_GHOST = NUM_GHOST - 1
            INDEX_GHOST(NUM_GHOST) = I
            GHOSTP(NUM_GHOST,1) = 2. * X_MID - CXY(I,1)
            GHOSTP(NUM_GHOST,2) = 2. * Y_MID - CXY(I,2)
            
C            POINT_CELL_G(IPOINT_COM) = POINT_CELL_G(IPOINT_COM) + 1
C            NEW_VER_CELL = POINT_CELL_G(IPOINT_COM)
C            INDEX_POINT_G(IPOINT_COM,NEW_VER_CELL) = NUM_GHOST
          Endif
       Enddo
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
         Do K1 = -1, NUM_GHOST, -1
            If(INDEX_GHOST(K1) .EQ. I) Then
               IDEF = 0
               Do K3 = 1, K
                  If(K1 .EQ. I_STENCIL(I,K3)) IDEF = IDEF + 1
               Enddo
               If(IDEF .EQ. 0) Then
                  K = K + 1
                  I_STENCIL(I,K) = K1
               Endif   
            Endif
         Enddo      
         NUM_CLASS1 = K                        ! Norber cells with common edges
C         If(K .GE. 25) Goto 100
         Do J = 1, CELL_POLYGEN(I)
            IPOINT = CELL_SIDE(I,J,3)
            Do ICELL = 1, POINT_CELL_G(IPOINT) 
               ITEMP = INDEX_POINT_G(IPOINT,ICELL)
               ITRUE = 1
               Do L = 1, K
                  If(ITEMP .EQ. I_STENCIL(I,L)) ITRUE = 0  
               Enddo
               If(ITRUE .EQ. 1 .AND. ITEMP .NE. I) Then
                  K = K + 1
                  I_STENCIL(I,K) = ITEMP
               Endif
            Enddo
         Enddo
C         Do K1 = -1, NUM_GHOST, -1
C            Do K2 = 1, K
C               If(INDEX_GHOST(K1) .EQ. I_STENCIL(I,K2)) Then
C                  K = K + 1
C                  I_STENCIL(I,K) = K1
C               Endif
C            Enddo   
C         Enddo
         NUM_CLASS2 = K                        ! Norber cells with common vertex
C         If(K .GE. 25) Goto 100
         Do J = 1, NUM_CLASS1
            II = I_STENCIL(I,J)
            If(II .GT. 0) Then
              Do JJ = 1, CELL_POLYGEN(II)
                 ITEMP = CELL_SIDE(II,JJ,2)
                 If(ITEMP .NE. -999) Then
                    ITRUE = 1
                    Do L = 1, K
                       If(ITEMP .EQ. I_STENCIL(I,L)) ITRUE = 0
                    Enddo
                    If(ITRUE .EQ. 1 .AND. ITEMP .NE. I) Then
                       K = K + 1
                       If(K .GE. 30) Goto 100
                       I_STENCIL(I,K) = ITEMP
                    Endif
                 Endif
              Enddo
            Endif  
         Enddo
         NUM_CLASS3 = K                         !Norber cells to class1
C         If(K .GE. 25) Goto 100
         Do J = NUM_CLASS1 + 1, NUM_CLASS2
            II = I_STENCIL(I,J)
            If(II .GT. 0) Then
              Do JJ = 1, CELL_POLYGEN(II)
                 ITEMP = CELL_SIDE(II,JJ,2)
                 If(ITEMP .NE. -999) Then
                    ITRUE = 1
                    Do L = 1, K
                       If(ITEMP .EQ. I_STENCIL(I,L)) ITRUE = 0
                    Enddo
                    If(ITRUE .EQ. 1 .AND. ITEMP .NE. I) Then
                       K = K + 1
                       If(K .GE. 30) Goto 100
                       I_STENCIL(I,K) = ITEMP
                    Endif
                 Endif
              Enddo
            Endif   
         Enddo
         Do K1 = -1, NUM_GHOST, -1
            Do K2 = 1, K
               If(INDEX_GHOST(K1) .EQ. I_STENCIL(I,K2)) Then
                  IDEF = 0
                  Do K3 = 1, K
                     If(K1 .EQ. I_STENCIL(I,K3)) IDEF = IDEF + 1
                  Enddo
                  If(IDEF .EQ. 0) Then      
                     K = K + 1
                     If(K .GE. 30) Goto 100
                     I_STENCIL(I,K) = K1
                  Endif  
               Endif
            Enddo   
         Enddo
         NUM_CLASS4 = K                        !Norber cells to class2
C         If(K .GE. 25) Goto 100 
         Do J = NUM_CLASS2 +1 , NUM_CLASS4
            II = I_STENCIL(I,J)
            If(II .GT. 0) Then
              Do JJ = 1, CELL_POLYGEN(II)
                 ITEMP = CELL_SIDE(II,JJ,2)
                 If(ITEMP .NE. -999) Then
                    ITRUE = 1
                    Do L = 1, K
                       If(ITEMP .EQ. I_STENCIL(I,L)) ITRUE = 0
                    Enddo
                    If(ITRUE .EQ. 1 .AND. ITEMP .NE. I) Then
                       K = K + 1
                       If(K .GE. 30) Goto 100
                       I_STENCIL(I,K) = ITEMP
                    Endif
                 Endif
              Enddo
            Endif  
         Enddo
         Do K1 = -1, NUM_GHOST, -1
            Do K2 = 1, K
               If(INDEX_GHOST(K1) .EQ. I_STENCIL(I,K2)) Then
                  IDEF = 0
                  Do K3 = 1, K
                     If(K1 .EQ. I_STENCIL(I,K3)) IDEF = IDEF + 1
                  Enddo
                  If(IDEF .EQ. 0) Then      
                     K = K + 1
                     If(K .GE. 30) Goto 100
                     I_STENCIL(I,K) = K1
                  Endif  
               Endif
            Enddo   
         Enddo
         NUM_CLASS5 = K                       !Norber cells to class3+class4
         Do J = NUM_CLASS4 +1 , NUM_CLASS5
            II = I_STENCIL(I,J)
            If(II .GT. 0) Then
              Do JJ = 1, CELL_POLYGEN(II)
                 ITEMP = CELL_SIDE(II,JJ,2)
                 If(ITEMP .NE. -999) Then
                    ITRUE = 1
                    Do L = 1, K
                       If(ITEMP .EQ. I_STENCIL(I,L)) ITRUE = 0
                    Enddo
                    If(ITRUE .EQ. 1 .AND. ITEMP .NE. I) Then
                       K = K + 1
                       If(K .GE. 30) Goto 100
                       I_STENCIL(I,K) = ITEMP
                    Endif
                 Endif
              Enddo
            Endif  
         Enddo
         Do K1 = -1, NUM_GHOST, -1
            Do K2 = 1, K
               If(INDEX_GHOST(K1) .EQ. I_STENCIL(I,K2)) Then
                  IDEF = 0
                  Do K3 = 1, K
                     If(K1 .EQ. I_STENCIL(I,K3)) IDEF = IDEF + 1
                  Enddo
                  If(IDEF .EQ. 0) Then      
                     K = K + 1
                     If(K .GE. 30) Goto 100
                     I_STENCIL(I,K) = K1
                  Endif  
               Endif
            Enddo   
         Enddo
         NUM_CLASS6 = K                       !Norber cells to class3+class5
         Do J = NUM_CLASS5 +1 , NUM_CLASS6
            II = I_STENCIL(I,J)
            If(II .GT. 0) Then
              Do JJ = 1, CELL_POLYGEN(II)
                 ITEMP = CELL_SIDE(II,JJ,2)
                 If(ITEMP .NE. -999) Then
                    ITRUE = 1
                    Do L = 1, K
                       If(ITEMP .EQ. I_STENCIL(I,L)) ITRUE = 0
                    Enddo
                    If(ITRUE .EQ. 1 .AND. ITEMP .NE. I) Then
                       K = K + 1
                       If(K .GE. 30) Goto 100
                       I_STENCIL(I,K) = ITEMP
                    Endif
                 Endif
              Enddo
            Endif  
         Enddo
         Do K1 = -1, NUM_GHOST, -1
            Do K2 = 1, K
               If(INDEX_GHOST(K1) .EQ. I_STENCIL(I,K2)) Then
                  IDEF = 0
                  Do K3 = 1, K
                     If(K1 .EQ. I_STENCIL(I,K3)) IDEF = IDEF + 1
                  Enddo
                  If(IDEF .EQ. 0) Then      
                     K = K + 1
                     If(K .GE. 30) Goto 100
                     I_STENCIL(I,K) = K1
                  Endif  
               Endif
            Enddo   
         Enddo
         NUM_CLASS7 = K                       !Norber cells to class3+class4
         Do J = NUM_CLASS6 +1 , NUM_CLASS7
            II = I_STENCIL(I,J)
            If(II .GT. 0) Then
              Do JJ = 1, CELL_POLYGEN(II)
                 ITEMP = CELL_SIDE(II,JJ,2)
                 If(ITEMP .NE. -999) Then
                    ITRUE = 1
                    Do L = 1, K
                       If(ITEMP .EQ. I_STENCIL(I,L)) ITRUE = 0
                    Enddo
                    If(ITRUE .EQ. 1 .AND. ITEMP .NE. I) Then
                       K = K + 1
                       If(K .GE. 30) Goto 100
                       I_STENCIL(I,K) = ITEMP
                    Endif
                 Endif
              Enddo
           Endif   
         Enddo
         Do K1 = -1, NUM_GHOST, -1
            Do K2 = 1, K
               If(INDEX_GHOST(K1) .EQ. I_STENCIL(I,K2)) Then
                  IDEF = 0
                  Do K3 = 1, K
                     If(K1 .EQ. I_STENCIL(I,K3)) IDEF = IDEF + 1
                  Enddo
                  If(IDEF .EQ. 0) Then      
                     K = K + 1
                     If(K .GE. 30) Goto 100
                     I_STENCIL(I,K) = K1
                  Endif  
               Endif
            Enddo   
         Enddo
         NUM_CLASS8 = K                       !Norber cells to class3+class4
100   Continue
      Enddo
C      PRINT*, (I_STENCIL(16481,K),K=1,NUM_STENCIL)
C      STOP
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
            If(I_STENCIL(I,J) .GT. 0) Then
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
            Else
            A_LOCAL(J,1) = (GHOSTP(I_STENCIL(I,J),1)-CXY(I,1))*SCALING_X
            
            A_LOCAL(J,2) = (GHOSTP(I_STENCIL(I,J),2)-CXY(I,2))*SCALING_Y
            
            A_LOCAL(J,3) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) ** 2) *
     &                      SCALING_X ** 2.
     
            A_LOCAL(J,4) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) * 
     &                     (GHOSTP(I_STENCIL(I,J),2)-CXY(I,2))) * 
     &                      SCALING_X * SCALING_Y 
     
            A_LOCAL(J,5) = ((GHOSTP(I_STENCIL(I,J),2)-CXY(I,2)) ** 2.)* 
     &                      SCALING_Y ** 2.       
     
            A_LOCAL(J,6) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) ** 3.) *
     &                      SCALING_X ** 3.        
                           
            A_LOCAL(J,7) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) ** 2. * 
     &                     (GHOSTP(I_STENCIL(I,J),2)-CXY(I,2))) * 
     &                     SCALING_X ** 2. * SCALING_Y 
     
            A_LOCAL(J,8) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) * 
     &                     (GHOSTP(I_STENCIL(I,J),2)-CXY(I,2)) ** 2.) *
     &                      SCALING_X ** 1. * SCALING_Y ** 2. 
     
            A_LOCAL(J,9) = ((GHOSTP(I_STENCIL(I,J),2)-CXY(I,2)) ** 3.) *
     &                      SCALING_Y ** 3.        
     
            A_LOCAL(J,10) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) ** 4) *
     &                      SCALING_X ** 4. 
                           
            A_LOCAL(J,11) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) ** 3. * 
     &                     (GHOSTP(I_STENCIL(I,J),2)-CXY(I,2))) * 
     &                      SCALING_X ** 3. * SCALING_Y ** 1.0  
     
            A_LOCAL(J,12) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) ** 2. * 
     &                     (GHOSTP(I_STENCIL(I,J),2)-CXY(I,2)) ** 2.) *
     &                      SCALING_X ** 2. * SCALING_Y ** 2.0 
     
            A_LOCAL(J,13) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) * 
     &                     (GHOSTP(I_STENCIL(I,J),2)-CXY(I,2)) ** 3.) * 
     &                      SCALING_X ** 1. * SCALING_Y ** 3.0
     
            A_LOCAL(J,14) = ((GHOSTP(I_STENCIL(I,J),2)-CXY(I,2)) ** 4.)*
     &                      SCALING_Y ** 4.        

            A_LOCAL(J,15) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) ** 5.)*
     &                      SCALING_X ** 5.        

            A_LOCAL(J,16) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) ** 4. * 
     &                     (GHOSTP(I_STENCIL(I,J),2) -CXY(I,2))) * 
     &                      SCALING_X ** 4. * SCALING_Y ** 1.0  

            A_LOCAL(J,17) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) ** 3. * 
     &                     (GHOSTP(I_STENCIL(I,J),2)-CXY(I,2)) ** 2.) * 
     &                      SCALING_X ** 3. * SCALING_Y ** 2.0  

            A_LOCAL(J,18) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) ** 2. * 
     &                     (GHOSTP(I_STENCIL(I,J),2)-CXY(I,2)) ** 3.) * 
     &                      SCALING_X ** 2. * SCALING_Y ** 3.0 

            A_LOCAL(J,19) = ((GHOSTP(I_STENCIL(I,J),1)-CXY(I,1)) ** 1. * 
     &                     (GHOSTP(I_STENCIL(I,J),2)-CXY(I,2)) ** 4.) * 
     &                      SCALING_X ** 1. * SCALING_Y ** 4.0  

            A_LOCAL(J,20) = ((GHOSTP(I_STENCIL(I,J),2)-CXY(I,2)) ** 5.)*
     &                      SCALING_Y ** 5. 
            Endif      
           
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
            CELL_ENO(15,K) = CELL_ENO(15,K)*SCALING_X**5.
            CELL_ENO(16,K) = CELL_ENO(16,K)*SCALING_X**4.* SCALING_Y
            CELL_ENO(17,K) = CELL_ENO(17,K)*SCALING_X**3.* SCALING_Y**2.
            CELL_ENO(18,K) = CELL_ENO(18,K)*SCALING_X**2.* SCALING_Y**3.
            CELL_ENO(19,K) = CELL_ENO(19,K)*SCALING_X**1.* SCALING_Y**4.
            CELL_ENO(20,K) = CELL_ENO(20,K)*SCALING_Y**5.
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
            EDGE_FIT(J,15) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),5,0)
            EDGE_FIT(J,16) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),4,1)
            EDGE_FIT(J,17) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),3,2)
            EDGE_FIT(J,18) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),2,3)
            EDGE_FIT(J,19) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),1,4)
            EDGE_FIT(J,20) = FIT(X1,Y1,X0,Y0,CELL_CUV(I,J,6),
     &                        CELL_CUV(I,J,7),CELL_CUV(I,J,8),0,5)
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
C      PRINT*, (A_ENO(16481,3,K),K=0,NUM_STENCIL)
C      PRINT*, 'OK'
C      PRINT*, (A_ENO(16560,3,K),K=0,NUM_STENCIL)
C      STOP
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
