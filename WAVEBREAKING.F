C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (01/05/2016)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                Wave breaking model                                        #
c                                                                           #
c############################################################################
      Subroutine WAVEBREAKING
      Include './Include/OCERM_INF'
	Common/DYNBLK/AS(IJM,KB,IPOLYGEN),AB(IJM,KB),AT(IJM,KB),
     &	          AP(IJM,KB),BB(IJM,KB),X(IJM,KB),XINI(IJM,KB)
      Parameter(WB_ANGLE_I = 2.*3.14*30./360.,
     &           WB_ANGLE_E = 2.*3.14*10./360.)
      Dimension WB_MASK2(IJM)
      Do I = 1, IJM
         Do K = 1, KBM
            VIS_BW(I,K) = 0.0
         Enddo
      Enddo
      Do I = 1, IJM
         If(CCM(I) .EQ. 1.0) Then
            Do I1 = 1, N_SOURCE
	         II = IGEN(I1)
	         If(I .EQ. II) Goto 10001
	      Enddo
C---------- the local still water depth
            If(HC(I) .GT. 0.0) Then
               DD = HC(I)
            Else
               DD = DC(I)
            Endif
c---------- the flow direction   
            WAVE_DIR_X = UAVE(I) / Sqrt(UAVE(I)**2.+VAVE(I)**2.)
            WAVE_DIR_Y = VAVE(I) / Sqrt(UAVE(I)**2.+VAVE(I)**2.)
c---------- define the wave crest
            SCOPE = Dmax1(DD,3.*Sqrt(AREA(I)))
c            SCOPE = 2. * Sqrt(AREA(I))
            Do I2 = 1, IJM
               DIS = Sqrt((CXY(I2,1) - CXY(I,1)) ** 2. +
     &                    (CXY(I2,2) - CXY(I,2)) ** 2.)
               If(DIS .LE. SCOPE .AND. I2 .NE. I) Then
                  CELL_DIR_X = (CXY(I2,1) - CXY(I,1)) / DIS
                  CELL_DIR_Y = (CXY(I2,2) - CXY(I,2)) / DIS
                  CROSS_DIR = WAVE_DIR_X * CELL_DIR_X + 
     &                        WAVE_DIR_Y * CELL_DIR_Y
                  If(Abs(CROSS_DIR) .GE. 0.90) Then
                     If(ELF(I2) .GT. ELF(I)) Goto 1001
                  Endif
               Endif
            Enddo
C----------  wave breaking model 
            SCOPE = Dmax1(2.0 * DD,5.*Sqrt(AREA(I)))
            SLOPE_MAX = 0.0
            WB_LENGTH_MAX = 0.0
            WB_LENGTH_MIN = 2. * DD
            Do I2 = 1, IJM
               DIS = Sqrt((CXY(I2,1) - CXY(I,1)) ** 2. +
     &                    (CXY(I2,2) - CXY(I,2)) ** 2.)
               If(DIS .LE. SCOPE .AND. I2 .NE. I) Then
                  CELL_DIR_X = (CXY(I2,1) - CXY(I,1)) / DIS
                  CELL_DIR_Y = (CXY(I2,2) - CXY(I,2)) / DIS
                  CROSS_DIR = WAVE_DIR_X * CELL_DIR_X + 
     &                        WAVE_DIR_Y * CELL_DIR_Y
                  If(CROSS_DIR .GE. 0.90) Then
                     WB_SLOPE_X = 0.0
                     WB_SLOPE_Y = 0.0
                     WB_SLOPE = 0.0
                     Do J = 1, CELL_POLYGEN(I2)
		                If(CFM(CELL_SIDE(I2,J,1)) .EQ. 1.0) Then
                           WB_SLOPE_X = WB_SLOPE_X + 
     &                       WIX(I2,J)*(ELF(CELL_SIDE(I2,J,2))-ELF(I2))
                           WB_SLOPE_Y = WB_SLOPE_Y + 
     &                       WIY(I2,J)*(ELF(CELL_SIDE(I2,J,2))-ELF(I2))
                        Endif
                     Enddo
                     WB_SLOPE = Sqrt(WB_SLOPE_X**2. + WB_SLOPE_Y**2.) 
                     If(WB_SLOPE .GT. SLOPE_MAX) Then
                        SLOPE_MAX = WB_SLOPE
                        I_WB = I2
                     Endif 
                     If(WB_SLOPE .GE. Tan(WB_ANGLE_I)) Then
                        DIS = Sqrt((CXY(I2,1) - CXY(I,1)) ** 2. +
     &                             (CXY(I2,2) - CXY(I,2)) ** 2.)
                        If(DIS .GT. WB_LENGTH_MAX) WB_LENGTH_MAX = DIS
                        If(DIS .LT. WB_LENGTH_MIN) WB_LENGTH_MIN = DIS
                     Endif 
                  Endif
               Endif
            Enddo
C----------  wave breaking when the local slope is enough steep
            If(SLOPE_MAX .GE. Tan(WB_ANGLE_I) ) Then 
               WB_SCOPE = WB_LENGTH_MAX - WB_LENGTH_MIN 
               WB_SCOPE = Dmax1(WB_SCOPE,2.*Sqrt(AREA(I)))
               If(WB_SCOPE .LT. 0.40 * DD .AND. 
     &            WB_SCOPE .GE. 0.35 * DD) Then
                  Do I2 = 1, IJM
                     DIS = Sqrt((CXY(I2,1) - CXY(I,1)) ** 2. +
     &                          (CXY(I2,2) - CXY(I,2)) ** 2.)
                     If(DIS .LE. SCOPE .AND. I2 .NE. I) Then
                        CELL_DIR_X = (CXY(I2,1) - CXY(I,1)) / DIS
                        CELL_DIR_Y = (CXY(I2,2) - CXY(I,2)) / DIS
                        CROSS_DIR = WAVE_DIR_X * CELL_DIR_X + 
     &                              WAVE_DIR_Y * CELL_DIR_Y
                        If(CROSS_DIR .GE. 0.90) Then
                           WB_SLOPE_X = 0.0
                           WB_SLOPE_Y = 0.0
                           WB_SLOPE = 0.0
                           Do J = 1, CELL_POLYGEN(I2)
		                      If(CFM(CELL_SIDE(I2,J,1)) .EQ. 1.0) Then
                                 WB_SLOPE_X = WB_SLOPE_X + 
     &                        WIX(I2,J)*(ELF(CELL_SIDE(I2,J,2))-ELF(I2))
                                 WB_SLOPE_Y = WB_SLOPE_Y + 
     &                        WIY(I2,J)*(ELF(CELL_SIDE(I2,J,2))-ELF(I2))
                              Endif
                           Enddo
c-----        local slope                      
                           WB_SLOPE=Sqrt(WB_SLOPE_X**2.+WB_SLOPE_Y**2.) 
c-----        wave breaking zone
                           If(WB_SLOPE .GE. Tan(WB_ANGLE_I)) Then 
c                           If(DIS .GE. WB_LENGTH_MIN .AND.
c     &                        DIS .LT. WB_LENGTH_MAX) Then 
                              Do K = 1, KBM
                                 Do J = 1, CELL_POLYGEN(I2)
		                            If(CFM(CELL_SIDE(I2,J,1)) .EQ. 1.0) 
     &                                 AS(I2,K,J) = 0.0
                                 Enddo
                                 AP(I2,K) = 1.0
	                           AB(I2,K) = 0.0
	                           AT(I2,K) = 0.0
	                           BB(I2,K) = 0.0
C                  VIS_BW(I2,K) = 0.1 * (DC(I2) * (1. + ZZ(K))) ** 2. *
C     &                          Abs(WB) / DC(I2)

                                 C(I2,K) = 1.0
                             Enddo
                             WB_MASK2(I2) = 1.0
                          Endif
                       Endif
                     Endif
                  Enddo
               Else If(WB_SCOPE .GE. 0.4 * DD) Then
c                If(WB_LENGTH_MIN/WB_LENGTH_MAX .GT. 0.2) 
c     &             WB_UP = WB_LENGTH_MIN
c                If(WB_LENGTH_MIN/WB_LENGTH_MAX .LE. 0.2) 
c     &             WB_UP = 0.0
                 Do I2 = 1, IJM
                   DIS = Sqrt((CXY(I2,1) - CXY(I,1)) ** 2. +
     &                        (CXY(I2,2) - CXY(I,2)) ** 2.)
                   If(DIS .LE. SCOPE .AND. I2 .NE. I) Then
                      CELL_DIR_X = (CXY(I2,1) - CXY(I,1)) / DIS
                      CELL_DIR_Y = (CXY(I2,2) - CXY(I,2)) / DIS
                      CROSS_DIR = WAVE_DIR_X * CELL_DIR_X + 
     &                            WAVE_DIR_Y * CELL_DIR_Y
                      If(CROSS_DIR .GE. 0.90) Then
                         WB_SLOPE_X = 0.0
                         WB_SLOPE_Y = 0.0
                         WB_SLOPE = 0.0
                         Do J = 1, CELL_POLYGEN(I2)
		                    If(CFM(CELL_SIDE(I2,J,1)) .EQ. 1.0) Then
                              WB_SLOPE_X = WB_SLOPE_X + 
     &                        WIX(I2,J)*(ELF(CELL_SIDE(I2,J,2))-ELF(I2))
                              WB_SLOPE_Y = WB_SLOPE_Y + 
     &                        WIY(I2,J)*(ELF(CELL_SIDE(I2,J,2))-ELF(I2))
                            Endif
                         Enddo
c-----        local slope                      
                         WB_SLOPE = Sqrt(WB_SLOPE_X**2.+WB_SLOPE_Y**2.) 
                         WB = (ELF(I2) - EL(I2)) / DTI
c-----        wave breaking zone
                         If(DIS .GE. 0.0 .AND. DIS .LE. 
     &                      Dmax1(1.5*DD,5.*Sqrt(AREA(I)))) Then
c                         If(WB_SLOPE .GE. Tan(WB_ANGLE_I)) Then 
                            Do K = 1, KBM
                               Do J = 1, CELL_POLYGEN(I2)
		                          If(CFM(CELL_SIDE(I2,J,1)) .EQ. 1.0) 
     &                               AS(I2,K,J) = 0.0
                               Enddo
                               AP(I2,K) = 1.0
	                         AB(I2,K) = 0.0
	                         AT(I2,K) = 0.0
	                         BB(I2,K) = 0.0
C                  VIS_BW(I2,K) = 0.1 * (DC(I2) * (1. + ZZ(K))) ** 2. *
C     &                          Abs(WB) / DC(I2)
C     &             VIS_BW(I2,K) = 0.1*Sqrt(9.8*DD)*DD*(1.-DIS/(0.5*DD))
C                    IF(DIS .LE. 2.*Sqrt(AREA(I))) 
C     &             VIS_BW(I2,K) = 0.5*Sqrt(9.8*DD)*DD* (1. + ZZ(K))
                               C(I2,K) = 1.0
                            Enddo
                            WB_MASK2(I2) = 1.0
                         Endif
                      Endif
C-------------------  BACK
                      If(CROSS_DIR .LE. -0.90) Then
                         WB_SLOPE_X = 0.0
                         WB_SLOPE_Y = 0.0
                         WB_SLOPE = 0.0
                         Do J = 1, CELL_POLYGEN(I2)
		                    If(CFM(CELL_SIDE(I2,J,1)) .EQ. 1.0) Then
                              WB_SLOPE_X = WB_SLOPE_X + 
     &                        WIX(I2,J)*(ELF(CELL_SIDE(I2,J,2))-ELF(I2))
                              WB_SLOPE_Y = WB_SLOPE_Y + 
     &                        WIY(I2,J)*(ELF(CELL_SIDE(I2,J,2))-ELF(I2))
                            Endif
                         Enddo
c-----        local slope                      
                         WB_SLOPE = Sqrt(WB_SLOPE_X**2.+WB_SLOPE_Y**2.) 
                         WB = (ELF(I2) - EL(I2)) / DTI
c-----        wave breaking zone
                         
                         If(DIS .GE. 0.0 .AND. DIS .LE. 
     &                      Dmax1(1.0*DD,5.*Sqrt(AREA(I)))) Then
C                         If(WB_SLOPE .GE. Tan(WB_ANGLE_I)) Then 
                            Do K = 1, KBM
                               Do J = 1, CELL_POLYGEN(I2)
		                          If(CFM(CELL_SIDE(I2,J,1)) .EQ. 1.0) 
     &                               AS(I2,K,J) = 0.0
                               Enddo
                               AP(I2,K) = 1.0
	                         AB(I2,K) = 0.0
	                         AT(I2,K) = 0.0
	                         BB(I2,K) = 0.0
C                  VIS_BW(I2,K) = 0.1 * (DC(I2) * (1. + ZZ(K))) ** 2. *
C     &                          Abs(WB) / DC(I2)
C                   VIS_BW(I2,K) = 0.5 * Sqrt(9.8*DD)*DD* (1. + ZZ(K))

c                               C(I2,K) = 1.0
                            Enddo
                            WB_MASK2(I2) = 1.0
                         Endif
                      Endif
                      
                   Endif
                Enddo 
c                if(scope .ne. 0.0) then
                WB = (ELF(I) - EL(I)) / DTI

                Do K = 1, KBM
                   Do J = 1, CELL_POLYGEN(I)
		              If(CFM(CELL_SIDE(I,J,1)) .EQ. 1.0) 
     &                   AS(I,K,J) = 0.0
                   Enddo
                   AP(I,K) = 1.0
	             AB(I,K) = 0.0
	             AT(I,K) = 0.0
	             BB(I,K) = 0.0
C                   VIS_BW(I,K) = 0.1 * (DC(I) * (1. + ZZ(K))) ** 2. *
C     &                          Abs(WB) / DC(I)
C                   VIS_BW(I,K) = 1.0 * Sqrt(9.8*DD) * DD* (1. + ZZ(K))
                   C(I,K) = 1.0
                Enddo
                WB_MASK2(I) = 1.0
c                endif
               Endif
            Endif  
C--------   define when to stop wave breaking
            If(SLOPE_MAX .GE. Tan(WB_ANGLE_E) .AND.
     &         SLOPE_MAX .LT. Tan(WB_ANGLE_I)) Then
                 Do I2 = 1, IJM
                   DIS = Sqrt((CXY(I2,1) - CXY(I,1)) ** 2. +
     &                        (CXY(I2,2) - CXY(I,2)) ** 2.)
                   If(DIS .LE. SCOPE .AND. I2 .NE. I) Then
                      CELL_DIR_X = (CXY(I2,1) - CXY(I,1)) / DIS
                      CELL_DIR_Y = (CXY(I2,2) - CXY(I,2)) / DIS
                      CROSS_DIR = WAVE_DIR_X * CELL_DIR_X + 
     &                            WAVE_DIR_Y * CELL_DIR_Y
                      If((CROSS_DIR) .GE. 0.90) Then
                         If(WB_MASK(I2) .EQ. 1.0) Then
                       
c-----test                       
                         WB_SLOPE_X = 0.0
                         WB_SLOPE_Y = 0.0
                         WB_SLOPE = 0.0
                         Do J = 1, CELL_POLYGEN(I2)
		                    If(CFM(CELL_SIDE(I2,J,1)) .EQ. 1.0) Then
                              WB_SLOPE_X = WB_SLOPE_X + 
     &                        WIX(I2,J)*(ELF(CELL_SIDE(I2,J,2))-ELF(I2))
                              WB_SLOPE_Y = WB_SLOPE_Y + 
     &                        WIY(I2,J)*(ELF(CELL_SIDE(I2,J,2))-ELF(I2))
                            Endif
                         Enddo
                         Do K = 1, KBM
                            Do J = 1, CELL_POLYGEN(I2)
		                        If(CFM(CELL_SIDE(I2,J,1)) .EQ. 1.0) 
     &                             AS(I2,K,J) = 0.0
                            Enddo
                            AP(I2,K) = 1.0
	                      AB(I2,K) = 0.0
	                      AT(I2,K) = 0.0
	                      BB(I2,K) = 0.0
c                  VIS_BW(I2,K) = 0.1 * (DC(I2) * (1. + ZZ(K))) ** 2. *
c     &                          Abs(WB) / DC(I2)

                            C(I2,K) = 1.0
                         Enddo
                         WB_MASK2(I2) = 1.0
                       Endif
                    Endif
                  Endif   
                Enddo 
            Endif
            If(SLOPE_MAX .LT. Tan(WB_ANGLE_E)) Then
                Do I2 = 1, IJM
                   DIS = Sqrt((CXY(I2,1) - CXY(I,1)) ** 2. +
     &                        (CXY(I2,2) - CXY(I,2)) ** 2.)
                   If(DIS .LE. SCOPE .AND. I2 .NE. I) Then
                      CELL_DIR_X = (CXY(I2,1) - CXY(I,1)) / DIS
                      CELL_DIR_Y = (CXY(I2,2) - CXY(I,2)) / DIS
                      CROSS_DIR = WAVE_DIR_X * CELL_DIR_X + 
     &                            WAVE_DIR_Y * CELL_DIR_Y
                      If(Abs(CROSS_DIR) .GE. 0.90) Then
                         WB_MASK2(I2) = 0.0
                      Endif
                   Endif
                Enddo 
            Endif
1001        Continue 
10001       Continue            
           Endif
      Enddo
      Do I = 1, IJM
         WB_MASK(I) = WB_MASK2(I)
      Enddo          
      Return
      End
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
