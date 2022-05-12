C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (11/01/2010)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c      computes the astrotide constants: the amplitude factor and the       #
c      initial phase angles                                                 #
c                                                                           #
c############################################################################

      Subroutine  GETCOR(KY,KKM,KDAY,AMP0,PHASE0,FREQ,NUMT)
	Double precision PHASE0(NUMT),FREQ(NUMT),PHASE(NUMT),AMP0(NUMT)
      Double precision A(11),B(11),FF(11),UU(11),PI,tharm(11)
      Double precision S,HH,P,XN,PP,TAO,C,D
      Integer KC(11,6) ,MONTH(12)
C===========================================================================c
c	    the tidal constitite M2 S2 N2 K2 K1 O1 P1 Q1 M4 MS4 M6            c
C===========================================================================c
      Data MONTH/-1,30,58,89,119,150,180,211,242,272,303,333/
      Data ((KC(I,J),J=1,6),I=1,11)
     $       /2,0,0,0,0,0,2,2,-2,0,0,0,2,-1,0,1,0,0,2,2,0,0,0,0,
     $        1,1,0,0,0,1,1,-1,0,0,0,-1,1,1,-2,0,0,-1,1,-2,0,1,0,
     $        -1,4,0,0,0,0,0,4,2,-2,0,0,0,6,0,0,0,0,0/
      PI = 2. * ACos(0.)
      NADDI = MONTH(KKM) + KDAY + Int((KY-1901)/4)
      IF(Mod(KY,4) .EQ. 0 .AND. KKM .GT. 2) NADDI = NADDI + 1

      S = 277.02 + 129.3848 * (KY - 1900) + 13.1764 * (NADDI)
      HH = 280.19 - 0.2387 * (KY - 1900) + 0.9857 * (NADDI)
      P = 334.39 + 40.6625 * (KY - 1900) + 0.1114 * (NADDI)
      XN = 100.84 + 19.3282 * (KY - 1900) + 0.053 * (NADDI)
      PP = 281.22 + 0.0172 * (KY - 1900) + 0.00005 * 0.053 * (NADDI)
      TAO = -S + HH
C===========================================================================c
C	Attention in all coefficient of "XN" equal to 0                       c
C===========================================================================c
      Do I = 1, 11
        PHASE(I) = KC(I,1) * TAO + KC(I,2) * S + KC(I,3) * HH + 
     &	         KC(I,4) * P + KC(I,5) * PP + KC(I,6) * 90.
        FREQ(I) = 14.49205211 * KC(I,1) + 0.54901653 * KC(I,2) +
     &           0.04106864 * KC(I,3) + 0.00646183 * KC(I,4) +  
     &           0.00000196 * KC(I,5)
	Enddo
      C = P * PI / 180.0
      D = XN * PI / 180.0
c---------------------------------------------------------------------------c
      A(1) = 0.0005 * Cos(-2*D) - 0.0373 * Cos(-D) + 0.0006 * Cos(2*C) +
     &      0.0002 * Cos(2*C+D) + Cos(0*D)
      A(4) = -0.0128 * Cos(-D) + 0.2980 * Cos(D) + 0.0324*Cos(2*D) +
     &	   Cos(0*D)
      A(5)= 0.0002 * Cos(-2*C-D) + 0.0001 * Cos(-2*D) - 0.0198 * Cos(-D)
     &      + 0.1356 * Cos(D) - 0.0029 * Cos(2*D) + Cos(0*D)
      A(6) =- 0.0058 * Cos(-2*D) + 0.1885 * Cos(-D) + 0.0002 * 
     &   	  Cos(2*C-D) - 0.0064 * Cos(2*C) - 0.0010 * Cos(2*C+D) + 
     &      Cos(0*D)     
      A(7) = 0.0008 * Cos(-2*D) - 0.0112 * Cos(-D) - 0.0015 * Cos(2*C) -
     &       0.0003 * Cos(2*C+D) + Cos(0*D)
C----------------------------------------------------------------------------c
      B(1) = -0.0005 * Sin(-2*D) - 0.0373 * Sin(-D) + 0.0006 * Sin(2*C)
     &      + 0.0002 * Sin(2*C+D) + Sin(0*D)
      B(4) = -0.0128 * Sin(-D) + 0.2980 * Sin(D) + 0.0324 * Sin(2*D) +
     &	   Sin(0*D)
      B(5) = 0.0002 * Sin(-2*C-D) + 0.0001 * Sin(-2*D) - 0.0198 * 
     &	   Sin(-D) + 0.1356 * Sin(D) - 0.0029 * Sin(2*D) + Sin(0*D)
      B(6) = -0.0058 * Sin(-2*D) + 0.1885 * Sin(-D) + 0.0002 * 
     &	   Sin(2*C-D) - 0.0064 * Sin(2*C) - 0.0010 * Sin(2*C+D) + 
     &       Sin(0*D)     
      B(7) = 0.0008 * Sin(-2*D) - 0.0112 * Sin(-D) - 0.0015 * Sin(2*C) -
     &       0.0003 * Sin(2*C+D) + Sin(0*D)
c----------------------------------------------------------------------------c
      FF(1) = Sqrt(A(1) **2 + B(1) **2)
      FF(2) = 1.
      FF(3) = FF(1)
      FF(4) = Sqrt(A(4) ** 2 + B(4) **2)
      FF(5) = Sqrt(A(5) ** 2 + B(5) **2)
      FF(6) = Sqrt(A(6) ** 2 + B(6) **2)
      FF(7) = Sqrt(A(7) ** 2 + B(7) **2)
      FF(8) = FF(6)
      FF(9) = FF(1) ** 2
      FF(10) = FF(1)
      FF(11) = FF(1) ** 3
      UU(1) = Atan(B(1) / A(1))
      UU(2) = 0.
      UU(3) = UU(1)
      UU(4) = Atan(B(4) / A(4))
      UU(5) = Atan(B(5) / A(5))
      UU(6) = Atan(B(6) / A(6))
      UU(7) = Atan(B(7) / A(7))
      UU(8) = UU(6)
      UU(9) = UU(1) * 2
      UU(10) = UU(1)
      UU(11) = UU(1) * 3
      Do I11 = 1, 11
         THARM(i11) = FREQ(i11)
      Enddo
	  THARM(3) = FREQ(5)
	  THARM(4) = FREQ(6)
	  THARM(5) = FREQ(3)
	  THARM(6) = FREQ(4)
	Do I11 = 1, 11
	   FREQ(I11) = THARM(I11)
	Enddo
      Do I11 = 1, 11
         PHASE0(I11) = PHASE(I11) * PI / 180.0 + UU(I11)
      Enddo
	PHASE0(3) = PHASE(5) * PI / 180.0 + UU(5)
	PHASE0(4) = PHASE(6) * PI / 180.0 + UU(6)
	PHASE0(5) = PHASE(3) * PI / 180.0 + UU(3)
	PHASE0(6) = PHASE(4) * PI / 180.0 + UU(4)
      Do I11 = 1, 11
         AMP0(I11) = FF(I11)
      Enddo
	AMP0(3) = FF(5)
	AMP0(4) = FF(6)
	AMP0(5) = FF(3)
	AMP0(6) = FF(4)
	Return
      End 
