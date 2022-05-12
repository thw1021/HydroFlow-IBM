C==================================================================================C
C                         RADOM SIGNALS GENERATING                                 C
C----------------------------------------------------------------------------------C
      Subroutine TurGenerator(TIME1,NSTEPS)
 
      Include './Include/OCERM_INF'
      Include './Include/VORGEN_INF'

	Parameter (J1 = 10000, K1 =100)
	
      Double precision TIME1, WIDTH
      Integer NSTEPS

C-----   TEMPARAL ARRAY TO STORE DATAS

      Dimension U1(J1,K1),V1(J1,K1),W1(J1,K1),
     &            U2(J1,K1),V2(J1,K1),W2(J1,K1),
     &            U2S(J1,K1),V2S(J1,K1),W2S(J1,K1),
     &            U3(J1,K1),V3(J1,K1),W3(J1,K1),
     &            U4(J1,K1),V4(J1,K1),W4(J1,K1)
     
C-----   RADOM NUMBER GENERATING FUNCTIONS
      
      Integer NP1, NP2, NP3, NP4
      Real NRND1, GFUN, HFUN
      WIDTH = ZR - ZL
C-----  Structure convective velocity 

      P1_T = P1_LX / C1
      P2_T = P2_LX / C2
      P2S_T = P2S_LX / C2S
      P3_T = P3_LX / C3
      P4_T = P4_LX / C4
      
c-----  Number of structures per mode N(p)

      NP1 = Int(WIDTH * (P1_UP-P1_LOW) / (4.*P1_LY*P1_LZ)) * IPLANE
      NP2 = Int(WIDTH * (P2_UP-P2_LOW) / (4.*P2_LY*P2_LZ)) * IPLANE
      NP2S = Int(WIDTH * (P2S_UP-P2S_LOW) / (4.*P2S_LY*P2S_LZ)) * IPLANE
      NP3 = Int(WIDTH * (P3_UP-P3_LOW) / (4.*P3_LY*P3_LZ)) * IPLANE
      NP4 = Int(WIDTH * (P4_UP-P4_LOW) / (4.*P4_LY*P4_LZ)) * IPLANE

C-------------------------------------------------------------------------C
C     PICKUP RANDOMLY THE SPATIAL AND TEMPARAL PARAMETERS OF STRUCTURES   C
C-------------------------------------------------------------------------C 
C----  P1 MODE STEL-JHSK
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(K,T_NOR)
!$OMP DO      
      Do K = 1, NP1
     
         !Randomly pick locations and instantaneous time
         
         T_NOR = (TIME1 - SP1_DEFINE(K,1) - P1_T) / P1_T
         If(T_NOR .GE. 1.0) Then
         
         !Randomly select  a new structure 
          
C            SP1_DEFINE(K,1) = NRND1(RP1(1)) * (P1_T) + TIME1
            SP1_DEFINE(K,1) = -NRND1(RT1(K)) * (P1_T) + TIME1
            SP1_DEFINE(K,2) = NRND1(RP1(2)) * (P1_UP-P1_LOW) + P1_LOW
            SP1_DEFINE(K,3) = NRND1(RP1(3)) * WIDTH + ZL

           ! Random sign of u
           
           SR1(K,1) = NRAB1(SR1S(K,1),1,NSTEPS)
           If(Mod(SR1(K,1),2.) .EQ. 0) Then
              SR1(K,1) = -1
           Else
              SR1(K,1) = 1
           Endif  
     
           ! Random sign of v
           
           SR1(K,2) = NRAB1(SR1S(K,2),1,NSTEPS)
           If(Mod(SR1(K,2),2.) .EQ. 0) Then
              SR1(K,2) = -1
           Else
              SR1(K,2) = 1
           Endif 
           
            ! Random sign of w
            
           SR1(K,3) = NRAB1(SR1S(K,3),1,NSTEPS)
           If(Mod(SR1(K,3),2.) .EQ. 0) Then
              SR1(K,3) = -1
           Else
              SR1(K,3) = 1
           Endif 

         Endif
      Enddo 
!$OMP END DO     
C----  P2 MODE LEGS-JHSK
!$OMP DO      
      Do K = 1, NP2
      
         ! Randomly pick locations and instantaneous time
         
         T_NOR = (TIME1 - SP2_DEFINE(K,1) - P2_T) / P2_T
         If(T_NOR .GE. 1.0) Then
         
         ! Randomly select a new structure  
         
C           SP2_DEFINE(K,1) = NRND1(RP2(1)) * (P2_T) + TIME1
           SP2_DEFINE(K,1) = -NRND1(RT2(K)) * (P2_T) + TIME1
           SP2_DEFINE(K,2) = NRND1(RP2(2)) * (P2_UP-P2_LOW) + P2_LOW
           SP2_DEFINE(K,3) = NRND1(RP2(3)) * WIDTH + ZL

           ! Random sign of u
           
           SR2(K,1) = NRAB1(SR2S(K,1),1,NSTEPS)
           If(Mod(SR2(K,1),2.) .EQ. 0) Then
              SR2(K,1) = -1
           Else
              SR2(K,1) = 1
           Endif  
     
           ! Random sign of v
           
           SR2(K,2) = NRAB1(SR2S(K,2),1,NSTEPS)
           If(Mod(SR2(K,2),2.) .EQ. 0) Then
              SR2(K,2) = -1
           Else
              SR2(K,2) = 1
           Endif 
           
            ! Random sign of w
            
           SR2(K,3) = NRAB1(SR2S(K,3),1,NSTEPS)
           If(Mod(SR2(K,3),2.) .EQ. 0) Then
              SR2(K,3) = -1
           Else
              SR2(K,3) = 1
           Endif 
           
          Endif 

      Enddo 
!$OMP END DO            
C----  P2S MODE HEAD-TRANS
!$OMP DO
      Do K = 1, NP2S
      
         !Randomly pick locations and instantaneous time
         
         T_NOR = (TIME1 - SP2S_DEFINE(K,1) - P2S_T) / P2S_T
         If(T_NOR .GE. 1.0) Then
         
           ! Randomly select structure as a new one 
           
C           SP2S_DEFINE(K,1) = NRND1(RP2S(1)) * (P2S_T) + TIME1
           SP2S_DEFINE(K,1) = -NRND1(RT2S(K)) * (P2S_T) + TIME1
           SP2S_DEFINE(K,2) = NRND1(RP2S(2)) * (P2S_UP-P2S_LOW)+P2S_LOW
           SP2S_DEFINE(K,3) = NRND1(RP2S(3)) * WIDTH + ZL

           ! Random sign of u
           
           SR22(K,1) = NRAB1(SR2SS(K,1),1,NSTEPS)
           If(Mod(SR22(K,1),2.) .EQ. 0) Then
              SR22(K,1) = -1
           Else
              SR22(K,1) = 1
           Endif  
     
           ! Random sign of v
           
           SR22(K,2) = NRAB1(SR2SS(K,2),1,NSTEPS)
           If(Mod(SR22(K,2),2.) .EQ. 0) Then
              SR22(K,2) = -1
           Else
              SR22(K,2) = 1
           Endif 
           
            ! Random sign of w
            
           SR22(K,3) = NRAB1(SR2SS(K,3),1,NSTEPS)
           If(Mod(SR22(K,3),2.) .EQ. 0) Then
              SR22(K,3) = -1
           Else
              SR22(K,3) = 1
           Endif 
           
         Endif  

      Enddo 
!$OMP END DO      
C----  P3 MODE WAKE-GAUSS
!$OMP DO
      Do K = 1, NP3
         ! Randomly pick locations and instantaneous time
         
         T_NOR = (TIME1 - SP3_DEFINE(K,1) - P3_T) / P3_T
         
         If(T_NOR .GE. 1.0) Then
          
           ! Randomly select structure as a new one 
           
C           SP3_DEFINE(K,1) = NRND1(RT3(K)) * (P3_T) + T
C           SP3_DEFINE(K,2) = NRND1(RY3(K)) * (P3_UP-P3_LOW) + P3_LOW
C           SP3_DEFINE(K,3) = NRND1(RZ3(K)) * WIDTH + ZL
C           SP3_DEFINE(K,1) = NRND1(RP3(1)) * (P3_T) + TIME1
           SP3_DEFINE(K,1) = -NRND1(RT3(K)) * (P3_T) + TIME1
           SP3_DEFINE(K,2) = NRND1(RP3(2)) * (P3_UP-P3_LOW) + P3_LOW
           SP3_DEFINE(K,3) = NRND1(RP3(3)) * WIDTH + ZL
            
           ! Random sign of u
           
           SR3(K,1) = NRAB1(SR3S(K,1),1,NSTEPS)
           If(Mod(SR3(K,1),2.) .EQ. 0) Then
              SR3(K,1) = -1
           Else
              SR3(K,1) = 1
           Endif  
     
           ! Random sign of v
           
           SR3(K,2) = NRAB1(SR3S(K,2),1,NSTEPS)
           If(Mod(SR3(K,2),2.) .EQ. 0) Then
              SR3(K,2) = -1
           Else
              SR3(K,2) = 1
           Endif 
           
            ! Random sign of w
            
           SR3(K,3) = NRAB1(SR3S(K,3),1,NSTEPS)
           If(Mod(SR3(K,3),2.) .EQ. 0) Then
              SR3(K,3) = -1
           Else
              SR3(K,3) = 1
           Endif 
        Endif
      Enddo 
!$OMP END DO      
C----  P4 MODE WAKE-GAUSS
!$OMP DO
      Do K = 1, NP4
         ! Randomly pick locations and instantaneous time
         
         T_NOR = (TIME1 - SP4_DEFINE(K,1) - P4_T) / P4_T
         
         If(T_NOR .GE. 1.0) Then
         
           ! Randomly select structure as a new one 
C           SP4_DEFINE(K,1) = NRND1(RP4(1)) * (P4_T) + TIME1
           SP4_DEFINE(K,1) = -NRND1(RT4(K)) * (P4_T) + TIME1
           SP4_DEFINE(K,2) = NRND1(RP4(2)) * (P4_UP-P4_LOW) + P4_LOW
           SP4_DEFINE(K,3) = NRND1(RP4(3)) * WIDTH + ZL
           
           ! Random sign of u
           
           SR4(K,1) = NRAB1(SR4S(K,1),1,NSTEPS)
           If(Mod(SR4(K,1),2.) .EQ. 0) Then
              SR4(K,1) = -1
           Else
              SR4(K,1) = 1
           Endif 

           ! Random sign of v
           
           SR4(K,2) = NRAB1(SR4S(K,2),1,NSTEPS)
           If(Mod(SR4(K,2),2.) .EQ. 0) Then
              SR4(K,2) = -1
           Else
              SR4(K,2) = 1
           Endif 

           ! Random sign of w
           
           SR4(K,3) = NRAB1(SR4S(K,3),1,NSTEPS)
           If(Mod(SR4(K,3),2.) .EQ. 0) Then
              SR4(K,3) = -1
           Else
              SR4(K,3) = 1
           Endif  
           
        Endif  
      Enddo 
!$OMP END DO
!$OMP END PARALLEL

C---------------------   END TASK  ---------------------------------------C
C-------------------------------------------------------------------------C
C         GENERATING RANDOM SIGNALS                                       C
C-------------------------------------------------------------------------C
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(II,J11,J22,I1S,I1E,I2S,I2E,ID,
!$OMP&                        I22S,I22E,I3S,I3E,I4S,I4E,I,J,K,DW)
!$OMP DO
      Do II = 1, IPLANE

         J11 = (II - 1) * NUM_CELL / IPLANE + 1
         J22 = II * NUM_CELL / IPLANE

         I1S = (II - 1) * INT( 1. * NP1 / IPLANE) + 1
         I1E = II * INT(1. * NP1 / IPLANE)
         I2S = (II - 1) * INT( 1. * NP2 / IPLANE) + 1
         I2E = II * INT(1. * NP2 / IPLANE)
         I22S = (II - 1) * INT( 1. * NP2S / IPLANE) + 1
         I22E = II * INT(1. * NP2S / IPLANE)
         I3S = (II - 1) * INT( 1. * NP3 / IPLANE) + 1
         I3E = II * INT(1. * NP3 / IPLANE)
         I4S = (II - 1) * INT( 1. * NP4 / IPLANE) + 1
         I4E = II * INT(1. * NP4 / IPLANE)
        
         Do J = J11, J22
         ID = ID_CELL(J)

         Do K = 1, NUM_VER
       
         ! DETTING ZERO VALUES FOR SUMMING
         
         UDIS(J,K) = 0.0
         VDIS(J,K) = 0.0
         WDIS(J,K) = 0.0
         
         U1(J,K) = 0.0
         V1(J,K) = 0.0
         W1(J,K) = 0.0
         
         U2(J,K) = 0.0
         V2(J,K) = 0.0
         W2(J,K) = 0.0

         U2S(J,K) = 0.0
         V2S(J,K) = 0.0
         W2S(J,K) = 0.0
         
         U3(J,K) = 0.0
         V3(J,K) = 0.0
         W3(J,K) = 0.0
         
         U4(J,K) = 0.0
         V4(J,K) = 0.0
         W4(J,K) = 0.0
         
C------  P1 MODE STEL-JHSK  
         
         Do I = I1S, I1E
         
            T_NOR = (TIME1-SP1_DEFINE(I,1)-P1_T) / P1_T
         
            U1(J,K) = U1(J,K) + SR1(I,1) * GFUN(T_NOR) * 
     &                GFUN((YCEN(J,K)-SP1_DEFINE(I,2))/P1_LY) * 
     &                HFUN((ZCEN(J,K)-SP1_DEFINE(I,3))/P1_LZ)
     
            V1(J,K) = V1(J,K) - SR1(I,2) * GFUN(T_NOR) * 
     &                GFUN((YCEN(J,K)-SP1_DEFINE(I,2))/P1_LY) * 
     &                HFUN((ZCEN(J,K)-SP1_DEFINE(I,3))/P1_LZ)
     
            W1(J,K) = W1(J,K) + SR1(I,3) * GFUN(T_NOR) * 
     &                HFUN((YCEN(J,K)-SP1_DEFINE(I,2))/P1_LY) * 
     &                GFUN((ZCEN(J,K)-SP1_DEFINE(I,3))/P1_LZ)
     
         Enddo 
         
C------  P2 MODE LEGS-JHSK

         Do I = I2S, I2E
         
            T_NOR = (TIME1-SP2_DEFINE(I,1)-P2_T)/ P2_T
          
            U2(J,K) = U2(J,K) + SR2(I,1) * GFUN(T_NOR) * 
     &                GFUN((YCEN(J,K)-SP2_DEFINE(I,2))/P2_LY) * 
     &                HFUN((ZCEN(J,K)-SP2_DEFINE(I,3))/P2_LZ)
     
            V2(J,K) = V2(J,K) - SR2(I,2) * GFUN(T_NOR) * 
     &                GFUN((YCEN(J,K)-SP2_DEFINE(I,2))/P2_LY) * 
     &                HFUN((ZCEN(J,K)-SP2_DEFINE(I,3))/P2_LZ)
     
            W2(J,K) = W2(J,K) + SR2(I,3) * GFUN(T_NOR) * 
     &                HFUN((YCEN(J,K)-SP2_DEFINE(I,2))/P2_LY) * 
     &                GFUN((ZCEN(J,K)-SP2_DEFINE(I,3))/P2_LZ)

         Enddo 
           
C------  P2S MODE HEAD-TRANS

         Do I = I22S, I22E
         
            T_NOR = (TIME1-SP2S_DEFINE(I,1)-P2S_T)/P2S_T
          
            U2S(J,K) = U2S(J,K) - SR22(I,1) * GFUN(T_NOR) * 
     &                HFUN((YCEN(J,K)-SP2S_DEFINE(I,2))/P2S_LY) * 
     &                GFUN((ZCEN(J,K)-SP2S_DEFINE(I,3))/P2S_LZ)
    
            V2S(J,K) = V2S(J,K) + SR22(I,2) * T_NOR*GFUN(T_NOR) * 
     &                GFUN((YCEN(J,K)-SP2S_DEFINE(I,2))/P2S_LY) * 
     &                GFUN((ZCEN(J,K)-SP2S_DEFINE(I,3))/P2S_LZ)
     
            W2S(J,K) = W2S(J,K) + SR22(I,3) * GFUN(T_NOR) * 
     &                GFUN((YCEN(J,K)-SP2S_DEFINE(I,2))/P2S_LY) * 
     &                GFUN((ZCEN(J,K)-SP2S_DEFINE(I,3))/P2S_LZ) *
     &                (ZCEN(J,K)-SP2S_DEFINE(I,3))/P2S_LZ

         Enddo 
 
C------  P3 MODE WAKE-GAUSS

         Do I = I3S, I3E
         
            T_NOR = (TIME1 - SP3_DEFINE(I,1) - P3_T) / P3_T
            
            U3(J,K) = U3(J,K) + SR3(I,1) * GFUN(T_NOR) * 
     &                GFUN((YCEN(J,K)-SP3_DEFINE(I,2))/P3_LY) * 
     &                GFUN((ZCEN(J,K)-SP3_DEFINE(I,3))/P3_LZ)
     
            V3(J,K) = V3(J,K) + SR3(I,2) * GFUN(T_NOR) * 
     &                GFUN((YCEN(J,K)-SP3_DEFINE(I,2))/P3_LY) * 
     &                GFUN((ZCEN(J,K)-SP3_DEFINE(I,3))/P3_LZ)
     
            W3(J,K) = W3(J,K) + SR3(I,3) * GFUN(T_NOR) * 
     &                GFUN((YCEN(J,K)-SP3_DEFINE(I,2))/P3_LY) * 
     &                GFUN((ZCEN(J,K)-SP3_DEFINE(I,3))/P3_LZ) 
            
         Enddo 

C------  P4 MODE WAKE-GAUSS
        
         Do I = I4S, I4E
         
            T_NOR = (TIME1 - SP4_DEFINE(I,1) - P4_T) / P4_T
         
            U4(J,K) = U4(J,K) + SR4(I,1) * GFUN(T_NOR) * 
     &                GFUN((YCEN(J,K)-SP4_DEFINE(I,2))/P4_LY) * 
     &                GFUN((ZCEN(J,K)-SP4_DEFINE(I,3))/P4_LZ)
     
            V4(J,K) = V4(J,K) + SR4(I,2) * GFUN(T_NOR) * 
     &                GFUN((YCEN(J,K)-SP4_DEFINE(I,2))/P4_LY) * 
     &                GFUN((ZCEN(J,K)-SP4_DEFINE(I,3))/P4_LZ)
     
            W4(J,K) = W4(J,K) + SR4(I,3) * GFUN(T_NOR) * 
     &                GFUN((YCEN(J,K)-SP4_DEFINE(I,2))/P4_LY) * 
     &                GFUN((ZCEN(J,K)-SP4_DEFINE(I,3))/P4_LZ) 
     
         Enddo
        
C------  SUMMING ALL OF THE RANDOM SIGNALS  
c         DW = DC(ID) * (1. + ZZ(K))
c         If(DW .GE. DWINTERFACE) Then

            U_FLU(J,K) = U1(J,K) + U2(J,K) +
     &                   U3(J,K) + U4(J,K) +
     &                   U2S(J,K)
     
            V_FLU(J,K) = V1(J,K) + V2(J,K) +
     &                   V3(J,K) + V4(J,K) +
     &                   V2S(J,K) 
 
            W_FLU(J,K) = W1(J,K) + W2(J,K) +
     &                   W3(J,K) + W4(J,K) +
     &                   W2S(J,K)
c         Endif

c         U_FLU(J,K) = U1(J,K)/Sqrt(NP1*1.) + U2(J,K)/Sqrt(NP2*1.)+
c     &               U3(J,K)/Sqrt(NP3*1.) + U4(J,K)/Sqrt(NP4*1.) +
c     &               U2S(J,K)/Sqrt(NP2S*1.)
      
c         V_FLU(J,K) = V1(J,K)/Sqrt(NP1*1.) + V2(J,K)/Sqrt(NP2*1.)+
c     &               V3(J,K)/Sqrt(NP3*1.) + V4(J,K)/Sqrt(NP4*1.) +
c     &               V2S(J,K)/Sqrt(NP2S*1.) 
  
c         W_FLU(J,K) = W1(J,K)/Sqrt(NP1*1.) + W2(J,K)/Sqrt(NP2*1.)+
c     &               W3(J,K)/Sqrt(NP3*1.) + W4(J,K)/Sqrt(NP4*1.) +
c     &               W2S(J,K)/Sqrt(NP2S*1.) 
        
        Enddo
        Enddo
       
      Enddo
!$OMP END DO
!$OMP END PARALLEL 
      Return
      End
      
C----  Defining functions

      Function GFUN(T)
      Double precision T
          GFUN = Exp(-T **2. / (2./2./2.)) * 1.5
      Return    
      End   
         
      Function HFUN(T)
      Double precision T
          HFUN = (1. - Cos(2. * Acos(-1.0) * T)) / 
     &           (2. * Acos(-1.0) * T * Sqrt(0.214)+1.E-6) 
      Return
      End  
      
C----  Generating Random value ----------------------------------------c

	Real Function NRND1(R)
	Double precision S,U,V,R
	S=65536.0
	U=2053.0
	V=13849.0
	M=R/S
	R=R-M*S
	R=U*R+V
	M=R/S
	R=R-M*S
	NRND1=R/S
	Return
	End
	
C----   Generating Random Sign  ----------------------------------------c	

	Function NRAB1(R,A,B)
	Double precision R,S
	Integer A,B
	S=B-A+1.0
	K=Log(S-0.5)/Log(2.0)+1
	L=1
	Do 10 I=1,K
10	L=2*L
	K=1
	S=4.0*L
20	If (K.LE.1) Then
	  R=R+R+R+R+R
	  M=R/S
	  R=R-M*S
	  J=A+R/4.0
	  If (J.LE.B) Then
	    NRAB1=J
	    K=K+1
	  Endif
	  Goto 20
	Endif
	Return
	End
      
      
