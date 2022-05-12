c====================================================================================================c
c            Inilizating the parameters                                                              c
c====================================================================================================c      
      Subroutine INIVOR
      Include './Include/OCERM_INF'
      Include './Include/VORGEN_INF'
     
C----  TEMP FILES TO SAVE THE INSTANTANEOUS SIGNALS
     
      Dimension A11(200,2),A22(200,2),A33(200,2),
     &           A21(200,2),A31(200,2),A32(200,2) 
     
      Double precision  UNIT_LENGTH, TIME1 
      Real NRND1

C-----  Read Pre-Reynolds Stress from data sources  ---------------------c

      Open(1000,file='./Include/Reynolds Stress.dat',Status='Unknown') 
      
      Read(1000,*) I11
      Do I = 1, I11
          Read(1000,*) A11(I,1),A11(I,2)
C          A11(I,1) = A11(I,1) * BHPLUX / 2000. 
      Enddo

      Read(1000,*) I22
      Do I = 1, I22
          Read(1000,*) A22(I,1),A22(I,2)
C          A22(I,1) = A22(I,1) * BHPLUX / 2000.
      Enddo

      Read(1000,*) I33
      Do I = 1, I33
          Read(1000,*) A33(I,1),A33(I,2)
C          A33(I,1) = A33(I,1) * BHPLUX / 2000.
      Enddo

      Read(1000,*) I21
      Do I = 1, I21
          Read(1000,*) A21(I,1),A21(I,2)
C          A21(I,1) = A21(I,1) * BHPLUX / 2000.
      Enddo

      Read(1000,*) I31
      Do I = 1, I31
          Read(1000,*) A31(I,1),A31(I,2)
C          A31(I,1) = A31(I,1) * BHPLUX / 2000.
      Enddo

      Read(1000,*) I32
      Do I = 1, I32
          Read(1000,*) A32(I,1),A32(I,2)
C          A32(I,1) = A32(I,1) * BHPLUX / 2000.
      Enddo
      Close(1000)
C-----  Initializing values ---------------------------------------------C
      Open(1000,file='./Include/VorXYZ.dat',Status='Unknown')
c      do i = 1, 1440
c         write(1000, '(2i8,200f12.4)') i, 3, (hc(i)*(1.+zz(k)), cxy(i,2), K = 1, kbm)   
c      enddo
c      stop
      Do I = 1, NUM_CELL
         Read(1000,*) ID_CELL(I), ID_CELL_EDGE(I),
     &                (YCEN(I,K), ZCEN(I,K), K = 1, NUM_VER)
      Enddo
      Close(1000)
C-----  Cholesky decomposition of Reynolds stress  ---------------------c
      Do J = 1, NUM_CELL
      Do K = 1, NUM_VER
C-----  Interplate the spatial AIJ from the a priori database 
         
         YPLUS = YCEN(J,K) * VEL_STAR / 1.E-6
         
         ! A11
         Do I = 1, I11 - 1
            If(YPLUS .GE. A11(I,1) .AND. YPLUS .LT. A11(I+1,1)) Then
                R11(J,K) = (A11(I,2) * (A11(I+1,1) - YPLUS) + 
     &                     A11(I+1,2) * (YPLUS - A11(I,1))) /
     &                    (A11(I+1,1) - A11(I,1)) 
            Endif
         Enddo
         If(YPLUS .GE. A11(I11,1)) Then
            R11(J,K) = A11(I11,2)
         Endif
         R11(J,K) = R11(J,K) * VEL_STAR
         
         ! A22
         Do I = 1, I22 - 1
            If(YPLUS .GE. A22(I,1) .AND. YPLUS .LT. A22(I+1,1)) Then
               R22(J,K) = (A22(I,2) * (A22(I+1,1) - YPLUS) + 
     &                     A22(I+1,2) * (YPLUS - A22(I,1))) /
     &                    (A22(I+1,1) - A22(I,1)) 
            Endif
         Enddo
         If(YPLUS .GE. A22(I22,1)) Then
            R22(J,K) = A22(I22,2)
         Endif
         R22(J,K) = R22(J,K) * VEL_STAR
         
         ! A33
         Do I = 1, I33 - 1
            If(YPLUS .GE. A33(I,1) .AND. YPLUS .LT. A33(I+1,1)) Then
               R33(J,K) = (A33(I,2) * (A33(I+1,1) - YPLUS) + 
     &                     A33(I+1,2) * (YPLUS - A33(I,1))) /
     &                    (A33(I+1,1) - A33(I,1)) 
            Endif
         Enddo
         If(YPLUS .GE. A33(I33,1)) Then
            R33(J,K) = A33(I33,2)
         Endif
         R33(J,K) = R33(J,K) * VEL_STAR
         
         ! A21
         Do I = 1, I21 - 1
            If(YPLUS .GE. A21(I,1) .AND. YPLUS .LT. A21(I+1,1)) Then
               R21(J,K) = (A21(I,2) * (A21(I+1,1) - YPLUS) + 
     &                     A21(I+1,2) * (YPLUS - A21(I,1))) /
     &                    (A21(I+1,1) - A21(I,1)) 
            Endif
         Enddo
         If(YPLUS .GE. A21(I21,1)) Then
            R21(J,K) = A21(I21,2)
         Endif
         R21(J,K) = R21(J,K) * VEL_STAR
         
         ! A31
         Do I = 1, I31 - 1
            If(YPLUS .GE. A31(I,1) .AND. YPLUS .LT. A31(I+1,1)) Then
               R31(J,K) = (A31(I,2) * (A31(I+1,1) - YPLUS) + 
     &                     A31(I+1,2) * (YPLUS - A31(I,1))) /
     &                    (A31(I+1,1) - A31(I,1)) 
            Endif
         Enddo
         If(YPLUS .GE. A31(I31,1)) Then
            R31(J,K) = A31(I31,2)
         Endif
         R31(J,K) = R31(J,K) * VEL_STAR
         
         ! A31
         Do I = 1, I32 - 1
            If(YPLUS .GE. A32(I,1) .AND. YPLUS .LT. A32(I+1,1)) Then
               R32(J,K) = (A32(I,2) * (A32(I+1,1) - YPLUS) + 
     &                     A32(I+1,2) * (YPLUS - A32(I,1))) /
     &                    (A32(I+1,1) - A32(I,1)) 
            Endif
         Enddo
         If(YPLUS .GE. A32(I32,1)) Then
            R32(J,K) = A32(I31,2)
         Endif
         R32(J,K) = R32(J,K) * VEL_STAR
         
      Enddo
      Enddo
 
c---  INFORMATION OF TURBULENCE STRUCTURES  -----------------------------------------------c
      Do K = 1, 10000
         SP1_DEFINE(K,1) = -999.0
         SP2_DEFINE(K,1) = -999.0
         SP2S_DEFINE(K,1) = -999.0
         SP3_DEFINE(K,1) = -999.0
         SP4_DEFINE(K,1) = -999.0
      Enddo       
        
c-----  Characteristic spatial and temporal scales

      UNIT_LENGTH = 1.0 * 1.E-6 / VEL_STAR
      
C-----  Characteristic parameters 

      P1_LX =  P1_X_PLUX * UNIT_LENGTH
      P1_LY =  P1_Y_PLUX * UNIT_LENGTH
      P1_LZ =  P1_Z_PLUX * UNIT_LENGTH
      C1 = C1PLUX * VEL_STAR
      
      P2_LX =  P2_X_PLUX * UNIT_LENGTH
      P2_LY =  P2_Y_PLUX * UNIT_LENGTH
      P2_LZ =  P2_Z_PLUX * UNIT_LENGTH
      C2 = C2PLUX * VEL_STAR

      P2S_LX =  P2S_X_PLUX * UNIT_LENGTH
      P2S_LY =  P2S_Y_PLUX * UNIT_LENGTH
      P2S_LZ =  P2S_Z_PLUX * UNIT_LENGTH
      C2S = C2SPLUX * VEL_STAR

      P3_LX =  P3_X_PLUX * BH
      P3_LY =  P3_Y_PLUX * BH
      P3_LZ =  P3_Z_PLUX * BH
      C3 = C3PLUX * UNIT_VEL
      
      P4_LX =  P4_X_PLUX * BH
      P4_LY =  P4_Y_PLUX * BH
      P4_LZ =  P4_Z_PLUX * BH
      C4 = C4PLUX * UNIT_VEL
      
      P1_LOW = P1_LOW_PLUX * UNIT_LENGTH
      P1_UP  = P1_UP_PLUX  * UNIT_LENGTH
      P2_LOW = P2_LOW_PLUX * UNIT_LENGTH
      P2_UP  = P2_UP_PLUX  * BH
      P2S_LOW = P2S_LOW_PLUX * BH
c      P2S_LOW = P2S_LOW_PLUX * UNIT_LENGTH
      P2S_UP  = P2S_UP_PLUX  * BH
      P3_LOW = P3_LOW_PLUX * BH
      P3_UP  = P3_UP_PLUX  * BH
      P4_LOW = P4_LOW_PLUX * BH
      P4_UP  = P4_UP_PLUX  * BH

C========================================================================C
C----     BEGINNING  ----------------------------------------------------C
C----  RANDOM SEEDS FOR EACH STRUCTURE  ---------------------------------C
   
      R = 3.0       ! INITIAL RANDOM SEED FOR VORTEX TIME
      
      Do K = 1, 10000
         
         SR1S(K,1) = 2 * K - 1
         SR1S(K,2) = 2 * (K + 10000) - 1
         SR1S(K,3) = 2 * (K + 20000) - 1
         
         SR2S(K,1) = 2 * (K + 30000) - 1
         SR2S(K,2) = 2 * (K + 40000) - 1
         SR2S(K,3) = 2 * (K + 50000) - 1

         SR2SS(K,1) = 2 * (K + 60000)- 1
         SR2SS(K,2) = 2 * (K + 70000) - 1
         SR2SS(K,3) = 2 * (K + 80000) - 1
 
         SR3S(K,1) = 2 * (K + 90000) - 1
         SR3S(K,2) = 2 * (K + 100000) - 1
         SR3S(K,3) = 2 * (K + 110000) - 1

         SR4S(K,1) = 2 * (K + 120000) - 1
         SR4S(K,2) = 2 * (K + 130000) - 1
         SR4S(K,3) = 2 * (K + 140000) - 1
         
      Enddo
      
      RP1(1) = 1
      RP1(2) = 2
      RP1(3) = 3
      
      RP2(1) = 4
      RP2(2) = 5
      RP2(3) = 6
      
      RP2S(1) = 7
      RP2S(2) = 8
      RP2S(3) = 9
      
      RP3(1) = 10 
      RP3(2) = 11
      RP3(3) = 12
      
      RP4(1) = 13
      RP4(2) = 14
      RP4(3) = 15
      
      R10 = 1
      R20 = 3
      R2S = 5
      R30 = 7
      R40 = 9
      
      Do I = 1, 10000
         RT1(I) = NRAB1(R10,1,10000)
         RT2(I) = NRAB1(R20,20001,30000)
         RT2S(I) = NRAB1(R2S,30001,40000)
         RT3(I) = NRAB1(R30,40001,50000)
         RT4(I) = NRAB1(R40,50001,60000)
      Enddo
      
c===========================================================================c      
      Return
      End
      
