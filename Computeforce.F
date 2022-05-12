      Subroutine Computeforce 
      Include './Include/OCERM_INF'
C     Parameter (DCYC = 0.133,WIDTH = 0.006)
C     Parameter (U0 = 0.088)


      Common/CF3/ TWVX(NEdge,KB),TWVY(NEdge,KB),TWVZ(NEdge,KB),
     &           TWV(NEdge,KB),PWVX(NEdge,KB),PWVY(NEdge,KB),
     &           PWV(NEdge,KB), DELA(NEdge,KB),W_NORMAL(NCYC,NEdge,2),
     &           TWVX_AVE(NEdge,KB),TWVY_AVE(NEdge,KB),
     &           TWVZ_AVE(NEdge,KB),TWV_AVE(NEdge,KB),TWV_DIR(NEdge,KB),
     &           PWVX_AVE(NEdge,KB),PWVY_AVE(NEdge,KB),PWV_AVE(NEdge,KB)
      Common/CF4/  FDRAG(3,NEdge,KB), FLIFT(3,NEdge,KB), 
     &            FDRAG_AVE(3,NCYC), FLIFT_AVE(3,NCYC),
     &            DIS_WALL(NCYC,NEdge),CELL_LENGTH(NCYC,NEdge)
      Common/CF5/  TWX_AVE(NCYC,NEdge,KB), TWY_AVE(NCYC,NEdge,KB), 
     &            TWZ_AVE(NCYC,NEdge,KB), PWX_AVE(NCYC,NEdge,KB),
     &            PWY_AVE(NCYC,NEdge,KB), PW_AVE(NCYC,NEdge,KB)

      Character * 4 FN

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K)      
C----    INFORMATIONS OF THE WALL CELL  -------------------------------------C
!$OMP DO
         Do I = 1, NUM_WALL
            Do J = 1, IWALL(I)
               W_NORMAL(I,J,1) = CELL_CUV(NUM_WALLCELL(I,J),1,7) !cos
               W_NORMAL(I,J,2) = CELL_CUV(NUM_WALLCELL(I,J),1,8) !sin
               CELL_LENGTH(I,J) = CELL_CUV(NUM_WALLCELL(I,J),1,6)
             DIS_WALL(I,J) = .5*AREA(NUM_WALLCELL(I,J))/CELL_LENGTH(I,J)
            Enddo
         Enddo
!$OMP END DO          
C----    WALL FORCE CALCULATION ---------------------------------------------C
       
         DO K=1,KB
            DO I=1,3
!$OMP DO                 
              DO J=1,NEdge
                 FDRAG(I,J,K) = 0.0
                 FLIFT(I,J,K) = 0.0  
              ENDDO
!$OMP END DO               
             ENDDO
         ENDDO

         UVC = 0.0
         
         DO I=1,3
!$OMP DO             
           DO J=1,NCYC
             FDRAG_AVE(I,J) = 0.0
             FLIFT_AVE(I,J) = 0.0
           ENDDO
!$OMP END DO              
         ENDDO
!$OMP END PARALLEL


!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP&PRIVATE(I,J,K,INDEXCELL,I1,I2,I3,I4,NP1,NP2,DD,UUU,VVV,WWW,
!$OMP&        TWV,TWVX,TWVY,TWVZ,PWV,PWVX,PWVY,FDRAG,FLIFT)
         TTT = 0.0
!$OMP DO         
         Do J = 1, NUM_WALL
            Do I = 1, IWALL(J)
               INDEXCELL = NUM_WALLCELL(J,I)   !index of points of the cell
               I1 = CELL_SIDE(INDEXCELL,1,3)
               I2 = CELL_SIDE(INDEXCELL,2,3)
               I3 = CELL_SIDE(INDEXCELL,3,3)
               I4 = CELL_SIDE(INDEXCELL,4,3)
               
               NP1 = CELL_SIDE(INDEXCELL,1,3)
               NP2 = CELL_SIDE(INDEXCELL,1,4)!该单元第1条边的两个顶点编号
               DD = (HP(NP1)+HP(NP2))/2.+ !该边中点总水深
     &              (ELFV(NP1)+ELFV(NP2))/2.

               Do K = 1, KBM
               UUU = ((URV(I1,K)+URV(I2,K)+URV(I3,K)+URV(I4,K))/4.+
     &         (URV(I1,K+1)+URV(I2,K+1)+URV(I3,K+1)+URV(I4,K+1))/4.)/2.
               VVV = ((VRV(I1,K)+VRV(I2,K)+VRV(I3,K)+VRV(I4,K))/ 4.+
     &         (VRV(I1,K+1)+VRV(I2,K+1)+VRV(I3,K+1)+VRV(I4,K+1))/4.)/2. 
               WWW = ((WRV(I1,K)+WRV(I2,K)+WRV(I3,K)+WRV(I4,K))/ 4.+
     &         (WRV(I1,K+1)+WRV(I2,K+1)+WRV(I3,K+1)+WRV(I4,K+1))/4.)/2. !单元中心速度
      
               TWV(I,K) = 1.E-3*(VVV * W_NORMAL(J,I,1) -  !牛顿内摩擦定律求切应力
     &                           UUU * W_NORMAL(J,I,2))/DIS_WALL(J,I)
              
               TWVX(I,K) = -TWV(I,K) *W_NORMAL(J,I,2)
               TWVY(I,K) =  TWV(I,K) *W_NORMAL(J,I,1) 
               TWVZ(I,K) = 1.E-3 * WWW / DIS_WALL(J,I)
               
               
               PWV(I,K) = ((PNV(NP1,K)+PNV(NP2,K))/2. + 9800.0*
     &            (-DD * ZZ(K)))   !plus of hydrostatic and non-hydrostatic pressure
C     &            *CELL_LENGTH(J,I)*(HP(INDEXCELL)+ELFV(INDEXCELL))*
C     &            DZ(K)
              
               PWVX(I,K) = PWV(I,K) * W_NORMAL(J,I,1)
               PWVY(I,K) = PWV(I,K) * W_NORMAL(J,I,2)
               
               Enddo  
           Enddo
           Do I = 1, IWALL(J)
              INDEXCELL = NUM_WALLCELL(J,I)
              NP1 = CELL_SIDE(INDEXCELL,1,3)
              NP2 = CELL_SIDE(INDEXCELL,1,4)
              DD = (HP(NP1)+HP(NP2))/2.+
     &             (ELFV(NP1)+ELFV(NP2))/2.
              Do K = 1, KBM
                 FDRAG(1,I,K) = (TWVX(I,K) + PWVX(I,K)) *          ! TOTAL DRAG
     &                          CELL_LENGTH(J,I) * DD * DZ(K)
                 FLIFT(1,I,K) = (TWVY(I,K) + PWVY(I,K)) *          ! TOTAL LIFT
     &                          CELL_LENGTH(J,I) * DD * DZ(K)
     
                 FDRAG(2,I,K) = (TWVX(I,K)) *                      ! FRICTION
     &                          CELL_LENGTH(J,I) * DD * DZ(K)
                 FLIFT(2,I,K) = (TWVY(I,K)) *                      ! FRICTION
     &                          CELL_LENGTH(J,I) * DD * DZ(K)
     
                 FDRAG(3,I,K) = (PWVX(I,K)) *                      ! PN
     &                          CELL_LENGTH(J,I) * DD * DZ(K)
                 FLIFT(3,I,K) = (PWVY(I,K)) *                      ! PN
     &                          CELL_LENGTH(J,I) * DD * DZ(K)

              Enddo  
           Enddo
           Do I = 1, IWALL(J)
              FDRAG(1,I,KB) = 0.0
              FLIFT(1,I,KB) = 0.0
              FDRAG(2,I,KB) = 0.0
              FLIFT(2,I,KB) = 0.0
              FDRAG(3,I,KB) = 0.0
              FLIFT(3,I,KB) = 0.0
              
              Do K = 1, KBM
                 FDRAG(1,I,KB) = FDRAG(1,I,KB) + FDRAG(1,I,K)
                 FLIFT(1,I,KB) = FLIFT(1,I,KB) + FLIFT(1,I,K)
                 FDRAG(2,I,KB) = FDRAG(2,I,KB) + FDRAG(2,I,K)
                 FLIFT(2,I,KB) = FLIFT(2,I,KB) + FLIFT(2,I,K)
                 FDRAG(3,I,KB) = FDRAG(3,I,KB) + FDRAG(3,I,K)
                 FLIFT(3,I,KB) = FLIFT(3,I,KB) + FLIFT(3,I,K)
              Enddo           
           Enddo
           Do I = 1, IWALL(J)
              FDRAG_AVE(1,J) = FDRAG_AVE(1,J) + FDRAG(1,I,KB)
              FLIFT_AVE(1,J) = FLIFT_AVE(1,J) + FLIFT(1,I,KB)
              FDRAG_AVE(2,J) = FDRAG_AVE(2,J) + FDRAG(2,I,KB)
              FLIFT_AVE(2,J) = FLIFT_AVE(2,J) + FLIFT(2,I,KB)
              FDRAG_AVE(3,J) = FDRAG_AVE(3,J) + FDRAG(3,I,KB)
              FLIFT_AVE(3,J) = FLIFT_AVE(3,J) + FLIFT(3,I,KB)
           Enddo
c          FDRAG_AVE(1,J) = FDRAG_AVE(1,J)  !/ (U0**2.*500.*WIDTH*DCYC)
c          FLIFT_AVE(1,J) = FLIFT_AVE(1,J)  !/ (U0**2.*500.*WIDTH*DCYC)
c          FDRAG_AVE(2,J) = FDRAG_AVE(2,J)  !/ (U0**2.*500.*WIDTH*DCYC)
c          FLIFT_AVE(2,J) = FLIFT_AVE(2,J)  !/ (U0**2.*500.*WIDTH*DCYC)
c          FDRAG_AVE(3,J) = FDRAG_AVE(3,J)  !/ (U0**2.*500.*WIDTH*DCYC)
c          FLIFT_AVE(3,J) = FLIFT_AVE(3,J)  !/ (U0**2.*500.*WIDTH*DCYC)

C           FDRAG_AVE(J) = FDRAG_AVE(J) / (U0**2.*500.*WIDTH*DCYC)
C           FLIFT_AVE(J) = FLIFT_AVE(J) / (U0**2.*500.*WIDTH*DCYC)
      Enddo   
!$OMP END DO      
!$OMP END PARALLEL      
      Return
      End
