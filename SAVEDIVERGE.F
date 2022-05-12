C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (06/01/2022)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #     
C                 Save results after divergence of calculation              #                                           
C                                                                           #
C############################################################################
      Subroutine SAVEDIVERGE
      
      Include './Include/OCERM_INF'
      
      Character*10 FN
      Integer IJK,IJCELL
      
      
      Write (FN,'(F10.3)')  TIME * 24
      IJK = IJP * KB
      IJCELL = IJM * KBM
      Write (ITF,'(A8,A4,A10,A4,A4,I8,A4,I8,A10,A20,A20)')
     &      'ZONE T=', '"',FN,'"',' I=',IJK,' K=',IJCELL,',F=FEPOINT',
     &      ', ET=BRICK',', D=(FECONNECT)'
     
      Do K = 1, KB
      Do I = 1, IJP
          ZP = HP(I)+ELFV(I)
          If(ZP .LE. 0.0) Then
              ZP = 0.001
              ZZZ = Z(K) * ZP + (ZP - HP(I))
          Else
              ZZZ = Z(K)*(HP(I)+ELFV(I))+ELFV(I)
          Endif	       
          Write(ITF,'(20E16.5)') 
     &          PXY(I,1), PXY(I,2), ZZZ, ELFV(I), PNV(I,K),
     &          URV(I,K), VRV(I,K), WRV(I,K), VORTICITYV(I,K),
     &          VORXV(I,K), VORYV(I,K), VORZV(I,K),OMIGAV(I,K), 
     &          KMV(I,K),TKERV(I,K),TDISSRV(I,K),POREV(I,K)
      Enddo
      Enddo

      
      End Subroutine