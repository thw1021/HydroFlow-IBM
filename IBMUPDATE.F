C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (12/05/2022)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                                                                           #
c         Update coordinates of immersed boundary (moving boundary)         #                                           
c                                                                           #
c############################################################################
      Subroutine IBMUPDATE
      
      Include './Include/OCERM_INF'

      Double precision RN,OM,BCOS,BSIN
      
      
C===========================================================================C      
C         Update position of boundary
C===========================================================================C
      
      nn = 1
      OM = 0.0125 * 3.1415926 / 180  
      
      Do J = 1,IBN(nn)
          RN = Sqrt((BXY(nn,J,1) - PXY(1020,1))**2+
     &          (BXY(nn,J,2) - PXY(1020,2))**2)
          BCOS = (BXY(nn,J,1) - PXY(1020,1)) /  RN
          BSIN = (BXY(nn,J,2) - PXY(1020,2)) /  RN
          BXY(nn,J,1) = RN * (BCOS*COS(OM)-BSIN*SIN(OM)) + PXY(1020,1)
          BXY(nn,J,2) = RN * (BSIN*COS(OM)+BCOS*SIN(OM)) + PXY(1020,2)
      Enddo

 
      End Subroutine