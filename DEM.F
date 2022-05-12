C####################################################################################
c                                                                                   #
c                   SUBROUTINE DEM PROGRAM                                          # 
C                 VERSION 1.0 (29/12/2014)                                          #
C                 AUTHORIZED BY ZHANG JINGXIN                                       #
C                               SHANGHAI JIAO TONG UNIVERSITY                       #
C                               SHANGHAI, CHINA                                     #
c-----------------------------------------------------------------------------------#
c      Simulation of particles motion by Discreted Element Method                   #
c      The DEM model is cooperated with hydrodynamic model to carry out             #
c      the numerical simulation of sediment transportation                          #
c      Model basis ---- Lagrangian method                                           #
c      Submodel ---- soft particle model                                            #
c      Submodel ---- solid particle model                                           #
c      Variables index                                                              #
c                                                                                   #
c      Subroutine index                                                             #
c                                                                                   #
c####################################################################################
      Subroutine DEMMODEL
      Include './Include/OCERM_INF'
c-----------------------------------------------------------------------------------c
c      If(DEM .EQ. 'DEMDARD') Call HARDMODULE
c      If(DEM .EQ. 'DEMSOFT') Call SOFTMODULE 

c-----------------------------------------------------------------------------------c


c-----------------------------------------------------------------------------------c

      
c-----------------------------------------------------------------------------------c
      Return
      End
