c-----------------------------------------------------------------------------------c
c                   SUBROUTINE PROGRAM                                              #
C                 VERSION 1.0 (13/04/2015)                                          #
C                 AUTHORIZED BY CHE HAIOU                                           #
C                               SHANGHAI JIAO TONG UNIVERSITY                       #
C                               SHANGHAI, CHINA                                     #
c-----------------------------------------------------------------------------------#
c                       DEM Normal vector  cells' surface                           #
c-----------------------------------------------------------------------------------c
      Subroutine Normal
      Include './Include/OCERM_INF'
      Parameter(NSCHEME=1)
      Double Precision X1DEM,Y1DEM,ZS1DEM,ZB1DEM
      Double Precision X2DEM,Y2DEM,ZS2DEM,ZB2DEM
      Double Precision X3DEM,Y3DEM,ZS3DEM,ZB3DEM
      Double Precision X4DEM,Y4DEM,ZS4DEM,ZB4DEM
      Double Precision XX1DEM,XY1DEM,XZ1DEM,XN1DEM
      Double Precision XX2DEM,XY2DEM,XZ2DEM,XN2DEM               

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(J,X1DEM,Y1DEM,ZS1DEM,
!$OMP&         ZB1DEM,X2DEM,Y2DEM,ZS2DEM,ZB2DEM,X3DEM,Y3DEM,
!$OMP&         ZS3DEM,ZB3DEM,X4DEM,Y4DEM,ZS4DEM,ZB4DEM,XX1DEM,
!$OMP&         XY1DEM,XZ1DEM,XN1DEM,XX2DEM,XY2DEM,XZ2DEM,XN2DEM)      
!$OMP DO
      Do J = 1 , IJM  
       IF(NSCHEME .EQ. 1 .OR. Cell_POLYGEN(J) .EQ. 3)THEN   
        X1DEM = PXY(CELL_INDEX(J,1),1)
        X2DEM = PXY(CELL_INDEX(J,2),1)                                  
        X3DEM = PXY(CELL_INDEX(J,3),1)
        X4DEM = PXY(CELL_INDEX(J,Cell_POLYGEN(J)),1)
        
        Y1DEM = PXY(CELL_INDEX(J,1),2)
        Y2DEM = PXY(CELL_INDEX(J,2),2)
        Y3DEM = PXY(CELL_INDEX(J,3),2)
        Y4DEM = PXY(CELL_INDEX(J,Cell_POLYGEN(J)),2)
        
        ZS1DEM= ELFV(CELL_INDEX(J,1))  
        ZS2DEM= ELFV(CELL_INDEX(J,2))
        ZS3DEM= ELFV(CELL_INDEX(J,3))
        ZS4DEM= ELFV(CELL_INDEX(J,Cell_POLYGEN(J)))
        
        ZB1DEM= - HP(CELL_INDEX(J,1))
        ZB2DEM= - HP(CELL_INDEX(J,2))
        ZB3DEM= - HP(CELL_INDEX(J,3))        
        ZB4DEM= - HP(CELL_INDEX(J,Cell_POLYGEN(J))) 
c-----
        XX1DEM = (Y3DEM - Y1DEM) * (ZS4DEM - ZS2DEM) - 
     &           (ZS3DEM - ZS1DEM) * (Y4DEM - Y2DEM)
        XY1DEM = (ZS3DEM - ZS1DEM) * (X4DEM - X2DEM) - 
     &           (X3DEM - X1DEM) * (ZS4DEM - ZS2DEM)
        XZ1DEM = (X3DEM - X1DEM) * (Y4DEM - Y2DEM) - 
     &           (Y3DEM - Y1DEM) * (X4DEM - X2DEM)
        XX2DEM = (Y4DEM - Y2DEM) * (ZB3DEM - ZB1DEM) - 
     &           (ZB4DEM - ZB2DEM) * (Y3DEM - Y1DEM)
        XY2DEM = (ZB4DEM - ZB2DEM) * (X3DEM - X1DEM) - 
     &           (X4DEM - X2DEM) * (ZB3DEM - ZB1DEM)
        XZ2DEM = (X4DEM - X2DEM) * (Y3DEM - Y1DEM) - 
     &           (Y4DEM - Y2DEM) * (X3DEM - X1DEM)
c----- 
        XN1DEM = Sqrt(XX1DEM**2 + XY1DEM**2 + XZ1DEM**2)
        XN2DEM = Sqrt(XX2DEM**2 + XY2DEM**2 + XZ2DEM**2)
c-----       
        DEMCOSA(J,1) = XX1DEM / XN1DEM
        DEMCOSA(J,2) = XY1DEM / XN1DEM
        DEMCOSA(J,3) = XZ1DEM / XN1DEM
        DEMCOSB(J,1) = XX2DEM / XN2DEM
        DEMCOSB(J,2) = XY2DEM / XN2DEM
        DEMCOSB(J,3) = XZ2DEM / XN2DEM 
       
       ENDIF 
        
       IF(NSCHEME .EQ. 2 .AND. Cell_POLYGEN(J) .EQ. 4)THEN 

        X1DEM = PXY(CELL_INDEX(J,1),1)
        X2DEM = PXY(CELL_INDEX(J,2),1)                                  
        X3DEM = PXY(CELL_INDEX(J,3),1)
        X4DEM = PXY(CELL_INDEX(J,4),1)
        
        Y1DEM = PXY(CELL_INDEX(J,1),2)
        Y2DEM = PXY(CELL_INDEX(J,2),2)
        Y3DEM = PXY(CELL_INDEX(J,3),2)
        Y4DEM = PXY(CELL_INDEX(J,4),2)
        
        ZS1DEM= ELFV(CELL_INDEX(J,1))  
        ZS2DEM= ELFV(CELL_INDEX(J,2))
        ZS3DEM= ELFV(CELL_INDEX(J,3))
        ZS4DEM= ELFV(CELL_INDEX(J,4))
        
        ZB1DEM= - HP(CELL_INDEX(J,1))
        ZB2DEM= - HP(CELL_INDEX(J,2))
        ZB3DEM= - HP(CELL_INDEX(J,3))        
        ZB4DEM= - HP(CELL_INDEX(J,4))       
        
        XX1DEM = (Y2DEM - Y1DEM) * (ZS3DEM - ZS2DEM) - 
     &           (Y3DEM - Y2DEM) * (ZS2DEM - ZS1DEM)       
        XY1DEM = (X3DEM - X2DEM) * (ZS2DEM - ZS1DEM) - 
     &           (X2DEM - X1DEM) * (ZS3DEM - ZS2DEM)        
        XZ1DEM = (X2DEM - X1DEM) * (Y3DEM - Y2DEM) - 
     &           (X3DEM - X2DEM) * (Y2DEM - Y1DEM)        
        XX2DEM = (Y2DEM - Y1DEM) * (ZB3DEM - ZB2DEM) - 
     &           (Y3DEM - Y2DEM) * (ZB2DEM - ZB1DEM)        
        XY2DEM = (X3DEM - X2DEM) * (ZB2DEM - ZB1DEM) - 
     &           (X2DEM - X1DEM) * (ZB3DEM - ZB2DEM)        
        XZ2DEM = (X2DEM - X1DEM) * (Y3DEM - Y2DEM) - 
     &           (X3DEM - X2DEM) * (Y2DEM - Y1DEM)
        
c----- 
        XN1DEM = Sqrt(XX1DEM**2 + XY1DEM**2 + XZ1DEM**2)
        XN2DEM = Sqrt(XX2DEM**2 + XY2DEM**2 + XZ2DEM**2)
c-----       
        DEMCOSA(J,1) = XX1DEM / XN1DEM
        DEMCOSA(J,2) = XY1DEM / XN1DEM
        DEMCOSA(J,3) = XZ1DEM / XN1DEM
        DEMCOSB(J,1) = XX2DEM / XN2DEM
        DEMCOSB(J,2) = XY2DEM / XN2DEM
        DEMCOSB(J,3) = XZ2DEM / XN2DEM   
        
c----- 
        XX1DEM = (Y3DEM - Y1DEM) * (ZS4DEM - ZS3DEM) - 
     &           (Y4DEM - Y3DEM) * (ZS3DEM - ZS1DEM)       
        XY1DEM = (X4DEM - X3DEM) * (ZS3DEM - ZS1DEM) - 
     &           (X3DEM - X1DEM) * (ZS4DEM - ZS3DEM)        
        XZ1DEM = (X3DEM - X1DEM) * (Y4DEM - Y3DEM) - 
     &           (X4DEM - X3DEM) * (Y3DEM - Y1DEM)        
        XX2DEM = (Y3DEM - Y1DEM) * (ZB4DEM - ZB3DEM) - 
     &           (Y4DEM - Y3DEM) * (ZB3DEM - ZB1DEM)     
        XY2DEM = (X4DEM - X3DEM) * (ZB3DEM - ZB1DEM) - 
     &           (X3DEM - X1DEM) * (ZB4DEM - ZB3DEM)        
        XZ2DEM = (X3DEM - X1DEM) * (Y4DEM - Y3DEM) - 
     &           (X4DEM - X3DEM) * (Y3DEM - Y1DEM)      
        
c----- 
        XN1DEM = Sqrt(XX1DEM**2 + XY1DEM**2 + XZ1DEM**2)
        XN2DEM = Sqrt(XX2DEM**2 + XY2DEM**2 + XZ2DEM**2)
c-----       
        DEMCOSA(J,4) = XX1DEM / XN1DEM
        DEMCOSA(J,5) = XY1DEM / XN1DEM
        DEMCOSA(J,6) = XZ1DEM / XN1DEM
        DEMCOSB(J,4) = XX2DEM / XN2DEM
        DEMCOSB(J,5) = XY2DEM / XN2DEM
        DEMCOSB(J,6) = XZ2DEM / XN2DEM         
       
       ENDIF
       
      Enddo
!$OMP END DO
!$OMP END PARALLEL               
      Return
      End
