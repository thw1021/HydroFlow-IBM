C############################################################################
C                                                                           #
C                   SUBROUTINE PROGRAM                                      #                                   
C                 VERSION 1.0 (12/05/2022)                                  #
C                 AUTHORIZED BY LUO PENGXUAN                                #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
C---------------------------------------------------------------------------#
C                                                                           #     
C      Read data of IBM and modify parameters of mesh geometry in 3D        #                                           
C                                                                           #
C############################################################################
      Subroutine IBMINIT
      
      Include './Include/OCERM_INF'
      
      Dimension COM(80)
      Integer ISNBR(IJM)
      Integer I2,I3,NUM
      Character*10 BOUNDARY_ID
      
      
C===========================================================================C      
C         Initialize the arrays 
C===========================================================================C
      
      Do K = 1,KBM
          Do J = 1,IJM
              IDC(I,K) = 0
              ZRC(I,K) = 0
          Enddo
      Enddo
      
      Do nn = 1,NUMIB     
          Do I=1,IJM
              BXY(nn,I,1) = 0
              BXY(nn,I,2) = 0
              BXY(nn,I,3) = 0
          Enddo
      Enddo
      
      Do I = 1,IJM
          ISNBR(I) = 0
          Do N = 1,30
              NBR(I,N) = 0
          Enddo
      Enddo
      
      Do K = 1,KB
          Do I = 1,IJP
              IDP(I,K) = 0
              ZRP(I,K) = 0
		    SAU(I,K) = 0
	 	    SAV(I,K) = 0
		    SAW(I,K) = 0
          Enddo
      Enddo    
      
C===========================================================================C      
C         Read data of immersed boundary method 
C===========================================================================C

      Open (IBMSET,file='./Include/IBM/IBM_SETTING.dat')
      
      Do N = 1,3
          Read (IBMSET,'(100A1)') (COM(I),I = 1,80)
      Enddo
      
      Read (IBMSET,1000) (COM(I),I = 1,80)
      Write(IUPRT,1000) (COM(I),I = 1,80)
      Read (IBMSET,'(3I16)') IBMSCHEME,IBMDOMAIN,NUMIB
      Write(IUPRT,'(3I16)') IBMSCHEME,IBMDOMAIN,NUMIB
      Read (IBMSET,1000) (COM(I),I = 1,80)
      Write(IUPRT,1000) (COM(I),I = 1,80)
      Read (IBMSET,'(3I16)') IBFSCHEME,IBGCSCHEME,IBGCMODE
      Write(IUPRT,'(3I16)') IBFSCHEME,IBGCSCHEME,IBGCMODE
      Read (IBMSET,1000) (COM(I),I = 1,80)
      Write(IUPRT,1000) (COM(I),I = 1,80)
      Read (IBMSET,'(2F16.8,I16)') RITP,RIFL,INFCELL
      Write (IUPRT,'(2F16.8,I16)') RITP,RIFL,INFCELL
      
      Close (IBMSET)
          
      If (IBMDOMAIN.EQ.2) Then 
          Open(IBM2D,File = './Include/IBM/IBM2D.dat')
          Do nn = 1,NUMIB     
              Read(IBM2D,*) BOUNDARY_ID       ! INDEX of boundary
              Read(IBM2D,*) IBN(nn)      ! Number of IBN points
              Do J = 1,IBN(nn)
                  Read(IBM2D,*) BXY(nn,J,1),BXY(nn,J,2)
              Enddo
          Enddo
          Close(IBM2D)
      Endif
          
      If (IBMDOMAIN.EQ.3) Then        
          Open(IBM3D,File = './Include/IBM/IBM3D.dat')
          Do nn = 1,NUMIB     ! Iterate each boundary
              Read(IBM3D,*) BOUNDARY_ID       ! Index of boundary
              Read(IBM3D,*) IBN(nn)      ! Number of IBN points
              Do J = 1,IBN(nn)
                  Read(IBM3D,*) BXY(nn,J,1),BXY(nn,J,2),BXY(nn,J,3)
              Enddo            
          Enddo
          Close(IBM3D)
      Endif

C===========================================================================C      
C         Set neighbor cells
C===========================================================================C
      
      Do I = 1,IJM
          ISNBR(I) = 1
	    NUM = 0
	    Do J = 1,CELL_POLYGEN(I)
              If(CELL_SIDE(I,J,2).GT.0) Then
                  ISNBR(CELL_SIDE(I,J,2)) = 1
		        NUM = NUM + 1
                  NBR(I,NUM) = CELL_SIDE(I,J,2)
              Endif
          Enddo
          
	    Do I2 = 1,IJM
		    If(ISNBR(I2).EQ.1) Then
		    Do J = 1,CELL_POLYGEN(I2)
			    If((CELL_SIDE(I2,J,2).GT.0).AND.
     &               (ISNBR(CELL_SIDE(I2,J,2)).EQ.0)) Then
				    ISNBR(CELL_SIDE(I2,J,2)) = 2
				    NUM = NUM + 1
                      NBR(I,NUM) = CELL_SIDE(I2,J,2)
			    Endif
              Enddo
		    Endif
	    Enddo
	  
	    Do I3 = 1,IJM
		    If(ISNBR(I3).EQ.2) Then
		    Do J = 1,CELL_POLYGEN(I3)
			    If((CELL_SIDE(I3,J,2).GT.0).AND.
     &               (ISNBR(CELL_SIDE(I3,J,2)).EQ.0)) Then
				    ISNBR(CELL_SIDE(I3,J,2)) = 3
				    NUM = NUM + 1
                      NBR(I,NUM) = CELL_SIDE(I3,J,2)
			    Endif
              Enddo
		    Endif
	    Enddo
	   	 
	    Do N = 1,IJM
		    ISNBR(N) = 0
          Enddo  
      Enddo

      
1000  Format (80A1)
      
      
      End Subroutine