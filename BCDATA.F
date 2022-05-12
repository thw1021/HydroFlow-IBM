C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (11/01/2010)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                     reading boundary condition                            #
c                                                                           #
c############################################################################

      Subroutine BCDATA
      Include './Include/OCERM_INF'
      Dimension COM(80)
C===========================================================================C
C             Elevation boundary condition                                  c
C===========================================================================C
      If (NUMEBC .NE. 0) Then
          Open (100, File = './Include/'//NAMEEBC, Status='unknown')
		Read (100,5800) (IEBC(N), IEBCINX(N), N = 1, NUMEBC)
          Write (IUPRT,5800) (IEBC(N), IEBCINX(N), N = 1, NUMEBC)
          Do I = 1, 100000
             Read (100,6400,End=10) TIME
		   Read (100,6400,End=10) (EBDRY(N),N = 1,NUMEBC)
	       
		   Write (IUPRT,6400) TIME
             Write (IUPRT,6400) (EBDRY(N),N = 1,NUMEBC)

             Write (IUT90,6500) TIME
             Write (IUT90,6500) (EBDRY(N),N = 1,NUMEBC)
	    Enddo
10	    Continue
		Close(100)
	Endif
C===========================================================================C
C             Velocity boundary condition                                   c
C===========================================================================C
      If (NUMVBC .NE. 0) Then
          Open (100, File = './Include/'//NAMEVBC, Status='unknown')
		Do I = 1, NUMVBC
		   Read (100,5900) IVBC(I), IVBCINX(I), 
     &		               (VVDIST(I,K), K =1, KBM)
             Write (IUPRT,5900) IVBC(I), IVBCINX(I), 
     &		                  (VVDIST(I,K), K =1, KBM)
	    Enddo
          Do I = 1, 100000
             Read (100,6400,End=20) TIME
             Read (100,6400,End=20) (UBDRY(N),N = 1,NUMVBC)
             Read (100,6400,End=20) (VBDRY(N),N = 1,NUMVBC)

		   Write (IUPRT,6400) TIME
             Write (IUPRT,6400) (UBDRY(N),N = 1,NUMVBC)
             Write (IUPRT,6400) (VBDRY(N),N = 1,NUMVBC)

             Write (IUT91,6500) TIME
             Write (IUT91,6500) (UBDRY(N),N = 1,NUMVBC)
             Write (IUT91,6500) (VBDRY(N),N = 1,NUMVBC)
	    Enddo
20	    Continue
	    Close(100)	   
	Endif
C===========================================================================C
C             River/dam and onshore intake/outfall discharge boundary       c
C===========================================================================C
      If (NUMQBC .NE. 0) Then
          Open (100, File = './Include/'//NAMEQBC, Status='unknown')
		Do I = 1, NUMQBC
		   Read (100,5900) IQBC(I), IQBCINX(I),
     &                       (VQDIST(I,K), K = 1, KBM)
     
             Write (IUPRT,5900) IQBC(I), IQBCINX(I),
     &                          (VQDIST(I,K), K = 1, KBM)
	    Enddo
          Do I = 1, 100000
             Read (100,6400,End=30) TIME
             Read (100,6400,End=30) (QDIS(N),N = 1,NUMQBC)
	
             Write (IUPRT,6400) TIME
             Write (IUPRT,6400) (QDIS(N),N = 1,NUMQBC)

             Write (IUT92,6500) TIME
             Write (IUT92,6500) (QDIS(N),N = 1,NUMQBC)
	    Enddo
30	    Continue
	    Close(100)
	Endif
C===========================================================================C
C             Offshore intake/outfall(diffuser) boundary                    c
C===========================================================================C
      If (NUMDBC .NE. 0) Then
          Open (100, File = './Include/'//NAMEDBC, Status='unknown')
		  Do I = 1, NUMDBC
		     Read (100,5910) IDBC(I), (VDDIST(I,K), K = 1, KBM)
             Write (IUPRT,5910) IDBC(I), (VDDIST(I,K), K = 1, KBM)
	    Enddo
          Do I = 1, 100000
             Read (100,6400,End=40) TIME
             Read (100,6400,End=40) (QDIFF(N), N = 1,NUMDBC)
             Read (100,6400,End=40) (QU(N), N = 1, NUMDBC)
             Read (100,6400,End=40) (QV(N), N = 1, NUMDBC)

             Write (IUPRT,6400) TIME
             Write (IUPRT,6400) (QDIFF(N),N = 1,NUMDBC)
             Write (IUPRT,6400) (QU(N),N = 1,NUMDBC)
             Write (IUPRT,6400) (QV(N),N = 1,NUMDBC)

             Write (IUT93,6500) TIME
             Write (IUT93,6500) (QDIFF(N),N = 1,NUMDBC)
             Write (IUT93,6500) (QU(N),N = 1,NUMDBC)
             Write (IUT93,6500) (QV(N),N = 1,NUMDBC)
	    Enddo
40	    Continue
	    Close(100)	   
	Endif
C===========================================================================C
C             Astrotide boundary                                            c
C===========================================================================C
      If (NUMAST.NE.0) Then
          Open (100, File = './Include/'//NAMEAST, Status='unknown')
C		Read (100,5800) (IABC(N), N = 1, NUMAST)
C          Write (IUPRT,5800) (IABC(N), N = 1, NUMAST)
          Do I = 1, NUMAST
             Read (100,6410) IABC(I), IABCINX(I), EMEAN(I)
	       Read (100,6400,Err=50) (AMP(I,J), J= 1, IAST)
	       Read (100,6400,Err=50) (PHASE(I,J), J= 1, IAST)

		   Write (IUPRT,6410) IABC(I), IABCINX(I), EMEAN(I)            
		   Write (IUPRT,6400) (AMP(I,J), J = 1, IAST)
		   Write (IUPRT,6400) (PHASE(I,J), J = 1, IAST)

		   Write (IUT94,6410) IABC(I), IABCINX(I), EMEAN(I)            
		   Write (IUT94,6400) (AMP(I,J), J = 1, IAST)
		   Write (IUT94,6400) (PHASE(I,J), J = 1, IAST)
	    Enddo
50	    Continue
	    Close(100)
          Call GETCOR(IYR,IMO,IDA,AMP0,PHASE0,FREQ,11)
   	Endif
C===========================================================================C
C             Viscosity boundary                                            c
C===========================================================================C
C----- One equation (SA)
	If (VERTMIX .EQ. 'SAMODEL   ' .OR. 
     *    DES .EQ. 'SADES  ' .OR. DES .EQ. 'SADDES ' .OR.
     *    DES .EQ. 'SAIDDES' .OR. DES .EQ. 'SAZDES ') Then
         If (NUMEBC .NE. 0) Then
 		  Open (100, File = './Include/VIS.EBC', Status='unknown')
            Do I = 1, 100000
               Read (100,6400,End=60) TIME
		     Read (100,6400,End=60) (VISE(N), N = 1, NUMEBC)
	       
		     Write (IUPRT,6400) TIME
               Write (IUPRT,6400) (VISE(N), N = 1, NUMEBC)

               Write (IUT96,6500) TIME
               Write (IUT96,6500) (VISE(N), N = 1, NUMEBC)
	      Enddo
60	      Continue
		  Close(100)
	    Endif
	    If (NUMQBC .NE. 0) Then
             Open (100, File = './Include/VIS.QBC', Status='unknown')
             Do I = 1, 100000
                Read (100,6400,End=70) TIME
	          Do N = 1, NUMQBC
                   Read (100,6450,End=70) (VISQ(N,K), K = 1, KBM)
	          Enddo
	
                Write (IUPRT,6400) TIME
	          Do N = 1, NUMQBC
                   Write (IUPRT,6500) (VISQ(N,K), K = 1, KBM)
	          Enddo

                Write (IUT97,6500) TIME
	          Do N = 1, NUMQBC
                   Write (IUT97,6500) (VISQ(N,K), K = 1, KBM)
	          Enddo
	       Enddo
70	       Continue
	       Close(100)
	    Endif
	Endif
c---- Two equation (SST k-w)	
	If (VERTMIX .EQ. 'SSTMODEL  ' .OR. DES .EQ. 'SSTDES ') Then
         If (NUMEBC .NE. 0) Then
 		  Open (100, File = './Include/VIS.EBC', Status='unknown')
            Do I = 1, 100000
               Read (100,6400,End=80) TIME
		       Read (100,6400,End=80) (TKEE(N), N = 1, NUMEBC)
		       Read (100,6400,End=80) (TDISSE(N), N = 1, NUMEBC)
	       
		       Write (IUPRT,6400) TIME
               Write (IUPRT,6400) (TKEE(N), N = 1, NUMEBC)
               Write (IUPRT,6400) (TDISSE(N), N = 1, NUMEBC)

               Write (IUT96,6500) TIME
               Write (IUT96,6500) (TKEE(N), N = 1, NUMEBC)
               Write (IUT96,6500) (TDISSE(N), N = 1, NUMEBC)
               
	      Enddo
80	      Continue
		  Close(100)
	    Endif
	    If (NUMQBC .NE. 0) Then
             Open (100, File = './Include/VIS.QBC', Status='unknown')
             Do I = 1, 100000
                Read (100,6400,End=90) TIME
	          Do N = 1, NUMQBC
                   Read (100,6450,End=90) (TKEQ(N,K), K = 1, KBM)
                   Read (100,6450,End=90) (TDISSQ(N,K), K = 1, KBM)
	          Enddo
	
                Write (IUPRT,6400) TIME
	          Do N = 1, NUMQBC
                   Write (IUPRT,6500) (TKEQ(N,K), K = 1, KBM)
                   Write (IUPRT,6500) (TDISSQ(N,K), K = 1, KBM)
	          Enddo

                Write (IUT97,6500) TIME
	          Do N = 1, NUMQBC
                   Write (IUT97,6500) (TKEQ(N,K), K = 1, KBM)
                   Write (IUT97,6500) (TDISSQ(N,K), K = 1, KBM)
	          Enddo
	       Enddo
90	       Continue
	       Close(100)
	    Endif
	Endif
C===========================================================================C
C             Concentration boundary                                        c
C===========================================================================C
	If (WAQ .NE. 'NEGLECT') Then
         If (NUMEBC .NE. 0) Then
 		  Open (100, File = './Include/CBC.EBC', Status='unknown')
            Do I = 1, 100000
               Read (100,6400,End=100) TIME
		       Read (100,6400,End=100) (CEBC(N), N = 1, NUMEBC)
	       
		       Write (IUPRT,6400) TIME
               Write (IUPRT,6400) (CEBC(N), N = 1, NUMEBC)

               Write (IUT98,6500) TIME
               Write (IUT98,6500) (CEBC(N), N = 1, NUMEBC)
	      Enddo
100	      Continue
		  Close(100)
	    Endif
	    If (NUMQBC .NE. 0) Then
             Open (100, File = './Include/CBC.QBC', Status='unknown')
             Do I = 1, 100000
                Read (100,6400,End=101) TIME
	          Do N = 1, NUMQBC
                   Read (100,6450,End=101) (CQBC(N,K), K = 1, KBM)
	          Enddo
	
                Write (IUPRT,6400) TIME
	          Do N = 1, NUMQBC
                   Write (IUPRT,6500) (CQBC(N,K), K = 1, KBM)
	          Enddo

                Write (IUT99,6500) TIME
	          Do N = 1, NUMQBC
                   Write (IUT99,6500) (CQBC(N,K), K = 1, KBM)
	          Enddo
	       Enddo
101	       Continue
	       Close(100)
	    Endif
	    If (NUMDBC .NE. 0) Then
             Open (100, File = './Include/CBC.DBC', Status='unknown')
             Do I = 1, 100000
                Read (100,6400,End=102) TIME
                Read (100,6450,End=102) (CDBC(N), N = 1, NUMDBC)
	
                Write (IUPRT,6400) TIME
                Write (IUPRT,6500) (CDBC(N), N = 1, NUMDBC)
                  
                Write (IUT101,6500) TIME
                Write (IUT101,6500) (CDBC(N), N = 1, NUMDBC)
                
	       Enddo
102	       Continue
	       Close(100)
	    Endif
	Endif
C===========================================================================C
c======  reading windfield  ================================================C
CC	If(WIND.NE.'NEGLECT') Then
CC	   Open(1111,File='windfield.dat',status='unknown')
CC	   Do 200 K = 1, 100000
CC            Read (1111,6600,End=210) TIME
CC	      Do I = 1, IM
CC	      Do J = 1, JM
CC               Read (1111,6601) I1, J1, WINDU(I,J), WINDV(I,J)
CC	      Enddo
CC	      Enddo
CC            Write (IUT95,6600) TIME
CC	      Do J = 1, JM
CC	      Do I = 1, IM
CC               Write (IUT95,6602) WINDU(I,J), WINDV(I,J)
CC	      Enddo
CC	      Enddo
CC200      Continue
CC210      Continue
CC	   Close(1111)
CC	Endif
c===============   end reading   ============================================
        Return
5000    Format (80A1)
5100    Format (/1X,80A1/)
5200    Format (' ksl = ',I5,/)
5300    Format (//' number of standard levels in run_data',I5,' (iksl)'/
     *    '           do not equal'/
     *    ' number of standard levels in GEFM_INC ',I5,' (ksl)'/
     *    ' please correct this problem and try again'//)
5400   Format (I5,1X,A20)
5500   Format (2I5,4F10.5)
5600   Format (4I5,1F10.5)
5700   Format (11F10.5)
5800   Format (8I14)
5900   Format (2I14,200F14.6)
5910   Format ( I14,200F14.6)
6000   Format (4I5,20F5.1)
6100   Format (4I5,/,20F5.1)
6200   Format (2I5,20F5.1)
6300   Format (2I5,/,20F5.1)
6400   Format (8F14.6)
6450   Format (8E14.6)
6410   Format ( 2I14,8F14.6)
6402   Format (F14.6)
6401   Format (8F9.2)
6500   Format (8E14.7)
6600   Format (F14.6)
6601	 Format (2I5,2F14.3)
6602   Format (2F14.3) 
       End
