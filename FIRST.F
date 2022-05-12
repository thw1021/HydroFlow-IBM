C############################################################################
c                                                                           #
c                   SUBROUTINE PROGRAM                                      # 
C                 VERSION 1.0 (11/07/2009)                                  #
C                 AUTHORIZED BY ZHANG JINGXIN                               #
C                               SHANGHAI JIAO TONG UNIVERSITY               #
C                               SHANGHAI, CHINA                             #
c---------------------------------------------------------------------------#
c                     initial values specifying                             #
c                                                                           #
c############################################################################
      Subroutine FIRST
      Include './Include/OCERM_INF'
C===========================================================================C
C             read elevation boundary condition                             c
c===========================================================================c
      If (NUMEBC .NE. 0) Then
          Rewind IUT90
          Do I = 1, 100000
		   Read (IUT90,5100,Err=170) T2E
             Read (IUT90,5100) (DEBDRY(N,2), N = 1, NUMEBC)
             If (THOUR.LT.T2E) Go To 10
             T1E = T2E
             Do N = 1, NUMEBC
                DEBDRY(N,1) = DEBDRY(N,2)
             Enddo
	    Enddo
10        Continue
      End If
c===========================================================================c
c              read velocity boundary condition                             c
c===========================================================================c
      If (NUMVBC .NE. 0) Then
          Rewind IUT91
          Do I = 1, 100000
		   Read (IUT91,5100,Err=170) T2V
             Read (IUT91,5100) (DUBDRY(N,2), N = 1, NUMVBC)
             Read (IUT91,5100) (DVBDRY(N,2), N = 1, NUMVBC)
             If (THOUR.LT.T2V) Go To 20
             T1V = T2V
             Do N = 1, NUMVBC
                DUBDRY(N,1) = DUBDRY(N,2)
                DVBDRY(N,1) = DVBDRY(N,2)
             Enddo
	    Enddo
20        Continue
      End If
c===========================================================================c
c              read discharge boundary condition                            c
c===========================================================================c
      If (NUMQBC .NE. 0) Then
          Rewind IUT92
          Do I = 1, 100000
		     Read (IUT92,5100,Err=170) T2Q
             Read (IUT92,5100) (DQDIS(N,2), N = 1, NUMQBC)
             If (THOUR.LT.T2Q) Go To 30
             T1Q = T2Q
             Do N = 1, NUMQBC
                DQDIS(N,1) = DQDIS(N,2)
	       Enddo
	    Enddo
30	    Continue
	Endif
c===========================================================================c
c              read offshore discharge boundary condition                   c
c===========================================================================c
      If (NUMDBC .NE. 0) Then
        Rewind IUT93
        Do I = 1, 100000
           Read (IUT93,5100,Err=170) T2D
           Read (IUT93,5100) (DQDIFF(N,2), N = 1, NUMDBC)
           Read (IUT93,5100) (DQU(N,2), N = 1, NUMDBC)
           Read (IUT93,5100) (DQV(N,2), N = 1, NUMDBC)
           If (THOUR.LT.T2D) Go To 40
           T1D = T2D
           Do N = 1, NUMDBC
              DQDIFF(N,1) = DQDIFF(N,2)
              DQU(N,1) = DQU(N,2)
              DQV(N,1) = DQV(N,2)
           Enddo
	  Enddo
40      Continue
      End If
c===========================================================================c
c              read astrotidal  boundary condition                          c
c===========================================================================c
      If (NUMAST .NE. 0) Then
          Rewind IUT94
          Do I = 1, NUMAST
             Read (IUT94,5110) IABC(I), IABCINX(I), EMEAN(I)
             Read (IUT94,5100) (AMP(I,J), J = 1, IAST)
             Read (IUT94,5100) (PHASE(I,J), J = 1, IAST)
	    Enddo
50        Continue
      End If
C---------------------------------------------------------------------------C
C   S-A Read viscosity / onshore boundary condition at the first time step  c
c---------------------------------------------------------------------------c
	If (VERTMIX .EQ. 'SAMODEL   ' .OR. 
     *    DES .EQ. 'SADES  ' .OR. DES .EQ. 'SADDES ' .OR.
     *    DES .EQ. 'SAIDDES' .OR. DES .EQ. 'SAZDES ') Then
         If (NUMEBC .NE. 0) Then
            Rewind IUT96
            Do I = 1, 100000
		       Read (IUT96,5100,Err=170) T2VISE
               Read (IUT96,5100) (VISEBDRY(N,2), N = 1, NUMEBC)
               If (THOUR .LT. T2VISE) Go To 60
               T1VISE = T2VISE
               Do N = 1, NUMEBC
                  VISEBDRY(N,1) = VISEBDRY(N,2)
               Enddo
	      Enddo
60          Continue
         End If
         If (NUMQBC .NE. 0) Then
            Rewind IUT97
            Do I = 1, 100000
		     Read (IUT97,5100,Err=170) T2VISQ
	         Do N = 1, NUMQBC
                  Read (IUT97,5100) (VISQBDRY(N,2,K), K = 1, KBM)
	         Enddo
               If (THOUR .LT. T2VISQ) Go To 70
               T1VISQ = T2VISQ
               Do N = 1, NUMQBC
	            Do K = 1, KBM
                     VISQBDRY(N,1,K) = VISQBDRY(N,2,K)
	            Enddo
	         Enddo
	      Enddo
70	      Continue
	   Endif
	Endif
C----------------------------------------------------------------------------C
C SST k-w Read viscosity / onshore boundary condition at the first time step c
c----------------------------------------------------------------------------c
	If (VERTMIX .EQ. 'SSTMODEL  ' .OR. DES .EQ. 'SSTDES ') Then
         If (NUMEBC .NE. 0) Then
            Rewind IUT96
            Do I = 1, 100000
		       Read (IUT96,5100,Err=170) T2VISE
               Read (IUT96,5100) (TKEEBDRY(N,2), N = 1, NUMEBC)
               Read (IUT96,5100) (TDISSEBDRY(N,2), N = 1, NUMEBC)
               
               If (THOUR .LT. T2VISE) Go To 80
               T1VISE = T2VISE
               Do N = 1, NUMEBC
                  TKEEBDRY(N,1) = TKEEBDRY(N,2)
                  TDISSEBDRY(N,1) = TDISSEBDRY(N,2)                  
               Enddo
	      Enddo
80          Continue
         End If
         If (NUMQBC .NE. 0) Then
            Rewind IUT97
            Do I = 1, 100000
		       Read (IUT97,5100,Err=170) T2VISQ
	         Do N = 1, NUMQBC
                  Read (IUT97,5100) (TKEQBDRY(N,2,K), K = 1, KBM)
                  Read (IUT97,5100) (TDISSQBDRY(N,2,K), K = 1, KBM)
	         Enddo
               If (THOUR .LT. T2VISQ) Go To 90
               T1VISQ = T2VISQ
               Do N = 1, NUMQBC
	            Do K = 1, KBM
                     TKEQBDRY(N,1,K) = TKEQBDRY(N,2,K)
                     TDISSQBDRY(N,1,K) = TDISSQBDRY(N,2,K)
	            Enddo
	         Enddo
	      Enddo
90	      Continue
	   Endif
	Endif	
C----------------------------------------------------------------------------C
C Concentration Read boundary condition at the first time step               c
c----------------------------------------------------------------------------c
	If (WAQ .NE. 'NEGLECT') Then
         If (NUMEBC .NE. 0) Then
            Rewind IUT98
            Do I = 1, 100000
		       Read (IUT98,5100,Err=170) T2CEBC
               Read (IUT98,5100) (DCEBC(N,2), N = 1, NUMEBC)
               If (THOUR .LT. T2CEBC) Go To 100
               T1CEBC = T2CEBC
               Do N = 1, NUMEBC
                  DCEBC(N,1) = DCEBC(N,2)
               Enddo
	      Enddo
100         Continue
         End If
         If (NUMQBC .NE. 0) Then
            Rewind IUT99
            Do I = 1, 100000
		     Read (IUT99,5100,Err=170) T2CQBC
	         Do N = 1, NUMQBC
                  Read (IUT99,5100) (DCQBC(N,2,K), K = 1, KBM)
	         Enddo
               If (THOUR .LT. T2CQBC) Go To 101
               T1CQBC = T2CQBC
               Do N = 1, NUMQBC
	            Do K = 1, KBM
                     DCQBC(N,1,K) =DCQBC(N,2,K)
	            Enddo
	         Enddo
	      Enddo
101	      Continue
	   Endif
         If (NUMDBC .NE. 0) Then
            Rewind IUT101
            Do I = 1, 100000
   		       Read (IUT101,5100,Err=170) T2CDBC
               Read (IUT101,5100) (DCDBC(N,2), N = 1, NUMDBC) 
               If (THOUR .LT. T2CDBC) Go To 102
               T1CDBC = T2CDBC
               Do N = 1, NUMDBC
C	            Do K = 1, KBM
                  DCDBC(N,1) =DCDBC(N,2)
C	            Enddo
	         Enddo
	      Enddo
102	      Continue
	   Endif
	Endif
C============================================================================C
C                MASK FOR WAVE BREAKING                                      C
C============================================================================C
      Do I = 1, IJM
         WB_MASK(I) = 0.0
         Do K = 1, KBM
            VIS_BW(I,K) = 0.0
         Enddo
      Enddo	
C============================================================================C
      Return
170   Write (IUPRT,5000)
      Stop
5000  Format (//' there is insufficient temporal data for this run'/,
     &    '        revise input deck and resubmit '//)
5100  Format (8E14.7)
5110  Format ( 2I14, 8E14.7)

      End
