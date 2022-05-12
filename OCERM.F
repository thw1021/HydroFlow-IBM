C##########################################################################C
C                        MAIN PROGRAMM                                     C
C                   Version 1.0 (22/04/2009)                               C
C                 AUTHORIZED BY ZHANG JINGXIN                              C
C                               SHANGHAI JIAO TONG UNIVERSITY              C
C                               SHANGHAI, CHINA                            C
C--------------------------------------------------------------------------C         
C     NUMERICAL MODELLING OF ENVIRONMENTAL FLUID MECHANICS.                C
C     THE MODEL IS MAJORING IN THE TIDAL, WAVE, STORM, SEDIMENT AND        C
C     WATER QUALITY ISSUSES IN THE ESTUARY, COASTAL AND OCEAN.             C
C                                                                          C
c##########################################################################C
      Program OCERM
!$    Use OMP_LIB     
      Include './Include/OCERM_INF'
      Include './Include/VORGEN_INF'
      
      Dimension COM(80),TEMP1(IJP)
      Character*10 RESTAR,FN
      Logical Path_alive
      
!$    Call OMP_SET_NUM_THREADS(NTHREAD) 
      
C============================================================================C
C         Files index                                                        C
C============================================================================C
	ITF = 1
	ITSEDC = 2
	ITSEDB = 3
	ITTEM = 4
	ITSAL = 5
	ITC = 7
      
      IURUN = 10
      IUGRD = 11
	IUCUV = 12
      IUPRT = 13
      IUPRT_DEM = 131      
      IUHYD = 14
      IUPLT = 15
      IUWRS = 16
      IURRS = 17
      IUFLW = 18
      IUTSR = 19
      IUSTA = 20
      IUTSR2 = 21

      IDEM= 22
      IDEM_BIN = 25
      IDEM2= 23     
      
      IBMSET = 60
      IBMSP = 61
      IBM2D = 62
      IBM3D = 63
      IBMDATA = 66
      IBMB = 70
      IBMF = 71
      IBMS = 72
      IBMGC = 73
      IBMGCI = 74
      IBMP = 75
      IBMN = 76
      IBMT = 77
      IBMDIB = 78
      IBMDF = 79
      IBMCP = 80
      IBMCDCL = 81
      
      IUT90 = 90
      IUT91 = 91
      IUT92 = 92
      IUT93 = 93
      IUT94 = 94
      IUT95 = 95
      IUT96 = 96
      IUT97 = 97
      IUT98 = 98
      IUT99 = 99
      IUT101 =101 

      IDEMRRS = 102                         ! HOT START FOR DEM
      IDEMSTA = 103                         ! HOT START FOR DEM
      IDEMWRS = 104                         ! HOT START FOR DEM
      
      IAVEDAT  = 201
      IAVESTA  = 202
      IAVERES  = 203
      ITF_BIN  = 204
      ITC_BIN  = 205
      IUINS    = 612
      IVINS    = 613
      IWINS    = 614
	ITEMP = 1000

C============================================================================C
C                                                                            C
C     *_run: the input file which includes input data for time step,         C 
C          forcings, and boundary points and conditions;                     C
C     *_grid: the input file which includes the vertical and horizontal      C
C           grids.                                                           C
C     *_cuv: geometry information of then meshes                             C
C                                                                            C
C============================================================================C
    
      Inquire(directory="./Include/",exist = Path_alive)
      If (.not. Path_alive) call system(trim(adjustl("mkdir Include")))
      Inquire(directory="./Records/",exist = Path_alive)
      If (.not. Path_alive) call system(trim(adjustl("mkdir Records")))
      Inquire(directory="./Restart/",exist = Path_alive)
      If (.not. Path_alive) call system(trim(adjustl("mkdir Restart")))
      Inquire(directory="./Results/",exist = Path_alive)
      If (.not. Path_alive) call system(trim(adjustl("mkdir Results")))
      Inquire(directory="./Results/IBM",exist = Path_alive)
      If (.not. Path_alive) 
     &    call system(trim(adjustl("mkdir Results\IBM")))

      Open (IURUN, File = './Include/OCERM.RUN', status='unknown')
      Open (IUGRD, File = './Include/OCERM.GRD', status='unknown')
      Open (IUCUV, File= './Include/OCERM.CUV', status='unknown')
      
C============================================================================C
C                                                                            C
C     *_output_prt: the output used to check the input data and trace the    C
C                  model operations. The interval of output can be           C
C                  specified using the parameter;                            C
C     *_output_plt: the output of velocities, temperature, salinity, etc     C
C                  averaged over a specified period;                         C
C     *_output_tsr: the output of the time series of velocities, temperature C
C                  , salinity, etc over a specified interval;                C
C     *_output_flow: the output file used to the water quality model         C
C                                                                            C
C============================================================================C

      Open (IUPRT,File='./Records/OCERM.PRT',status='unknown')
      Open (IUPRT_DEM,File='./Records/OCERM_DEM.PRT',status='unknown')
      Open (IUPLT,File='./Records/OCERM.PLT',status='unknown')
      Open (IUTSR,File='./Results/OBS.TSR',status='unknown')
      Open (IUTSR2,File='./Results/FLUX.TSR',status='unknown')
      Open (IUFLW,File='./Records/OCERM.FLW',Form='unformatted'
     &         ,status='unknown')
      
C============================================================================C
C                                                                            C
C     *_restart: the input data used to restart the model run                C 
C     *_startup: the output data used to restart the model run.              C
C                  it needs to be copied to *_restart before re-             C 
C                  starting the model again.                                 C
c                                                                            C
C============================================================================C

      Open (IURRS,Form='unformatted',File='./Restart/OCERM.RES'
     &         ,status='unknown')
      Open (IUWRS,Form='unformatted',File='./Restart/OCERM.STA'
     &       ,status='unknown')
      
C============================================================================C
C                                                                            C
C     *_temp*: the temporal files used to storage the data during the        C
C             model run.                                                     C
C                                                                            C
C============================================================================C

      Open (IUT90,File='./Records/temp.t90',status='unknown')
      Open (IUT91,File='./Records/temp.t91',status='unknown')
      Open (IUT92,File='./Records/temp.t92',status='unknown')
      Open (IUT93,File='./Records/temp.t93',status='unknown')
      Open (IUT94,File='./Records/temp.t94',status='unknown')
      Open (IUT95,File='./Records/temp.t95',status='unknown')
      Open (IUT96,File='./Records/temp.t96',status='unknown')
      Open (IUT97,File='./Records/temp.t97',status='unknown')
      Open (IUT98,File='./Records/temp.t98',status='unknown')
      Open (IUT99,File='./Records/temp.t99',status='unknown')
      Open (IUT101,File='./Records/temp.101',status='unknown')
      
C============================================================================C
C                                                                            C
C     *_tram: the output simulated data at every mesh point over a           C
C             specified interval                                             C
C                                                                            C
C============================================================================C                      
     
      If (ISAVE_BIN==0) Then
	  Open (ITF, File = './Results/VEL_TRAM.DAT',status='UNKNOWN')
	  Write (ITF,*) 'TITLE="Shallow water equation solution" '
	  Write (ITF,*) 'VARIBLES= "X","Y","Z","U","V","W","ELF" '        
      Else if (ISAVE_BIN==1) Then
        Open (ITF_BIN,Form='unformatted',File ='./Results/VEL_TRAM.BIN')
      Endif
      
      If (IVEL_INS.NE.0) Then
        Open(IUINS,File='./Results/U-INS.DAT',status='UNKNOWN')  
        Open(IVINS,File='./Results/V-INS.DAT',status='UNKNOWN')
        Open(IWINS,File='./Results/W-INS.DAT',status='UNKNOWN')
      Endif
      
      If (IAVESAVE.NE.0) Then   
        Open(IAVEDAT, File = './Results/VEL_AVE.DAT',status='UNKNOWN')
        Write(IAVEDAT,*) 'TITLE="results of statistics" '
        Write(IAVEDAT,'(A90)')'VARIBLES= "X","Y","Z","AVEU","AVEV",
     &        "AVEW","AVEUT","Ruu","Ruv","Ruw","Rvv","Rvw","Rww" '
      Endif  
      
C============================================================================C
C         Read the *_run data and *_grid data                                C
C============================================================================C
      
      Read (IURUN,5000) (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
	Read (IURUN,5000) (COM(I),I = 1,80)
	Write (IUPRT,5100) (COM(I),I = 1,80)

	Read (IURUN,5050) HYDTYPE, WIND, WAVEDYN, SEDTRAN, TFIELD, THINDAM
     &            , SFIELD, WAQ, CORLIS, DES, DEM, POREMODULE, MODULEIBM
	Write (IUPRT,5060) HYDTYPE, WIND, WAVEDYN, SEDTRAN, TFIELD,THINDAM
     &            , SFIELD, WAQ, CORLIS, DES, DEM, POREMODULE, MODULEIBM
	Read (IURUN,5000) (COM(I),I = 1,80)
	Write (IUPRT,5100) (COM(I),I = 1,80)
	Read (IURUN,5070)  IYR, IMO, IDA, IHR
	Write (IUPRT,5080) IYR, IMO, IDA, IHR  
      Read (IURUN,5000) (COM(I),I = 1,80)
      Read (IURUN,5200) DTI, IRAMP, IHOTSTART
      DAYI = 1. / 86400.
      GRAV = 9.806

      Write (IUPRT,5300) DTI, IRAMP, IHOTSTART

      If (ISAVE_BIN==0 .AND. DEM .NE. 'NEGLECT') Then
        Open (IDEM,FILE= './Results/DEM_PARTICLE.DAT',STATUS='UNKNOWN')
        Write (IDEM,*) 'TITLE="DEM results of particles" '
        Write (IDEM,*) 'VARIBLES= "X","Y","Z","U","V","W" '         
      Endif      
      
      If (WAQ .NE. 'NEGLECT') Then
        If (ISAVE_BIN==0) Then  
        Open (ITC, FILE = './Results/C_TRAM.DAT',STATUS='UNKNOWN')
	  Write(ITC,*)'TITLE="Resolved mass transport equation solution" '
        Write (ITC,*) 'VARIBLES= "X","Y","Z","C" ' 
        Elseif (ISAVE_BIN==1) Then
        Open (ITC_BIN,Form='unformatted',FILE = './Results/C_TRAM.BIN')
        Endif    
      Endif
      
C============================================================================C
C                                                                            C
C     Type of run                                                            C
C     Prognostic: calculation with t and s                                   C
C     Diagnostic: calculation with t and s held fixed                        C
C     3-d:  type of momentum advection and bottom friction                   C
C     Linear    ; all momentum advection neglected                           C 
C     Non-linear; complete physics                                           C
C                                                                            C
C============================================================================C

      Read (IURUN,5000) (COM(I),I = 1,80)
      Read (IURUN,5400) NSTEPS, IPRINT, RESTAR, TOR, ADVECT
      Read (IURUN,5000) (COM(I),I = 1,80)
      Read (IURUN,5600) BFRIC, Z0B, WFBC
      Read (IURUN,5000) (COM(I),I = 1,80)
      Read (IURUN,5500) SGSTYPE, UNICON, INDEX_SGS
      Read (IURUN,5000) (COM(I),I = 1,80)
      Read (IURUN,5500) HORZMIX, HORCON, HPRNU
      Read (IURUN,5000) (COM(I),I = 1,80)
      Read (IURUN,5500) VERTMIX, VERCON, UMOL, VPRNU
      
C============================================================================C
C         Read history output times                                          C
C============================================================================C
      
      Read (IURUN,5000) (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
      Read (IURUN,5900) JHM, IAVGE
      Write (IUPRT,5900) JHM, IAVGE
      If (JHM.NE.0) Then
        DEI = 1. / Float(IAVGE)
        Read (IURUN,5800) (IHIST(I,2),I = 1,JHM)
        Write (IUPRT,5800) (IHIST(I,2),I = 1,JHM)
      Endif
      
C============================================================================C
C         Read  averaging interval                                           C
C============================================================================C
      
      Read (IURUN,5000) (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
      Read (IURUN,5900) ISKILL, ISAVE
      Write (IUPRT,5900) ISKILL,ISAVE
      If (ISKILL.EQ.0) Then
         SKILLI = 1.0
      Else
         SKILLI = 1. / Float(ISKILL)
      Endif
      
C============================================================================C
C         Define of the output of elevation and velocity                     C
C============================================================================C
      
      Read (IURUN,5000) (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
      Read (IURUN,5700) EPTS
      Write (IUPRT,5700) EPTS
      If (EPTS .NE. 0) Then
          Read (IURUN,5700) (INXE(I),I = 1,EPTS)
          Write (IUPRT,5700) (INXE(I),I = 1,EPTS)
      End If
      
      Read (IURUN,5000) (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
      Read (IURUN,5700) VPTS
      Write (IUPRT,5700) VPTS
      If (VPTS .NE. 0) Then
          Read (IURUN,5700) (INXV(I),I = 1, VPTS)
          Write (IUPRT,5700) (INXV(I),I = 1, VPTS)
      End If
      
      Read (IURUN,5000) (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
      Read (IURUN,5700) FPTS
      Write (IUPRT,5700) FPTS
      If (FPTS .NE. 0) Then
      Do N = 1, FPTS
          Read (IURUN,5700) ISFLX(N), JSFLX(N), NFLXE(N)
          Write (IUPRT,5700) ISFLX(N),JSFLX(N), NFLXE(N)
      Enddo
      Endif
      
C============================================================================C
C         Read defining of elevation boundary                                C
C============================================================================C
      
      Read (IURUN,5000) (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
      Read (IURUN,5700) NUMEBC
      Write (IUPRT,5700) NUMEBC
	If(NUMEBC .NE. 0) Then
          Read (IURUN,5750) NAMEEBC
          Write (IUPRT,5750) NAMEEBC
      Endif
      
C============================================================================C
C         Read defining of velocity boundary                                 C
C============================================================================C
      
      Read (IURUN,5000) (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
      Read (IURUN,5700) NUMVBC
      Write (IUPRT,5700) NUMVBC
	If(NUMVBC .NE. 0) Then
          Read (IURUN,5750) NAMEVBC
          Write (IUPRT,5750) NAMEVBC
      Endif
      
C============================================================================C
C         Read defining of discharge boundary                                C
C============================================================================C
      
      Read (IURUN,5000) (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
      Read (IURUN,5700) NUMQBC
      Write (IUPRT,5700) NUMQBC
	If(NUMQBC .NE. 0) Then
          Read (IURUN,5750) NAMEQBC
          Write (IUPRT,5750) NAMEQBC
      Endif
      
C============================================================================C
C         Read defining of offshore discharge                                C
C============================================================================C
      
      Read (IURUN,5000) (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
      Read (IURUN,5700) NUMDBC
      Write (IUPRT,5700) NUMDBC
	If(NUMDBC .NE. 0) Then
          Read (IURUN,5750) NAMEDBC
          Write (IUPRT,5750) NAMEDBC
      Endif
      
C============================================================================C
C         Read defining of astrotide boundary                                C
C============================================================================C
      
      Read (IURUN,5000) (COM(I),I = 1,80)
      Write (IUPRT,5100) (COM(I),I = 1,80)
      Read (IURUN,5700) NUMAST
      Write (IUPRT,5700) NUMAST
	If(NUMAST .NE. 0) Then
          Read (IURUN,5750) NAMEAST
          Write (IUPRT,5750) NAMEAST
      Endif
      
C============================================================================C
C         Result files of immersed boundary method                           C
C============================================================================C
      
      Open(IBMDATA,File = './Results/IBM/ibmdata.dat', status='unknown')
      Open(IBMB,File = './Results/IBM/ibn.dat', status='unknown')
      Open(IBMF,File = './Results/IBM/fluid.dat', status='unknown')
      Open(IBMS,File = './Results/IBM/solid.dat', status='unknown')
      Open(IBMGC,File = './Results/IBM/gc.dat', status='unknown')
      Open(IBMGCI,File = './Results/IBM/gcimage.dat', status='unknown')
      Open(IBMP,File = './Results/IBM/points.dat', status='unknown')
      Open(IBMN,File = './Results/IBM/normal.dat', status='unknown')
      Open(IBMT,File = './Results/IBM/tangent.dat', status='unknown')
      Open(IBMDIB,File = './Results/IBM/delta_ib.dat', status='unknown')
      Open(IBMDF,File = './Results/IBM/delta_f.dat', status='unknown')
      Open(IBMCP,File = './Results/IBM/Cp.dat', status='unknown')
      Open(IBMCDCL,File = './Results/IBM/CDCL.dat', status='unknown') 
      
C=============================================================================C
C         Mode of the simulation                                              C
C=============================================================================C
      
      If (TOR .NE. 'BAROTROPIC' .AND. TOR .NE. 'BAROCLINIC') Then
        Write (IUPRT,6100) TOR
        Stop
      Endif

      If (ADVECT .NE. 'LINEAR    ' .AND. ADVECT.NE.'NON-LINEAR') Then
        Write (IUPRT,6200) ADVECT
        Stop
      Endif

      If (HORZMIX .NE. 'CLOSURE   ' .AND. HORZMIX.NE.'CONSTANT  ') Then
        Write (IUPRT,6300) HORZMIX
        Stop
      Endif

      If (VERTMIX .NE. 'CLOSURE   ' .AND. VERTMIX .NE. 'CONSTANT  '
     &	.AND. VERTMIX .NE. 'SMAGMODEL ' .AND. 
     &    VERTMIX .NE. 'SAMODEL   '. AND.VERTMIX .NE. 'SSTMODEL  ') Then
        Write (IUPRT,6400) VERTMIX
        Stop
      Endif
      Write (IUPRT,6500) TOR
      Write (IUPRT,6600) ADVECT

      If (HORZMIX .EQ. 'CLOSURE   ') Then
        Write (IUPRT,6700) HORZMIX, HORCON, HPRNU
      Else
        Write (IUPRT,6800) HORZMIX, HORCON, HPRNU
      Endif

      If (VERTMIX .EQ. 'CLOSURE   ') Then
        Write (IUPRT,6900) VERTMIX, VERCON, UMOL, VPRNU
      Else
        Write (IUPRT,7000) VERTMIX, VERCON, UMOL, VPRNU
      Endif
      
C============================================================================C
C         Mesh information and boundary conditions setting                   C 
C============================================================================C
      
	HOURSTAR = 0.0
	If (RESTAR .EQ. 'COLD START') Then
          Call SETDOM                         !  domain information
          Call BCDATA                         !  input the boundary conditions
          Call ZEROES                         !  specifying the initial values
      Else
	    Call SETDOM
          Read (IURRS) NSTEP, DZR, Z, ZZ, DZ, DZZ, HC, HS, HP, DC, DS, 
     &        EL, ELFV, ELFM, COR, WTSURF, WTBOT, AAM, KM, KH, DT,
     &        KQ, U, V, W, UR, VR, WR, UV, VV, WV, UN, VN, CBC, VIS,
     &        T, SED, RMEAN, Z01, PN, QZ, CCM, CFM, INDEX_EDGE,
     &        HOURSTAR, TBX, TBY, TKE, TDISS, STRESS,
     &        VIS, TENSOR, VORTICITY,C,QZ,GRADP
          Close (IURRS)

	    Call BCDATA 
	    If(HYDTYPE .EQ. 'NONSTATIC') Then
              Do I = 1, IJM
              If(CCM(I) .EQ. 1.0) Then
	            Do K = 1, KBM
	               QZ(I,K) = WR(I,K) * DC(I)
	            Enddo
              Endif
          Enddo
          Endif         
      Endif

      If(DEM .NE. 'NEGLECT' .AND. DEMRESTART .NE. 'COLD START') Then    
          Open (IDEMRRS,Form='unformatted',File='./Restart/DEM.RES'
     &         ,status='unknown')
          
          Read (IDEMRRS) NNDEM,TTTDEM,XPDEM,YPDEM,ZPDEM,
     &        UPDEM,VPDEM,WPDEM,OMXDEM,OMYDEM,OMZDEM,
     &        PORE,DELTAT_PP1X,DELTAT_PP1Y,DELTAT_PP1Z,
     &        N_P2,N_PP2,UPSUMDEM,VPSUMDEM,WPSUMDEM,
     &        UFPDEM_PRE,VFPDEM_PRE,WFPDEM_PRE,
     &        UFDEM_PRE,VFDEM_PRE,WFDEM_PRE,
     &        N_SYNFA,DEL_HIS,NSYM_CONT
          Close (IDEMRRS)
        
          Call SEARCH
          Call IPPKPPINITIAL   
        
      Endif

C============================================================================C
C         Bed frictional coefficients input                                  C 
C============================================================================C

      If(Z0B .EQ. 1.0) Then
	    Open(ITEMP, File= '.\Include\ROUGH.GRD',Status='UNKNOWN')
	    Do K = 1, IJP
	        Read (ITEMP,*) PXY(K,1), PXY(K,2), TEMP1(K)
	    Enddo
	    Do I = 1, IJM
	        If(CELL_POLYGEN(I) .EQ. 4) Then
	        Z01(I) = (TEMP1(CELL_SIDE(I,1,3))+TEMP1(CELL_SIDE(I,2,3))+
     &		      TEMP1(CELL_SIDE(I,3,3))+TEMP1(CELL_SIDE(I,4,3)))/4.
	        Else
	        Z01(I) = (TEMP1(CELL_SIDE(I,1,3))+TEMP1(CELL_SIDE(I,2,3))+ 
     &		      TEMP1(CELL_SIDE(I,3,3))) / 3.
	        Endif
	   Enddo 
	   Close(ITEMP)
	Else
	    Do IJ = 1, IJM
	        Z01(IJ) = Z0B
	    Enddo
	Endif

      ISTART = NSTEP + 1
      TNDAYS = Float(NSTEPS) * DAYI * DTI
      IEND = NSTEP + NSTEPS
      TPRT = Float(IPRINT) * DAYI * DTI
      TAVG = Float(IAVGE) * DAYI * DTI
      TSKILL = Float(ISKILL) * DAYI * DTI

      Write (IUPRT,7100) ISTART,IEND
      Write (IUPRT,7200) TNDAYS
      Write (IUPRT,7300) TPRT,IPRINT,TAVG,IAVGE,TSKILL,ISKILL
      Write (IUPRT,7400)
      Write (IUPRT,7800) BFRIC, Z0B, WFBC
      Write (IUPRT,7900)
      Do K = 1, KBM
          Write (IUPRT,8000) K, Z(K), ZZ(K), DZ(K)
      Enddo
      Write (IUPRT,8000) K, Z(KB)
      THOUR = HOURSTAR
      
      Call FIRST 

c	Call Random_seed()
!$    begin1 = OMP_GET_WTIME()
c      call cpu_time(begin1)
 
      If (IAVESAVE==1 .AND. IAVEHOT==1) Then
           Open(IAVERES,Form='unformatted',File='./Restart/AVE.RES',
     &	      Status='unknown')
	     Read (IAVERES)NSTEPAVE, AVEURV, AVEVRV,AVEWRV, AVEUTRV,
     &                  RYNUU, RYNUV, RYNUW, RYNVV, RYNVW,RYNWW  
	     Close(IAVERES)          
      Endif
                     
C**************************************************************************
C                                                                         *
C             Begin numerical simulation                                  *
C                                                                         *
C**************************************************************************
                     
      Do 100 NSTEP = ISTART, IEND  
      Print*, 'NUM. step of numerical case is ', NSTEP, '   /',IEND
      
C*************************************************************************C
          
      ITENPCT = (IEND - ISTART+1)/10
      RAMP = Tanh(Float(NSTEP)/Float(IRAMP+1)) / Tanh(1.0)
      If(NSTEP.GT.IRAMP) RAMP=1.0
      TIME = Float(NSTEP-ISTART+1) * DAYI * DTI
      THOUR = TIME * 24. + HOURSTAR
      
      Call BCOND(10)                              !  boundary conditions

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
      Do I = 1, IJM
          ELF(I) = EL(I)
      Enddo
!$OMP END PARALLEL DO

C============================================================================C
C         DEM module
C============================================================================C

      If (DEM .NE. 'NEGLECT') Then
          Call DEMM             
      Endif
      Call PORECAL  
                
C============================================================================C
C         Computing wind stress on the surface
C============================================================================C
        
      If (WIND.NE.'NEGLECT') Then
c         Call WINDUV                              !  surface wind stress
      Endif

C============================================================================C
C         Computing baroclinic pressure gradient
C============================================================================C

c	   Call BAROCLINIC                           !  baroclinic pressure

C============================================================================C
C         Computing viscosity coefficients 
C============================================================================C

      If (HORZMIX .EQ. 'CLOSURE   ') Call SUBGRIDH   !  horizontal edd vis.
      If (VERTMIX .NE. 'CONSTANT  ') Call SUBGRIDV   !  vertical edd vis.
      If (SGSTYPE .EQ. 'INCLUDE   ') Call SGSMODEL   !  LES modle 
      If (DES .EQ. 'SADES  ' .OR.
     &    DES .EQ. 'SADDES ' .OR.
     &    DES .EQ. 'SAIDDES' .OR.
     &    DES .EQ. 'SAZDES ' )  Call DESSA   !  DES modle
      If (DES .EQ. 'SSTDES ')   Call DESSST          !  DES modle

      If(WAVE_BREAKING .EQ. 1.) Then
      Do I = 1, IJM
          If(CCM(I) .EQ. 1.0) Then
          Do K = 1, KBM
              AAM(I,K) = AAM(I,K) + VIS_BW(I,K)
              KM(I,K) = KM(I,K) + VIS_BW(I,K) 
          Enddo
          Endif
      Enddo
      Endif   
      
C============================================================================C
C         Computing bed friction coefficient 
C============================================================================C
      
      If(WAVEDYN .EQ. 'NEGLECT') Then
          If(WFBC .EQ. '      FUN1') Then
          If (VERTMIX .EQ. 'SSTMODEL  ') Then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(IJ,ZSTAR)  
          Do IJ = 1, IJM
              If (DC(IJ) .GT. 0.0) Then
              If (Sqrt(UR(IJ,KBM)**2.+VR(IJ,KBM)**2.).LT.1.0E-10 ) Then
                  CBC(IJ)=(0.41/Log((ZZ(KBM)-Z(KB))*DC(IJ)/Z01(IJ)))**2.
              Else  
                  ZSTAR = Sqrt(0.3)*
     &                    Sqrt(TKE(IJ,KBM)+1.E-20)*DC(IJ)*DZZ(KBM)/1.E-6
                  CBC(IJ) = 0.41*Sqrt(0.3)*
     &                      Sqrt(TKE(IJ,KBM)+1.E-20)/Log(9.81*ZSTAR)
              Endif    
              Endif    
          Enddo
!$OMP END PARALLEL DO                  
          Else    
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(IJ,Z0,CBCMIN)
          Do IJ = 1, IJM
              If (DC(IJ) .GT. 0.0) Then
                  Z0 = Z0B
                  CBCMIN = BFRIC
                  CBC(IJ)=(0.41/Log((ZZ(KBM)-Z(KB))*DC(IJ)/Z01(IJ)))**2.
! manning coefficients 	
c			   CBC(IJ) = 9.8 * Z01(IJ) **2. / DC(IJ) ** 0.333
c	           CBC(IJ)=Dmax1(CBCMIN,CBC(IJ))
              Endif
          Enddo
!$OMP END PARALLEL DO
          Endif    
          Endif  
c         Else
c             If(WFBC .EQ. '      FUN2') Call WALLDRAG 
c         Endif
      Endif
      
C============================================================================C
C         Feeding velocity fluctuation for DES                             
C============================================================================C
      
      If(DES .EQ. 'SAZDES ') Then 
          TIME1 = NSTEP * DTI
          Print*, des
          Stop
          Call TurGenerator(TIME1,NSTEPS)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(J,K)           
          Do J = 1, NUM_CELL
          Do K = 1, NUM_VER
              UDIS(J,K) = U_FLU(J,K) * R11(J,K)
              WDIS(J,K) = U_FLU(J,K) * R21(J,K) + V_FLU(J,K) * R22(J,K)
              VDIS(J,K) = U_FLU(J,K) * R31(J,K) + V_FLU(J,K) * R32(J,K) 
     &                  + W_FLU(J,K) * R33(J,K) 
          Enddo
          Enddo
!$OMP END PARALLEL DO
      Endif
      
C============================================================================C
C         Velocities and elevation solver
C============================================================================C

      Call ADVU                                !  momentum in x direction
      Call ADVV                                !  momentum in y direction
      If (HYDTYPE .EQ. 'NONSTATIC') Call ADVW  !  momentum in z direction 
c      If (TFIELD .NE. 'NEGLECT')  Call ADVT    !  momentum of tem.
c      If (SFIELD .NE. 'NEGLECT')  Call ADVSAL  !  momentum of sal.
c      If (SEDTRAN .NE. 'NEGLECT') Call ADVSED  !  momentum of sed.
      If (WAQ .NE. 'NEGLECT')     Call ADVC    !  momentum of mass.
      
      Call ELTION                              !  compute the elevation
      Call PROFV                               !  compute the velocities     
c      If (TFIELD .NE. 'NEGLECT')  Call PROFT   !  momentum of tem.
c      If (SFIELD .NE. 'NEGLECT')  Call PROFSAL !  momentum of sal.
c      If (SEDTRAN .NE. 'NEGLECT') Call PROFSED !  momentum of sed.
      If (WAQ .NE. 'NEGLECT')     Call PROFC   !  momentum of mass.
      
C============================================================================C
C         Hydrostatic and nodhydrodynamic mode
C============================================================================C
      
      If (HYDTYPE .NE. 'NONSTATIC') Then
          Call VERTVL                           !  compute qw  
          Call WREAL                            !  compute w 
      Else
          Call PROFW                            !  compute the DW
          Call DYN                              !  hydrodynamic model
          Call UPDATEFLOW                       !  update after DYN
          Call WREAL
      Endif
      
C============================================================================C
C         Immersed boundary method                                          
C============================================================================C
      
      If (MODULEIBM .NE. 'NEGLECT') Then
          Call IBM
      Endif   
      
C============================================================================C
C         Value updating 
C============================================================================C
      
      Call REUV                                !  compute u v 
      Call UVFN                                !  vel. on sides of cells
      
      Do I = 1, IJM
          EL(I) = ELF(I)
      Enddo

C============================================================================C
C         Divergence chekcing                                                
C============================================================================C

      UVMAX = 0
      IJMAX = 0
      Do I = 1, IJM
      If (CCM(I) .EQ. 1.0) Then
          UVTEMP = Sqrt(UAVE(I)** 2. + VAVE(I)** 2.)
          If(UVTEMP .GT. UVMAX) Then
              UVMAX = UVTEMP
              IJMAX = I
          Endif
          If(UVMAX .GT. 100.0) Goto 101
          
          IJNAN = 1             
          If (ISNAN(ELF(I))) Then
              IJNAN = I
              goto 101
          Endif
      Endif
      Enddo
	  
C============================================================================C
C         Dats saving                                                        
C============================================================================C

      Call ARCHIVE(DAYI,HOURSTAR)                !  save in disk
      Call IBMARCHIVE                        !  save data of IBM
      
      If (IAVESAVE == 1) Then
          If (IAVEHOT == 0) Then
              Call STATISTICS(DAYI,HOURSTAR)
          Else if (IAVEHOT==1 .AND. NSTEP .GT. NSTEPAVE) Then      
              Call STATISTICS(DAYI,HOURSTAR)  
          Endif  
      Endif
 
C*************************************************************************C
 
100   Continue
          
C**************************************************************************
C                                                                         *
C                 End numerical simulation                                *
C                                                                         *
C**************************************************************************

1001  Write(*,*) 'Numerical case is successfully completed !'
!$    END1 = OMP_GET_WTIME()
C      call cpu_time(end1)
      Print*, 'Total time of numrical case is', END1-BEGIN1
      
C============================================================================C
C         Boundary condition iteration module                                                
C============================================================================C
      
      I_iteration = 0
      If (I_iteration.EQ.1) Call BCITERATION 
      
C============================================================================C
C         Divergence analysis                                                   
C============================================================================C
      
101   Continue
	If(UVMAX .GT. 100.0) THEN
          Write(IUPRT,8400) TIME,IJMAX,UVMAX,CXY(IJMAX,1),CXY(IJMAX,2)
          Write (*,*) 'Numerical case has been diverged !'
          Call SAVEDIVERGE
      Else if (ISNAN(ELF(IJNAN))) Then
          Write (*,*) 'The water level of number :',IJNAN,'is NaN !'
          Write (*,*) 'Numerical case has been diverged !'
	Else
          Write (IUPRT,8300) TIME
      Endif
      
      Write (IUWRS) NSTEP-1, DZR, Z, ZZ, DZ, DZZ, HC, HS, HP, DC, DS, 
     &        EL, ELFV, ELFM, COR, WTSURF, WTBOT, AAM, KM, KH, DT,
     &        KQ, U, V, W, UR, VR, WR, UV, VV, WV, UN, VN, CBC, VIS,
     &        T, SED, RMEAN, Z01, PN, QZ, CCM, CFM, INDEX_EDGE,
     &        THOUR, TBX, TBY, TKE, TDISS, STRESS,
     &        VIS, TENSOR, VORTICITY,C,QZ,GRADP
C     *        NNDEM,TTTDEM,XPDEM,YPDEM,ZPDEM,
C     *        UPDEM,VPDEM,WPDEM,OMXDEM,OMYDEM,OMZDEM,
C     *        PORE,GRADP,
C     *        DELTAT_PP1X,DELTAT_PP1Y,DELTAT_PP1Z,N_P2,N_PP2 
      If(DEM .NE. 'NEGLECT') Then  
          Open (IDEMWRS,Form='unformatted',File='./Restart/DEM.STA'
     &       ,status='unknown')          
          
          Write (IDEMWRS) NNDEM,TTTDEM,XPDEM,YPDEM,ZPDEM,
     &        UPDEM,VPDEM,WPDEM,OMXDEM,OMYDEM,OMZDEM,
     &        PORE,DELTAT_PP1X,DELTAT_PP1Y,DELTAT_PP1Z,
     &        N_P2,N_PP2,UPSUMDEM,VPSUMDEM,WPSUMDEM,
     &        UFPDEM_PRE,VFPDEM_PRE,WFPDEM_PRE,
     &        UFDEM_PRE,VFDEM_PRE,WFDEM_PRE,
     &        N_SYNFA,DEL_HIS,NSYM_CONT 
          Close(IDEMWRS) 
      Endif
      
C============================================================================C
C         Close data files                                                  
C============================================================================C  
        
      Close(IUWRS)
      Close(IUT90)
      Close(IUT91)
      Close(IUT92)
      Close(IUT93)
      Close(IUT94)
      Close(101)      
      Close(IUINS)
      Close(IVINS)
      Close(IWINS)
      Close(IAVEDAT)
      Close(IUPRT_DEM)
      
      IF (ISAVE_BIN==0) THEN
          Close(ITF)
          If (WAQ .NE. 'NEGLECT') Close(ITC)
      Else if (ISAVE_BIN==1) Then
          Close(ITF_BIN)
          If (WAQ .NE. 'NEGLECT') Close(ITC_BIN)    
      Endif

      If (DEM .NE. 'NEGLECT') Then
          Close(IDEM)
c        Close(IDEM2)
      Endif  
      
      Close(IBMDATA)
      Close(IBMB)
      Close(IBMF)
      Close(IBMS)
      Close(IBMGC)
      Close(IBMGCI)
      Close(IBMP)
      Close(IBMN)
      Close(IBMT)
      Close(IBMDIB)
      Close(IBMDF)
      Close(IBMCP)
      Close(IBMCDCL)
      
C============================================================================C
C         Format setting                                                      
C============================================================================C      
      
5000  Format (100A1)
5050  Format (A9,12(2X,A7))
5060  Format (' the hydrodynamic type is     ',A10,//,'wind is   ',A7,
     &    //,'wave is   ',A7,
     &    //,'sedment is   ',A7//,'temperature is   ',A7,
     &    //,'thin dam defined is  ',A7,
     &    //,'salinity defined is  ',A7,
     &    //,'water quality defined is  ',A7,
     &    //,'corlis considered is  ',A7,
     &    //,'des considered is  ',A7, 
     &    //,'demmodule considered is  ',A7, 
     &    //,'poremodule considered is  ',A7//)
5070  Format (4I6)
5080  Format ('the primary simulation time is      ',
     &    I6,'.',I6,'.',I6,'.',I6)
5100  Format (/1X,80A1/)
5200  Format (1E10.3,I10,I10)
5300  Format (' the model time step is                  ',F10.4,
     &    ' seconds',//,' number of ramp time steps               ',I10
     &    //,' number of temporal saving  ', I10//)
5400  Format (2I10,1X,A10,1X,A10,1X,A10)
5500  Format (A10,3E10.3)
5600  Format (2E10.3,A10)
5700  Format (10I8)
5750  Format (A8)
5800  Format (10I8)
5900  Format (8I10)
6000  Format (6(2I5,1X,A4,I5))
6100  Format (//' type of run (tor=',A10,') is specified incorrectly! ',
     &    'fix and resubmit'//)
6200  Format (//'  type of advection (advect=',A10,') is specified ',
     &    'incorrectly!    fix and resubmit'//)
6300  Format (//'  type of horizontal mixing (horzmix=',A10,
     &    ') is specified incorrectly!    fix and resubmit'//)
6400  Format (//'  type of vertical mixing (vertmix=',A10,
     &    ') is specified incorrectly!    fix and resubmit'//)
6500  Format (/' this omega model run is',2X,A10/)
6600  Format (/' this simulation has ',A10,' momentun advection '/)
6700  Format (/' this simulation has ',A10,' horizontal mixing ',
     &    ' horcon = ',1PE10.3,' hprnu = ',1PE10.3/)
6800  Format (/' this simulation has ',A10,' horizontal mixing ',
     &    ' constant = ',1PE10.3,'m**2/s  hprnu = ',1PE10.3/)
6900  Format (/' this simulation has ',A10,' vertical mixing ',
     &    'vercon = ',1pe10.3,' umol = ',1PE10.3,
     &    ' vprnu = ',1PE10.3/)
7000  Format (/' this simulation has ',A10,' vertical mixing ',
     &    'vercon = ',1pe10.3,' constant = ',1PE10.3,
     &    'm**2/s  vprnu = ',1PE10.3/)
7100  Format (//' model starting up...istart = ',I6,' iend = ',I6/)
7200  Format (//' number of days in simulation = ',F6.2/)
7300  Format (//' tprt =   ',F10.3,'  iprint =  ',I10,//' tavg =   ',F10
     &    .3,'  iavge =   ',I10,//' tskill = ',F10.3,'  iskill =  ',I10
     &    //)
7400  Format (//' history taken at timesteps  start---stop ')
7500  Format (28X,I5,2X,I5)
7600  Format (//' quality parameters for aesop integrated over timesteps
     &      start      stop')
7700  Format (57X,I8,2X,I8)
7800  Format (//' bfric            =   ',F10.4,' nondimensional'/
     &    ' z0b              =   ',F10.4,' m'/' WALL TYPE  =   ',
     &    A10,' seconds'/)
7900  Format (/1X,' k',6X,'z',10X,'zz',8X,'dz',/)
8000  Format (I3,3F10.3)
8100  Format (1X,I5,F9.3,2F12.7)
8200  Format (1X,'  ts     ele      usurf    ang       ubot    ang
     & ssurf       sbot    ')
8300  Format (/2X,'  attenzione :job successfully completed; time = ',1
     &    P1E10.2,' days',//)
8400  Format (/2x,'  attenzione :job abnormal; time= ',E10.2,'IJ= ',I8,
     &         2X,'max velocity= ',E10.2,2x, 'position:',2F15.2)
     
      End
      
