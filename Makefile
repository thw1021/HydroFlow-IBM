# CIP LINUX:
  FCOMPILER1 = ifort #gfortran ifort #g77
  FCOMPILER2 = gfortran #gfortran ifort #g77 
#  FFALGS  = -g -shared-intel -mcmodel=large -openmp  #-p 
  FFLAGS1  = -qopenmp -O2 -fp-model source -ipo -extend-source 132 -shared-intel -mcmodel=large #-g -shared-intel -check bounds #-axSSE4.1   #-shared-intel -mcmodel=medium  #-p 
  FFLAGS2  = -O2 -fp-model source -ipo -extend-source 132  #-fopenmp #-O3  #-shared-intel -mcmodel=medium  #-p 
  FLINK  = ${FFLAGS1} #-lgfortran #-static #-lgfortran
#  DIR1 = /home/zjx/water_entry/new
  INSTALLDIR = #/home/zy/omp/ocerm/
#

  EXOBJS = ADVC.o ADVTK.o ADVTKD.o ADVU.o ADVV.o ADVVIS.o ADVW.o ARCHIVE.o ATRDE.o ACHIVEDEM.o \
           BCDATA.o BCOND.o BCDEM.o BRINV.o DEMM.o DESSA.o DESSST.o DYN.o ELTION.o FCOUPLING.o FIRST.o \
           GETCOR.o GRAD.o GSDL.o HARDMODULE.o IPPKPPDEM.o NGRN1.o OCERM.o PMOVE.o PROFC.o PROFTK.o \
           PROFTKD.o PROFV.o PROFVIS.o PROFW.o REUV.o SETDOM.o SGSMODEL.o SMOOTHING.o SMOOTHINGVER.o \
           SOLVE3DPOLCG.o SOLVEELFPOLCG.o SUBGRIDH.o SUBGRIDV.o TVDSCHEMEH.o TVDSCHEMEV.o UPDATEFLOW.o\
           UVFN.o VERTVEL.o WALLDRAG.o WAVEBREAKING.o WAVEGEN.o WREAL.o WSIGMA.o ZEROS.o GRADPC.o ARCHIVE_BIN.o\
           STATISTICS.o TURGEN.o INIVOR.o WENO.o BRINV1.o FCOUPLINGoriginal.o IPPKPPINITIAL.o NORMAL.o PORECAL.o\
           SEARCH.o SOFTMODULE0.o SOFTMODULE1.o SOFTMODULE2.o DEMSETTING.o IPJPKPDEM1.o IPJPKPDEM2.o BCDEM1.o \
           DEMSETREADING.o PARTICLE_INFO.o IBM.o IBMALDF.o IBMALGC.o IBMALIDC.o IBMALIDP.o IBMARCHIVE.o \
		   IBMINIT.o IBMUPDATE.o IBMVISUAL.o SAVEDIVERGE.o 
#
#  MYINCS =


ocerm:$(EXOBJS)
	$(FCOMPILER1) $(FLINK) $(EXOBJS) -o $(INSTALLDIR)HydroFlow

# SUFFIXES-DEFINITION:
%.o : %.F ./Include/OCERM_INF ./Include/VORGEN_INF
	$(FCOMPILER1) -c   $(FFLAGS1)  $<

update:
	find . -name "*.F" -exec touch {} \; 

clean:
	rm *.o  HydroFlow
         
#
