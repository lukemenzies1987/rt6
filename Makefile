SHELL=/bin/sh
NAME=	RT6
CMD=$(NAME)

SRCS= kindmod.f90 loopngrpintmod.f90 rt6.f90 nuclearvars.f90 onedimensionalvars.f90 twodimensionalvars.f90 clusterintmod.f90 DEMOspeciesgen.f90 functions.f90 functions1d.f90 writesdfforvoids.f90 getparameters.f90 writearraytofile.f90 writearraytofile2.f90 monitr.f90 groupingsetup.f90 groupingsetup1d.f90 initialsetup.f90 initialsetuporiginal.f90 initialsetup1d.f90 clusterdensity.f90 clusterdensityngrp.f90 initialvalues.f90 clusters.f90 printinputs.f90 speciesgeneration.f90 funcallhybridngrpss.f90 removeoldfile.f90 loop1.f90 loop1ngrp.f90 loop2.f90 loop2ngrp.f90 loop3.f90 sinkstrengthcorrection.f90 sinkstrengthcorrectionngrp.f90 speciesgenerationHe.f90 deallocatearrays.f90 deallocatearrays1d.f90 funcallhybrid.f90 funcallhybridngrp.f90 funcallhybrid1d.f90 tolerancechange.f90 tolerancechangengrp.f90 MDDFTterms.f90 jac.f90 s08_dlsode.f opkdmainmodified.f opkda1modified.f opkda2modified.f compressibility.f90 modifiedbessels.f90 clusterss.f90 #euler.f90

OBJS= kindmod.o loopngrpintmod.o clusterintmod.o nuclearvars.o onedimensionalvars.o twodimensionalvars.o functions.o DEMOspeciesgen.o functions1d.o removeoldfile.o writearraytofile.o writearraytofile2.o clusterdensity.o clusterdensityngrp.o getparameters.o monitr.o initialvalues.o loop1.o loop1ngrp.o loop2.o loop2ngrp.o loop3.o clusters.o rt6.o jac.o groupingsetup.o groupingsetup1d.o MDDFTterms.o initialsetup.o initialsetup1d.o printinputs.o funcallhybrid.o funcallhybridngrp.o funcallhybridngrpss.o funcallhybrid1d.o sinkstrengthcorrection.o clusterss.o sinkstrengthcorrectionngrp.o compressibility.o tolerancechange.o tolerancechangengrp.o speciesgeneration.o speciesgenerationHe.o deallocatearrays.o deallocatearrays1d.o #euler.o opkdmainmodified.o opkda1modified.o opkda2modified.o modifiedbessels.o

#------------------------------------------------------------------
#Nag fortran directory setup
#NAGFDIR   = /opt/NAG/fll6i25dcl
NAGFDIR   = /opt/NAG/fll6i26dcl
NAGLDIR   = ${NAGFDIR}/lib
NAGMODDIR = ${NAGFDIR}/nag_interface_blocks
FLINK     = ${NAGLDIR}/libnag_nag.a
#Additional libraries like lapack and blas if required.
LIBRARIES = $(FLINK) -L/home/mbgnklm3/Documents/ODEPACK/libodepack -lodepack -L/home/luke/Documents/lapack-3.6.0 -llapack -L/usr/lib -lblas 
DEBUG     =  #-g (debugger)
OPTIMIZE  =  -O3 # Recommended to run on 02 maximum 
LISTING   = # -listing
SPECIAL   = -fp-model precise -fp-model source #-ffpe-trap=invalid -p -fbounds-check (gfortran flags) #-p -traceback -fpe:0 -check all # 

# Change compiler to from ifort to gfortran if desired
COMPILE   = ifort $(SPECIAL) $(DEBUG) $(OPTIMIZE) -I$(NAGMODDIR) -c 
LOAD      = ifort $(SPECIAL) $(OPTIMIZE) $(LIBRARIES) -o $(CMD) $(LIBS)
ARCHIVE   = ar rv
CATALOG   = ar tv
#
.SUFFIXES:

.SUFFIXES: .o .f .f90

all: $(CMD)
	@echo
	@echo ">>  `date '+%a %d-%h-%y %r'`  `pwd`  `uname -mns`  $(LOGNAME)"

$(CMD):         $(OBJS) 
	$(LOAD) $(OBJS) -pg $(LIBRARIES) 
	chmod 755 $(CMD)
	rm *.o

.f.o:	
	$(COMPILE) $<

.f90.o:	
	$(COMPILE) $<
#-------------------------------------------------------------------------
clean:
	rm -f $(CMD) *.a *.o

tarfile:
	tar cfv $(NAME).tar makefile $(INPUT) $(INCL) $(SRCS) 

