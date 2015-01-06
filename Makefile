# GNU Makefile for FISOC, assuming ESMF and all components have been compiled

################################################################################
## This Makefile must be able to find the "esmf.mk" Makefile fragment in the  ##
## 'include' line below.  This Makefile uses the "ESMFMKFILE" environment     ##
## variable.                                                                  ##
################################################################################

SRCDIR = src
FFLAGS = 

ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif

include $(ESMFMKFILE)


################################################################################

.SUFFIXES: .f90
%.o : %.f90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP) -cpp  -o $@  $<

################################################################################

# dependencies 
#FISOC_caller: FISOC_caller.o FISOC_parent.o FISOC_coupler.o FISOC_proc.o FISOC_Elmer_dummy.o FISOC_ROMS_dummy.o  
#FISOC_caller: FISOC_caller.o FISOC_parent.o FISOC_coupler.o FISOC_proc.o FISOC_Elmer.o FISOC_ROMS.o  
#	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS)
#FISOC_caller.o: FISOC_parent.o  
#FISOC_parent.o: FISOC_coupler.o FISOC_proc.o FISOC_Elmer_dummy.o FISOC_ROMS_dummy.o  
#FISOC_parent.o: FISOC_coupler.o FISOC_proc.o FISOC_ROMS.o FISOC_Elmer.o  

FISOC_caller: $(SRCDIR)/FISOC_caller.o $(SRCDIR)/FISOC_parent.o $(SRCDIR)/FISOC_ISM.o  
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS)
#FISOC_caller.o: FISOC_parent.o  
$(SRCDIR)/FISOC_caller.o: $(SRCDIR)/FISOC_parent.o $(SRCDIR)/FISOC_caller.f90
$(SRCDIR)/FISOC_parent.o: $(SRCDIR)/FISOC_parent.f90 $(SRCDIR)/FISOC_ISM.o
$(SRCDIR)/FISOC_ISM.o:  $(SRCDIR)/FISOC_ISM.f90

################################################################################

.PHONY: clean

clean:
	rm -f FISOC_caller *.o *.mod $(SRCDIR)/*.o $(SRCDIR)/*.mod


