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

ifneq ($(origin FISOC_ISM), environment)
$(error Environment variable FISOC_ISM was not set.)
endif

ifneq ($(origin FISOC_OM), environment)
$(error Environment variable FISOC_OM was not set.)
endif

ifneq ($(origin ELMER_SO), environment)
$(error Environment variable ELMER_SO was not set.)
endif

ifneq ($(origin HOME), environment)
$(error Environment variable HOME was not set.)
endif

include $(ESMFMKFILE)

#cp  $(SRCDIR)/FISOC_ISM_dummy.f90  $(SRCDIR)/FISOC_ISM.f90

################################################################################

.SUFFIXES: .f90
%.o : %.f90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) -I$(ELMER_HOME)/share/elmersolver/include  $(ESMF_F90COMPILEFREENOCPP) -cpp  -o $@  $<

################################################################################


FISOC_caller: $(SRCDIR)/FISOC_caller.o $(SRCDIR)/FISOC_parent.o $(SRCDIR)/FISOC_OM.o  $(SRCDIR)/FISOC_ISM.o $(SRCDIR)/FISOC_ElmerWrapper.o $(SRCDIR)/FISOC_coupler.o $(ELMER_SO)
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS)  
	rm  $(SRCDIR)/FISOC_ISM.f90  $(SRCDIR)/FISOC_OM.f90 
$(SRCDIR)/FISOC_caller.o:       $(SRCDIR)/FISOC_parent.o $(SRCDIR)/FISOC_caller.f90

$(SRCDIR)/FISOC_parent.o:       $(SRCDIR)/FISOC_parent.f90 $(SRCDIR)/FISOC_OM.o  $(SRCDIR)/FISOC_ISM.o $(SRCDIR)/FISOC_coupler.o
$(SRCDIR)/FISOC_coupler.o:      $(SRCDIR)/FISOC_coupler.f90
$(SRCDIR)/FISOC_ISM.o:          $(SRCDIR)/FISOC_ISM.f90 $(SRCDIR)/FISOC_ElmerWrapper.o
$(SRCDIR)/FISOC_ISM.f90:        $(SRCDIR)/FISOC_ISM_$(FISOC_ISM).f90
	cp  $(SRCDIR)/FISOC_ISM_$(FISOC_ISM).f90  $(SRCDIR)/FISOC_ISM.f90
$(SRCDIR)/FISOC_OM.o:           $(SRCDIR)/FISOC_OM.f90
$(SRCDIR)/FISOC_OM.f90:         $(SRCDIR)/FISOC_OM_$(FISOC_OM).f90
	cp  $(SRCDIR)/FISOC_OM_$(FISOC_OM).f90  $(SRCDIR)/FISOC_OM.f90
$(SRCDIR)/FISOC_ElmerWrapper.o: $(SRCDIR)/FISOC_ElmerWrapper.f90

################################################################################

install: FISOC_caller
	cp FISOC_caller $(HOME)/bin

.PHONY: clean

clean:
	rm -f FISOC_caller *.o *.mod $(SRCDIR)/*.o $(SRCDIR)/*.mod
	rm -f $(SRCDIR)/FISOC_ISM.f90  $(SRCDIR)/FISOC_OM.f90


