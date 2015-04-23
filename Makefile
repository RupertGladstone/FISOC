# GNU Makefile for FISOC, assuming ESMF and all components have been compiled

################################################################################
## This Makefile must be able to find the "esmf.mk" Makefile fragment in the  ##
## 'include' line below.  This Makefile uses the "ESMFMKFILE" environment     ##
## variable.                                                                  ##
################################################################################

SRCDIR = src
FFLAGS = 

# check for presence of required env vars
ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif

ifneq ($(origin FISOC_OM), environment)
$(error Environment variable FISOC_OM was not set.)
endif

ifneq ($(origin FISOC_ISM), environment)
$(error Environment variable FISOC_ISM was not set.)
endif

ifneq ($(origin FISOC_INSTALL_DIR), environment)
 ifneq ($(origin HOME), environment)
  $(error Environment variables neither HOME nor FISOC_INSTALL_DIR was set.)
 else
  INSTALL_DIR = $(HOME)/bin
 endif
else
 INSTALL_DIR = $(FISOC_INSTALL_DIR)
endif

include $(ESMFMKFILE)

################################################################################

.SUFFIXES: .f90
%.o : %.f90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP) -I$(FISOC_ISM_INCLUDE)  -cpp  -o $@  $<

################################################################################

FISOC_caller: $(SRCDIR)/FISOC_caller.o $(SRCDIR)/FISOC_parent.o $(SRCDIR)/FISOC_OM.o  $(SRCDIR)/FISOC_OM_Wrapper_$(FISOC_OM).o $(SRCDIR)/FISOC_ISM.o $(SRCDIR)/FISOC_ISM_Wrapper_$(FISOC_ISM).o $(SRCDIR)/FISOC_coupler.o $(SRCDIR)/FISOC_utils.o
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS) -L$(FISOC_ISM_LIBPATH) $(FISOC_ISM_LIBS)

$(SRCDIR)/FISOC_caller.o:                   $(SRCDIR)/FISOC_parent.o $(SRCDIR)/FISOC_caller.f90
$(SRCDIR)/FISOC_parent.o:                   $(SRCDIR)/FISOC_parent.f90 $(SRCDIR)/FISOC_OM.o  $(SRCDIR)/FISOC_ISM.o $(SRCDIR)/FISOC_coupler.o
$(SRCDIR)/FISOC_coupler.o:                  $(SRCDIR)/FISOC_coupler.f90
$(SRCDIR)/FISOC_ISM.o:                      $(SRCDIR)/FISOC_ISM.f90 $(SRCDIR)/FISOC_ISM_Wrapper_$(FISOC_ISM).o $(SRCDIR)/FISOC_utils.o
$(SRCDIR)/FISOC_OM.o:                       $(SRCDIR)/FISOC_OM.f90 $(SRCDIR)/FISOC_OM_Wrapper_$(FISOC_OM).o  $(SRCDIR)/FISOC_utils.o
$(SRCDIR)/FISOC_ISM_Wrapper_$(FISOC_ISM).o: $(SRCDIR)/FISOC_ISM_Wrapper_$(FISOC_ISM).f90 $(SRCDIR)/FISOC_$(FISOC_ISM)_types.o  $(SRCDIR)/FISOC_utils.o 
$(SRCDIR)/FISOC_OM_Wrapper_$(FISOC_OM).o:   $(SRCDIR)/FISOC_OM_Wrapper_$(FISOC_OM).f90  $(SRCDIR)/FISOC_utils.o
$(SRCDIR)/FISOC_$(FISOC_ISM)_types.o:       $(SRCDIR)/FISOC_$(FISOC_ISM)_types.f90
$(SRCDIR)/FISOC_utils.o:                    $(SRCDIR)/FISOC_utils.f90


################################################################################

install: FISOC_caller
	cp FISOC_caller $(INSTALL_DIR)

.PHONY: clean

clean:
	rm -f FISOC_caller *.o *.mod $(SRCDIR)/*.o $(SRCDIR)/*.mod


