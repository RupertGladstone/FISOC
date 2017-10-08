# GNU Makefile for FISOC, assuming ESMF and all components have been compiled

################################################################################
## This Makefile must be able to find the "esmf.mk" Makefile fragment in the  ##
## 'include' line below.  This Makefile uses the "ESMFMKFILE" environment     ##
## variable.                                                                  ##
################################################################################

SRCDIR = src
FFLAGS += -fbacktrace -g -O0 -fbounds-check #-Wall
#FFLAGS += -O0 -g -fbacktrace -fcheck=all # -Wall
#FFLAGS += -fbacktrace -g -debug -DD  -O0 # -inline-debug-info
#FFLAGS += -g -check all -fpe0 -warn -traceback -debug extended
#FFLAGS += -O3 -xHost #-ipo


# check for presence of required env vars
ifneq ($(origin FISOC_EXE), environment)
 FISOC_EXE = FISOC_caller
endif

ifneq ($(origin ESMFMKFILE), environment)
 $(error Environment variable ESMFMKFILE was not set.)
endif

ifneq ($(origin FISOC_ISM), environment)
 $(error Environment variable FISOC_ISM was not set.)
endif

ifneq ($(origin FISOC_OM), environment)
 $(error Environment variable FISOC_OM was not set.)
endif

ifneq ($(origin FISOC_AM), environment)
 $(error Environment variable FISOC_AM was not set.)
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

$(info )
$(info ********************************************************************)
$(info *** Building FISOC, Framework for Ice Sheet Ocean model Coupling ***)
$(info *** Relevant parameters for the build are now listed             ***)
$(info ********************************************************************)
$(info )
$(info FISOC will be installed in [${INSTALL_DIR}]/[${FISOC_EXE}])
$(info )
$(info ESMFMKFILE        [${ESMFMKFILE}])
$(info CPPFLAGS          [${CPPFLAGS}])
$(info )
$(info FISOC_ISM         [${FISOC_ISM}])
$(info FISOC_ISM_LIBS    [${FISOC_ISM_LIBS}])
$(info FISOC_ISM_LIBPATH [${FISOC_ISM_LIBPATH}])
$(info FISOC_ISM_INCLUDE [${FISOC_ISM_INCLUDE}])
$(info )
$(info FISOC_OM          [${FISOC_OM}])
$(info FISOC_OM_LIBS     [${FISOC_OM_LIBS}])
$(info FISOC_OM_LIBPATH  [${FISOC_OM_LIBPATH}])
$(info FISOC_OM_INCLUDE  [${FISOC_OM_INCLUDE}])
$(info )
$(info FISOC_AM          [${FISOC_AM}])
$(info FISOC_AM_LIBS     [${FISOC_AM_LIBS}])
$(info FISOC_AM_LIBPATH  [${FISOC_AM_LIBPATH}])
$(info FISOC_AM_INCLUDE  [${FISOC_AM_INCLUDE}])
$(info )

################################################################################

.SUFFIXES: .f90
%.o : %.f90
	$(ESMF_F90COMPILER) $(FFLAGS) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP) -I$(FISOC_ISM_INCLUDE)  -I$(FISOC_OM_INCLUDE) -I$(FISOC_ISM_INCLUDE) $(CPPFLAGS) -cpp  -o $@  $<

################################################################################

$(FISOC_EXE): $(SRCDIR)/FISOC_caller.o $(SRCDIR)/FISOC_parent.o $(SRCDIR)/FISOC_OM.o  $(SRCDIR)/FISOC_OM_Wrapper_$(FISOC_OM).o $(SRCDIR)/FISOC_ISM.o $(SRCDIR)/FISOC_ISM_Wrapper_$(FISOC_ISM).o $(SRCDIR)/FISOC_AM.o $(SRCDIR)/FISOC_AM_Wrapper_$(FISOC_AM).o $(SRCDIR)/FISOC_coupler.o $(SRCDIR)/FISOC_utils.o
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS) -L$(FISOC_ISM_LIBPATH) -L$(FISOC_OM_LIBPATH) -L$(FISOC_AM_LIBPATH) $(FISOC_ISM_LIBS) $(FISOC_OM_LIBS) $(FISOC_AM_LIBS) $(FFLAGS)

$(SRCDIR)/FISOC_caller.o:                   $(SRCDIR)/FISOC_parent.o $(SRCDIR)/FISOC_caller.f90
$(SRCDIR)/FISOC_parent.o:                   $(SRCDIR)/FISOC_parent.f90 $(SRCDIR)/FISOC_OM.o  $(SRCDIR)/FISOC_ISM.o $(SRCDIR)/FISOC_AM.o $(SRCDIR)/FISOC_coupler.o
$(SRCDIR)/FISOC_coupler.o:                  $(SRCDIR)/FISOC_coupler.f90
$(SRCDIR)/FISOC_ISM.o:                      $(SRCDIR)/FISOC_ISM.f90 $(SRCDIR)/FISOC_ISM_Wrapper_$(FISOC_ISM).o $(SRCDIR)/FISOC_utils.o
$(SRCDIR)/FISOC_OM.o:                       $(SRCDIR)/FISOC_OM.f90 $(SRCDIR)/FISOC_OM_Wrapper_$(FISOC_OM).o  $(SRCDIR)/FISOC_utils.o
$(SRCDIR)/FISOC_AM.o:                       $(SRCDIR)/FISOC_AM.f90 $(SRCDIR)/FISOC_AM_Wrapper_$(FISOC_AM).o  $(SRCDIR)/FISOC_utils.o
$(SRCDIR)/FISOC_ISM_Wrapper_$(FISOC_ISM).o: $(SRCDIR)/FISOC_ISM_Wrapper_$(FISOC_ISM).f90 $(SRCDIR)/FISOC_types.o $(SRCDIR)/FISOC_utils.o 
$(SRCDIR)/FISOC_OM_Wrapper_$(FISOC_OM).o:   $(SRCDIR)/FISOC_OM_Wrapper_$(FISOC_OM).f90 $(SRCDIR)/FISOC_types.o $(SRCDIR)/FISOC_utils.o
$(SRCDIR)/FISOC_AM_Wrapper_$(FISOC_AM).o:   $(SRCDIR)/FISOC_AM_Wrapper_$(FISOC_AM).f90 $(SRCDIR)/FISOC_types.o $(SRCDIR)/FISOC_utils.o
$(SRCDIR)/FISOC_types.o:                    $(SRCDIR)/FISOC_types.f90
$(SRCDIR)/FISOC_utils.o:                    $(SRCDIR)/FISOC_utils.f90

################################################################################

install: $(FISOC_EXE)
	cp $(FISOC_EXE) $(INSTALL_DIR)/$(FISOC_EXE)

.PHONY: clean

clean:
	rm -f FISOC_caller* *.o *.mod $(SRCDIR)/*.o $(SRCDIR)/*.mod


