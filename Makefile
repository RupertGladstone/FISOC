# GNU Makefile for FISOC, assuming ESMF and all components have been compiled

################################################################################
## This Makefile must be able to find the "esmf.mk" Makefile fragment in the  ##
## 'include' line below.  This Makefile uses the "ESMFMKFILE" environment     ##
## variable.                                                                  ##
################################################################################

FISOC_EXE ?= FISOC_caller
SRCDIR = src

# check for presence of required env vars
ifneq ($(origin ESMFMKFILE), environment)
 $(error Environment variable ESMFMKFILE was not set.)
endif

ifneq ($(origin FISOC_ISM), environment)
 $(error Environment variable FISOC_ISM was not set.)
endif

ifneq ($(origin FISOC_AM), environment)
 $(error Environment variable FISOC_AM was not set.)
endif

ifneq ($(origin FISOC_OM), environment)
 $(error Environment variable FISOC_OM was not set.)
endif

ifneq ($(origin FISOC_OM_GEOM), environment)
 CPPFLAGS += -D FISOC_OM_GRID
else
 CPPFLAGS += -D $(FISOC_OM_GEOM)
endif

ifneq ($(origin FISOC_ISM_GEOM), environment)
 CPPFLAGS += -D FISOC_ISM_MESH
else
 CPPFLAGS += -D $(FISOC_ISM_GEOM)
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
	$(ESMF_F90COMPILER) $(FFLAGS) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP) -I$(FISOC_ISM_INCLUDE)  -I$(FISOC_OM_INCLUDE) $(CPPFLAGS) -cpp  -o $@  $<

info:
	$(info )
	$(info ******************************************************************)
	$(info *** FISOC, Framework for Ice Sheet Ocean model Coupling        ***)
	$(info *** Relevant variables for the build are now listed            ***)
	$(info ******************************************************************)
	$(info )
	$(info )
	$(info ******************************************************************)
	$(info Generic FISOC build information)
	$(info ESMFMKFILE        [${ESMFMKFILE}])
	$(info ESMF_F90COMPILER  [${ESMF_F90COMPILER}])
	$(info CPPFLAGS          [${CPPFLAGS}])
	$(info FFLAGS            [${FFLAGS}])
	$(info Additional ESMF Fortran compilation flags)
	$(info ...               [${ESMF_F90COMPILEOPTS}])
	$(info Additional ESMF Fortran link flags)
	$(info ...               [${ESMF_F90LINKOPTS}])
	$(info FISOC will be installed in)
	$(info ...               [${INSTALL_DIR}])
	$(info FISOC executable will be named )
	$(info ...               [${FISOC_EXE}])
	$(info )
	$(info )
ifeq ($(FISOC_ISM), dummy)
	$(info ******************************************************************)
	$(info Ice Sheet Model (ISM) is dummy, skipping ISM component)
else
	$(info ******************************************************************)
	$(info Ice Sheet Model (ISM) build information)
	$(info FISOC_ISM         [${FISOC_ISM}])
	$(info FISOC_ISM_LIBS    [${FISOC_ISM_LIBS}])
	$(info FISOC_ISM_LIBPATH [${FISOC_ISM_LIBPATH}])
	$(info FISOC_ISM_INCLUDE [${FISOC_ISM_INCLUDE}])
 ifeq ($(origin FISOC_ISM_GEOM), environment)
	$(info FISOC_ISM_GEOM    [${FISOC_ISM_GEOM}])
 endif
endif
	$(info )
	$(info )
ifeq ($(FISOC_OM), dummy)
	$(info ******************************************************************)
	$(info Ocean Model (OM) is dummy, skipping OM component)
else
	$(info ******************************************************************)
	$(info Ocean Model (OM) build information)
	$(info FISOC_OM          [${FISOC_OM}])
	$(info FISOC_OM_LIBS     [${FISOC_OM_LIBS}])
	$(info FISOC_OM_LIBPATH  [${FISOC_OM_LIBPATH}])
	$(info FISOC_OM_INCLUDE  [${FISOC_OM_INCLUDE}])
 ifeq ($(origin FISOC_OM_GEOM), environment)
	$(info FISOC_OM_GEOM     [${FISOC_OM_GEOM}])
 endif
endif
	$(info )
	$(info )
ifeq ($(FISOC_AM), dummy)
	$(info ******************************************************************)
	$(info Atmosphere Model (AM) is dummy, skipping AM component)
else
	$(info ******************************************************************)
	$(info Atmosphere Model (AM) build information)
	$(info FISOC_AM          [${FISOC_AM}])
	$(info FISOC_AM_LIBS     [${FISOC_AM_LIBS}])
	$(info FISOC_AM_LIBPATH  [${FISOC_AM_LIBPATH}])
	$(info FISOC_AM_INCLUDE  [${FISOC_AM_INCLUDE}])
 ifeq ($(origin FISOC_AM_GEOM), environment)
	$(info FISOC_AM_GEOM     [${FISOC_AM_GEOM}])
 endif
endif
	$(info )
	$(info )

################################################################################

$(FISOC_EXE): $(SRCDIR)/FISOC_caller.o $(SRCDIR)/FISOC_parent.o $(SRCDIR)/FISOC_OM.o  $(SRCDIR)/FISOC_OM_Wrapper_$(FISOC_OM).o $(SRCDIR)/FISOC_ISM.o $(SRCDIR)/FISOC_ISM_Wrapper_$(FISOC_ISM).o $(SRCDIR)/FISOC_coupler.o $(SRCDIR)/FISOC_utils.o
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS) -L$(FISOC_ISM_LIBPATH) -L$(FISOC_OM_LIBPATH) $(FISOC_ISM_LIBS) $(FISOC_OM_LIBS) $(FFLAGS)

$(SRCDIR)/FISOC_caller.o:                   $(SRCDIR)/FISOC_parent.o $(SRCDIR)/FISOC_caller.f90
$(SRCDIR)/FISOC_parent.o:                   $(SRCDIR)/FISOC_parent.f90 $(SRCDIR)/FISOC_OM.o  $(SRCDIR)/FISOC_ISM.o $(SRCDIR)/FISOC_coupler.o
$(SRCDIR)/FISOC_coupler.o:                  $(SRCDIR)/FISOC_coupler.f90
$(SRCDIR)/FISOC_ISM.o:                      $(SRCDIR)/FISOC_ISM.f90 $(SRCDIR)/FISOC_ISM_Wrapper_$(FISOC_ISM).o $(SRCDIR)/FISOC_utils.o
$(SRCDIR)/FISOC_OM.o:                       $(SRCDIR)/FISOC_OM.f90 $(SRCDIR)/FISOC_OM_Wrapper_$(FISOC_OM).o  $(SRCDIR)/FISOC_utils.o
$(SRCDIR)/FISOC_ISM_Wrapper_$(FISOC_ISM).o: $(SRCDIR)/FISOC_ISM_Wrapper_$(FISOC_ISM).f90 $(SRCDIR)/FISOC_types.o $(SRCDIR)/FISOC_utils.o 
$(SRCDIR)/FISOC_OM_Wrapper_$(FISOC_OM).o:   $(SRCDIR)/FISOC_OM_Wrapper_$(FISOC_OM).f90 $(SRCDIR)/FISOC_types.o $(SRCDIR)/FISOC_utils.o
$(SRCDIR)/FISOC_types.o:                    $(SRCDIR)/FISOC_types.f90
$(SRCDIR)/FISOC_utils.o:                    $(SRCDIR)/FISOC_utils.f90

################################################################################

install: info $(FISOC_EXE) 
	$(info )
	$(info Building FISOC complete)
	cp $(FISOC_EXE) $(INSTALL_DIR)/$(FISOC_EXE) 

.PHONY: clean

clean:
	rm -f FISOC_caller* *.o *.mod $(SRCDIR)/*.o $(SRCDIR)/*.mod


