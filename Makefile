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

ifneq ($(origin HOME), environment)
$(error Environment variable HOME was not set.)
endif

# check for presence of required FISOC OM env vars
ifneq ($(origin FISOC_OM), environment)
$(error Environment variable FISOC_OM was not set.)
endif

# check for presence of required FISOC ISM env vars
ifneq ($(origin FISOC_ISM), environment)
$(error Environment variable FISOC_ISM was not set.)
endif

ifeq ($(FISOC_ISM), Elmer)
 ifneq ($(origin ELMER_SO), environment)
  $(error Environment variable ELMER_SO was not set.) # required for linking to Elmer
 endif
 ifeq ($(origin ELMER_HOME), environment)
  ELMER_INCLUDE = -I$(ELMER_HOME)/share/elmersolver/include
 else
  $(error Environment variable ELMER_HOME was not set.)
 endif
else
 ELMER_INCLUDE = 
 ELMER_SO = 
endif

include $(ESMFMKFILE)

################################################################################

.SUFFIXES: .f90
%.o : %.f90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ELMER_INCLUDE) $(ESMF_F90COMPILEFREENOCPP) -cpp  -o $@  $<

################################################################################


FISOC_caller: $(SRCDIR)/FISOC_caller.o $(SRCDIR)/FISOC_parent.o $(SRCDIR)/FISOC_OM.o  $(SRCDIR)/FISOC_ISM.o $(SRCDIR)/FISOC_ElmerWrapper.o $(SRCDIR)/FISOC_coupler.o $(ELMER_SO)
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS)  
$(SRCDIR)/FISOC_caller.o:       $(SRCDIR)/FISOC_parent.o $(SRCDIR)/FISOC_caller.f90

$(SRCDIR)/FISOC_parent.o:       $(SRCDIR)/FISOC_parent.f90 $(SRCDIR)/FISOC_OM.o  $(SRCDIR)/FISOC_ISM.o $(SRCDIR)/FISOC_coupler.o
$(SRCDIR)/FISOC_coupler.o:      $(SRCDIR)/FISOC_coupler.f90
$(SRCDIR)/FISOC_ISM.o:          $(SRCDIR)/FISOC_ISM.f90 $(SRCDIR)/FISOC_$(FISOC_ISM)Wrapper.o
$(SRCDIR)/FISOC_OM.o:           $(SRCDIR)/FISOC_OM.f90
$(SRCDIR)/FISOC_$(FISOC_ISM)Wrapper.o:  $(SRCDIR)/FISOC_$(FISOC_ISM)Wrapper.f90

################################################################################

install: FISOC_caller
	cp FISOC_caller $(HOME)/bin

.PHONY: clean

clean:
	rm -f FISOC_caller *.o *.mod $(SRCDIR)/*.o $(SRCDIR)/*.mod


