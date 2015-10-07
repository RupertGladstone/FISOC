
FC          = gfortran
FFLAGS     += -fPIC
INSTALL_DIR = $(HOME)/bin
LIBDIR      = /usr/local/lib/

.SUFFIXES: .f90
%.o : %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

FISh_caller:      FISh_caller.o FISh.o libFISh.so
	$(FC) $(FFLAGS) -o $@ $^

libFISh.so:       FISh.o
	$(FC) $(FFLAGS) -shared -o $@ $^

FISh_caller.o:    FISh.o FISh_caller.f90
FISh.o:           FISh.f90

install: FISh_caller libFISh.so
	cp FISh_caller $(INSTALL_DIR)
	cp libFISh.so $(LIBDIR)

.PHONY: clean

clean:
	rm -f FISh_caller *.o *.mod *.so
