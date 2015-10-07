
FC = gfortran

.SUFFIXES: .f90
%.o : %.f90
	$(FC) -c -o $@ $<

FISh_caller:      FISh_caller.o FISh.o
	$(FC) -o $@ $^

FISh_caller.o:    FISh.o FISh_caller.f90
FISh.o:           FISh.f90

install: FISh_caller
	cp FISh_caller $(INSTALL_DIR)


.PHONY: clean

clean:
	rm -f FISh_caller *.o *.mod 
