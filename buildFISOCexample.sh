#!/usr/bin/env bash
                                                                                                            
export FISOC_MPI="yes"

#export FISOC_ISM="dummy"
#export FISOC_ISM_LIBS=""
#export FISOC_ISM_LIBPATH="$HOME"
#export FISOC_ISM_INCLUDE="$HOME"

export FISOC_ISM="FISh"
export FISOC_ISM_LIBS="-lFISh"
export FISOC_ISM_LIBPATH="/usr/local/lib"
export FISOC_ISM_INCLUDE="$HOME"

#export FISOC_ISM="Elmer"
#export FISOC_ISM_LIBS="-lelmersolver"
#export FISOC_ISM_INCLUDE="$ELMER_HOME/share/elmersolver/include"
#export FISOC_ISM_LIBPATH="$ELMER_HOME/lib/elmersolver/"

export FISOC_OM="dummy"
export FISOC_OM_LIBS=""
export FISOC_OM_INCLUDE="$HOME"
export FISOC_OM_LIBPATH="$HOME"

#export FISOC_OM="ROMS"
#export FISOC_OM_LIBS="-loceanM"
#export FISOC_OM_INCLUDE="${MY_ROMS_DIR}/Build"
#export FISOC_OM_LIBPATH="$HOME/lib/"

export ESMFMKFILE="$ESMF_DIR/DEFAULTINSTALLDIR/lib/libO/Linux.gfortran.64.openmpi.default/esmf.mk"

make clean
make install

rm PET*Log

# to run, for example:
# mpirun -np 4 FISOC_caller

