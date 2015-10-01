#!/usr/bin/env bash
                                                                                                            
export FISOC_MPI="yes"

export FISOC_ISM="dummy"
export FISOC_ISM_LIBS=""
export FISOC_ISM_LIBPATH="$HOME"
export FISOC_ISM_INCLUDE="$HOME"

export FISOC_OM="ROMS"
export FISOC_OM_LIBS="-loceanM"
export FISOC_OM_LIBPATH="/usr/local/lib"
export FISOC_OM_INCLUDE="/home/elmeruser/Source/ROMSIceShelf/Build"

export ESMFMKFILE="$ESMF_DIR/DEFAULTINSTALLDIR/lib/libO/Linux.gfortran.64.openmpi.default/esmf.mk"

make clean
make install
