#!/usr/bin/env bash

module purge

module unload openmpi intel-fc intel-cc
module load intel-mpi/5.1.0.079
#module load intel-fc/15.0.3.187
#module load intel-cc/15.0.3.187
#module load intel-mkl/15.0.3.187
module load intel-fc/16.0.2.181
module load intel-cc/16.0.2.181
module load intel-mkl/16.0.2.181
module load cmake
module load netcdf/4.3.3.1
#module load netcdf

export FISOC_EXE="FISOC_caller_IceOcean1ra"

#export LD_LIBRARY_PATH="/short/ks3/rmg581/elmer/install/lib/elmersolver/:$LD_LIBRARY_PATH"

#export ESMF_DIR="/short/gh8/lmj581/Software/esmf"
export ESMF_DIR="/short/ks3/rmg581/esmf_git"
export CPPFLAGS="$CPPFLAGS -D FISOC_MPI"
#export FISOC_MPI="yes"


export ELMER_HOME="/short/gh8/cxz581/Elmer_FISOC/install/"
#export ELMER_HOME="/short/ks3/rmg581/elmer/install"

export FISOC_ISM="Elmer"
export FISOC_ISM_LIBS="-lelmersolver"
export FISOC_ISM_INCLUDE="$ELMER_HOME/share/elmersolver/include"
export FISOC_ISM_LIBPATH="$ELMER_HOME/lib/elmersolver/"

export LD_LIBRARY_PATH="/short/gh8/cxz581/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$ELMER_HOME/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$ELMER_HOME/lib/elmersolver:$LD_LIBRARY_PATH"



export FISOC_AM="dummy"
export FISOC_AM_LIBS=""
export FISOC_AM_LIBPATH="$HOME"
export FISOC_AM_INCLUDE="$HOME"

#export FISOC_OM="dummy"

#export FISOC_OM_INCLUDE="$HOME"
#export FISOC_OM_LIBPATH="$HOME"

export MY_ROMS_DIR="/short/gh8/cxz581/ROMSIceShelf_devel_IceOcean1ra"
export FISOC_OM="ROMS"
export FISOC_OM_LIBS="-loceanM_isoPlus_o3_2"
export FISOC_OM_INCLUDE="${MY_ROMS_DIR}/Build"
export FISOC_OM_LIBPATH="/short/gh8/cxz581/lib"
# These ROMS_ preprocessor keywords correspond to a relevant subset of  
# the preprocessor keywords in the ROMS .in file.
#export CPPFLAGS="$CPPFLAGS -D ROMS_SPHERICAL"
#export CPPFLAGS="$CPPFLAGS -D ROMS_MASKING"
export CPPFLAGS="$CPPFLAGS -D ROMS_DSDT"
export CPPFLAGS="$CPPFLAGS -D ROMS_DDDT"
#export CPPFLAGS="$CPPFLAGS -D ROMS_DRAFT"
export CPPFLAGS="$CPPFLAGS -D ROMS_AVERAGES"

export ESMFMKFILE="$ESMF_DIR/DEFAULTINSTALLDIR/lib/libO/Linux.intel.64.intelmpi.default/esmf.mk"


make clean
make install

#rm PET*Log
export PATH="/home/581/cxz581/bin:$PATH"
# to run, for example:
#mpirun -np 16 FISOC_caller_MISOMIP1

