#!/bin/bash
# this is to make sure that sort behaves in a sane way
export LC_NUMERIC="C"
##############################
#
# Preparation, mainly checking,
# whether we have serial or 
# parallel run and setting
# names accordingly
#
# You could pass the number
# of cores (= # of MPI threads)
# as an argument to the script
#
##############################
# how many cores (use 1 for serial)
ncores=2
rounds=2
# just to be sure, to have a compatible shared object
elmerf90 getweights.f90 -o getweights.so

if [ "$ncores" -gt "1" ] 
then
    # preparation for parallel run
    echo "attempting parallel Elmer-run with $ncores partitions"
    meshfilename="footprint/partitioning.$ncores"
    if [ -e "$meshfilename" ]; then
	echo "running using existing $meshfilename"
    else
	echo "creating  $meshfilename"
	ElmerGrid 2 2  footprint -metis $ncores 4
    fi
    basename="calv_1200-1100m_par${ncores}"
    sifname="flat_coupled_par.sif"
else
    # preparation for serial run
    echo "attempting serial Elmer-run"
    basename="calv_1200-1100m_ser"
    sifname="flat_coupled.sif"
fi
echo ${sifname} > ELMERSOLVER_STARTINFO
###########################################
#
# Initialization run 
# is separated, as it doesn't need restart
#
###########################################
echo "Starting initialization run"
source calving_parallel_step0.sh
echo "Done initialization run"

#############################################################
#
# This is the loop over the separate calving runs 
#
#############################################################
for ((i=1;i<=$rounds;i++))
do
    echo "######################################"
    echo "#   Starting step ${i} out of ${rounds}"
    echo "######################################"
    j=$(expr $i - 1)

    #---------------------------------------------
    # THIS WOULD be the workflow, if we ran a general calving model
    #
    #     # prepare Elmer output DEM's for particle model
    # 
    #     # run (dummy) particle model
    #
    #     # convert output from particle model to new Elmer footprint mesh
    #
    #     # convert output from particle model to DEM's
    #           # in fact, you do not really have to alter
    #           #calving_parallel_getDEM.sh, as it doesn't
    #           # matter if the DEM exceeds the area of the
    #           # new footprint mesh
    #----------------------------------------------

    # INSTEAD we are taking the shortcut to directly morph the DEM
    # from previous step to the new one
    # calving is achieved by using the old footrpint, hence reseting
    # the calving front always to the initial x=1000 m
    echo "Preparing DEM's for step ${i}"
    source calving_parallel_getDEM.sh

    # running new time evolotion with Elmer/Ice
    echo "running Elmer/Ice for step ${i}"
    source calving_parallel_elmer.sh

    echo "Step ${i} done"
done
echo "######################################"
echo "#   All done -  Bye! "
echo "######################################"
#
# End loop calving runs
#
