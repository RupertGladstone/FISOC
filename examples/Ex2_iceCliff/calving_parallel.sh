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
rm -f surface.xyz bedrock.xyz ${sifname}
ln -s surface_step0.xyz surface.xyz
ln -s bedrock_step0.xyz bedrock.xyz
namerun="${basename}_step0"
echo "\$namerun = \"${namerun}\" " > ${sifname}
cat flat_init_blueprint.sif >> ${sifname}
if [ "$ncores" -gt "1" ] 
then    
    mpirun -np ${ncores} ElmerSolver_mpi > par${ncores}_step_0.log 
else
    ElmerSolver > ser_step_0.log
fi
echo "Step 0 done"
#
# Initialization finished
#
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
    # parallel run
    if [ "$ncores" -gt "1" ] 
    then
	griddatname="calv_1200-1100m_par${ncores}_step${j}_grid"
	# backing up existing file from previous run (you never know)
	if [ -e "${griddatname}.dat" ]; then
	    echo "moving existing ${griddatname}.dat to ${griddatname}.dat.bak"
	    mv ${griddatname}.dat ${griddatname}.dat.bak
	fi
	# stitching the previous surface data files together
	# (from ncores to single)
	touch ${griddatname}
	for ((k=1;k<=$ncores;k++)){
	    if [ $k -lt 10 ]; then
		cat ${griddatname}par000${k}.dat >> ${griddatname}.dat
	    elif [$k -gt 9] && [$k -lt 100]; then
		cat ${griddatname}par00${k}.dat >> ${griddatname}.dat
	    elif [$k -gt 99] && [$k -lt 1000]; then
		cat ${griddatname}par0${k}.dat >> ${griddatname}.dat
	    elif [ $k -gt 9999]; then
		echo "more than 9999 cores, really?"
		exit -1
	    else
		cat ${griddatname}par${k}.dat >> ${griddatname}.dat
	    fi
	}	
	# mesh partitioning
	if [ -e "$meshfilename" ]; then
	    echo "running using existing $meshfilename"
	else
	    echo "creating  $meshfilename"
	    ElmerGrid 2 2  footprint -metis $ncores 4
	fi
    # serial run
    else
	griddatname="calv_1200-1100m_ser_step${j}_grid"
    fi
    echo "sorting previous output file"
    awk '{x=$4;y=$5; if (x*x < 0.0001) {x=0.0};  if (y*y < 0.0001) {y=0.0}; print x, y, $6}' ${griddatname}.dat | sort -t ' ' -k1g,1 -k2g,2 > surface_step${i}.xyz

    awk '{x=$4;y=$5; if (x*x < 0.0001) {x=0.0};  if (y*y < 0.0001) {y=0.0}; print x, y, $8}' ${griddatname}.dat | sort -t ' ' -k1g,1 -k2g,2 > bedrock_step${i}.xyz
    rm -f surface.xyz bedrock.xyz ${sifname}
    ln -s surface_step${i}.xyz surface.xyz
    ln -s bedrock_step${i}.xyz bedrock.xyz
    namerun="${basename}_step${i}"
    namerestart="${basename}_step${j}"
    echo "\$namerun = \"${namerun}\" " > ${sifname}
    echo "\$namerestart = \"${namerestart}\" " >> ${sifname}
    cat flat_coupled_blueprint.sif >> ${sifname}
    if [ "$ncores" -gt "1" ] 
    then    
	echo "launching the parallel Elmer run with command:"
	echo "mpirun -np ${ncores} ElmerSolver_mpi > par${ncores}_step_${i}.log" 
        mpirun -np ${ncores} ElmerSolver_mpi > par${ncores}_step_${i}.log 
    else
	echo "launching the serial Elmer run with command:"
	echo "ElmerSolver > ser_step_${i}.log"
	ElmerSolver > ser_step_${i}.log
    fi
    echo "Step ${i} done"
done
echo "######################################"
echo "#   All done -  Bye! "
echo "######################################"
#
# End loop calving runs
#
