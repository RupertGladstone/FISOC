#!/bin/bash
#################################
# run new Elmer/Ice simulation
#################################
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
