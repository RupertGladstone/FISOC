#!/bin/bash
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
