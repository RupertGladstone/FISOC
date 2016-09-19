#!/bin/bash
#################################
# prepare output of previous run
# to act as digital elevation
# model (DEM) for current run
#################################
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
