#!/usr/bin/python

from shutil import copyfile
from subprocess import call
from sys import exit
import numpy as np
import os


print "\n----------------------------------------------------------------------------------------"
print "FISOC test script.  Uses example 5."
print "Runs on 2 procs. Takes a few minutes."
print "Check the Verification ascii file afterwards."
print "This is intended to be a fail-fast script rather than a robust or flexible script."
print "ROMS must have been compiled with the ICESHELF2D_TOY_GL application."
print "Full path to the ROMS .in file must be given in ROMSinFileName.asc "
print "----------------------------------------------------------------------------------------"

print "\nCopy the ROMS input file..."
try:
    f = open("./ROMSinFileName.asc", "r")
    ROMSinFileName = f.readlines()[0] 
    f.close()
    print "Found ROMS in file path: ", ROMSinFileName.strip(), "\n"
except Exception, e:
    exit(repr(e))
try:
    copyfile(ROMSinFileName.strip(), "./ocean_iceshelf2d_toy_gl.in")
except Exception, e:
    exit(repr(e))

print "Copy and compile the code for generating idealised ice sheet and bedrock geometry..."
try:
    copyfile("../../examples/Ex1_LongThinMIS/FISOC_Elmer_geometries.f90", "./FISOC_Elmer_geometries.f90")
except Exception, e:
    exit(repr(e))
try:    
    call(["elmerf90","FISOC_Elmer_geometries.f90", "-o FISOC_Elmer_geometries.so"])
except Exception, e:
    exit(repr(e))

print "\nCopy the required example 5 solver input files and mesh files for Elmer/Ice..."
try:
    copyfile("../../examples/Ex5_BensBox_gl/FISOC_Ex5a.sif", "./FISOC_Ex5a.sif")
    copyfile("../../examples/Ex5_BensBox_gl/FISOC_Ex5b.sif", "./FISOC_Ex5b.sif")
    copyfile("../../examples/Ex5_BensBox_gl/FISOC_Ex5.grd", "./FISOC_Ex5.grd")
except Exception, e:
    exit(repr(e))

print "\nBuild and partition the Elmer mesh (coarse resolution on 2 procs)..."
try:
    f_stdout = open("./ElmerMakeGrid_stdout.asc", "w")
    call(["ElmerGrid", "1", "2", "FISOC_Ex5.grd"],stdout = f_stdout)
    f_stdout.close()
except Exception, e:
    exit(repr(e))
try:
    f_stdout = open("./ElmerPartitionGrid_stdout.asc", "w")
    call(["ElmerGrid", "2", "2", "FISOC_Ex5", "-partition", "1", "2", "0", "2"],stdout = f_stdout)
    f_stdout.close()
except Exception, e:
    exit(repr(e))

print "\nCopy config file for FISOC Ex5b..."
try:
    copyfile("../../examples/Ex5_BensBox_gl/FISOC_config.rc", "./FISOC_config.rc")
except Exception, e:
    exit(repr(e))

print "\nRun example 5a, an Elmer only solve for temperature..."
try:
    f = open("./ELMERSOLVER_STARTINFO", "w")
    f.write("FISOC_Ex5a.sif\n")
    f.close()
except Exception, e:
    exit(repr(e))
try:
    f_stdout = open("./FISOC_Ex5a_stdout.asc", "w")
    call(["mpirun", "-np", "2", "ElmerSolver"],stdout = f_stdout)
#    call(["srun", "ElmerSolver"],stdout = f_stdout)
    f_stdout.close()
except Exception, e:
    exit(repr(e))

print "\nVerify metrics for Ex5a..."
try:
    target_data = np.genfromtxt("target_metrics_Ex5a.dat")
    target_temps = target_data[:,0]
    target_vol = target_data[0,1]
    data = np.genfromtxt("metrics_Ex5a.dat")
    temps = data[:,0]
    vol = data[0,1]
except Exception, e:
    exit(repr(e))
try:
    f = open("./Verification_Ex5.asc", "w")
    msg = "Verification for Ex5a\n"
    f.write(msg)
    msg = "Volume discrepancy: "+str((target_vol - vol)/target_vol)+"\n"
    f.write(msg)
    td = str(np.linalg.norm(target_temps-temps)/np.linalg.norm(target_temps))
    msg = "Temperature discrepancy: "+td+"\n"
    f.write(msg)
except Exception, e:
    exit(repr(e))

print "\nRun example 5b, a coupled solve..."
try:
    os.remove("PET0.FISOC.Log")
    os.remove("PET1.FISOC.Log")
except OSError, e:
    print repr(e)
    pass
except Exception, e:
    exit(repr(e))
try:
    f_stdout = open("./FISOC_Ex5b_stdout.asc", "w")
    call(["mpirun", "-np", "2", "FISOC_caller"],stdout = f_stdout)
#    call(["srun", "FISOC_caller"],stdout = f_stdout)
    f_stdout.close()
except Exception, e:
    exit(repr(e))

print "\nVerify metrics for Ex5b..."
try:
    target_data = np.genfromtxt("target_metrics_Ex5b.dat")
    target_temps = target_data[:,0]
    target_vols = target_data[:,1]
    data = np.genfromtxt("metrics_Ex5b.dat")
    temps = data[:,0]
    vols = data[:,1]
except Exception, e:
    exit(repr(e))
try:
    msg = "Verification for Ex5b\n"
    f.write(msg)
    vd = str(np.linalg.norm(target_vols-vols)/np.linalg.norm(target_vols))
    msg = "Volume discrepancy: "+str(vd)+"\n"
    f.write(msg)
    td = str(np.linalg.norm(target_temps-temps)/np.linalg.norm(target_temps))
    msg = "Temperature discrepancy: "+td+"\n"
    f.write(msg)
    f.close()
except Exception, e:
    exit(repr(e))

print "\n"
