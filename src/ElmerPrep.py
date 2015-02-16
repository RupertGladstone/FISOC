
import shutil
#import os
#print os.environ['ELMER_HOME']

# start using Elmer_HOME maybe...

ElmerSolver = "/home/elmeruser/Source/elmerfem-svn/fem/src/ElmerSolver.src"
#ElmerSolver = "/home/svali-user/Work/Source/trunk/fem/src/ElmerSolver.src"
Solver = "/home/elmeruser/Source/elmerfem-svn/fem/src/Solver.src"
#Solver = "/home/svali-user/Work/Source/trunk/fem/src/Solver.src"
FISOC_ElmerSolver = "FISOC_ElmerSolver.f90"
FISOC_Solver = "FISOC_Solver.f90"

shutil.copyfile(FISOC_ElmerSolver,ElmerSolver)
shutil.copyfile(FISOC_Solver,Solver)

#shutil.copyfile("testCaller.f90","/home/svali-user/Work/Source/trunk/fem/src/testCaller.f90")
#shutil.copyfile("FISOC_ElmerWrapper.f90","/home/svali-user/Work/Source/trunk/fem/src/FISOC_ElmerWrapper.f90")

print "should probably try to replace only the relevant code section, but for now replace the whole file..."
