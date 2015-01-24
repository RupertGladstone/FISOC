
import shutil
#import os
#print os.environ['ELMER_HOME']

ElmerSolver = "/home/svali-user/Work/Source/trunk/fem/src/ElmerSolver.src"
FISOC_ElmerSolver = "FISOC_ElmerSolver.f90"

shutil.copyfile(FISOC_ElmerSolver,ElmerSolver)

print "should probably try to replace on ly the relevant code section, but for now replace the whole file..."
