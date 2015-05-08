
import shutil
#import os
#print os.environ['ELMER_HOME']

#ElmerSolver = "/home/elmeruser/Source/elmerfem-svn/fem/src/ElmerSolver.src"
#ElmerSolver = "/home/svali-user/Work/Source/trunk/fem/src/ElmerSolver.src"
ElmerSolver = "/home/elmeruser/Source/elmerfem-git/elmerice/fem/src/ElmerSolver.F90"

#Solver = "/home/elmeruser/Source/elmerfem-svn/fem/src/Solver.src"
#Solver = "/home/svali-user/Work/Source/trunk/fem/src/Solver.src"
Solver = "/home/elmeruser/Source/elmerfem-git/elmerice/fem/src/Solver.F90"

FISOC_ElmerSolver = "FISOC_ElmerSolver.f90"
#FISOC_Solver = "FISOC_Solver.f90"

print "Copying "+FISOC_ElmerSolver+" to "+ElmerSolver
shutil.move(ElmerSolver,ElmerSolver+".old")
shutil.copyfile(FISOC_ElmerSolver,ElmerSolver)

searchline = 'PROGRAM Solver\n'
newline    = '   USE ElmerSolver_mod\n'

ff = open(Solver,'r')
lines = ff.readlines()
ff.close()

if searchline not in lines:
    print "WARNING: can\'t find the searchline in "+Solver
if newline not in lines:
    print "Adding new line to "+Solver
    ii = lines.index(searchline)
    print "line",ii
    lines.insert(ii+1,newline)
    shutil.move(Solver,Solver+".old")
    ff = open(Solver,'w')
    lines = ff.writelines(lines)
    ff.close()
else:
    print "No changes needed to ",Solver

print "You may need to re-compile Elmer now"
