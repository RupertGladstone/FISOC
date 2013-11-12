#!/usr/bin/python
import os, sys, numpy

# create DEM
print "\ncreating DEM..."
# first some parameters used in the DEM and the grid creation
ny   = 6
nz   = 8
ylim =  20000. # in metres

nx   = 50
xlim = 500000. # in metres
bedShape = "linear" 

# defines slope and sea level for linear bed
z0    = 500.
zlast = -500.

# initial thickness
H = 200.

#nx   = 40
#xlim = 1700000. # in metres
#bedShape = "MM_overdeepened" 

# used in determining floatation
rho_i = 910.0   # ice density in kg/m^3
rho_w = 1025.0  # ocean water density in kg/m^3

DEM_dir    = "FISOC_mesh1a_DEM"
surfFile   = DEM_dir+"/surf.xyz"
bedFile    = DEM_dir+"/bed.xyz" # called bed because extrudemesh requires it, but actualy lower surface
bedrockFile= DEM_dir+"/bedrock.dat" # this really is the bedrock!
gridName2D = "FISOC_mesh1a_2D"
gridName3D = "FISOC_mesh1a_3D"
sif_bed    = "PYELA_readBedrock.sif"

# valid bedshapes:
# linear - downsloping towards ocean
# MM_overdeepened - MISMIP polynomial overdeepened bed (designed for 1800km domain but seems fine with 1700km)
# b(x) = - [ 729 - 2184.8(x/750km)^2 + 1031.72(x/750km)^4 - 151.72(x/750km)^6 ]m


if DEM_dir not in os.listdir('.'):
    ret = os.mkdir(DEM_dir)
# the +1 is because nx and ny refer to the number of elements, 
# but here we are writing values at nodes
xAxis = numpy.linspace(0,xlim,nx+1)
yAxis = numpy.linspace(0,ylim,ny+1)
bedrock  = numpy.zeros([nx+1, ny+1])
surface  = numpy.zeros([nx+1, ny+1])
if bedShape == "linear":
    bedrock[:,0] = numpy.linspace(z0,zlast,nx+1)
elif bedShape == "MM_overdeepened":
    bedrock[:,0] =  -  729. - 2184.8*(xAxis[:]/750000.)**2 + 1031.72*(xAxis[:]/750000.)**4 - 151.72*(xAxis[:]/750000.)**6 
else:
    sys.exit("ERROR: bedShape not recognised")
bedrock[:,:] = numpy.transpose(numpy.tile(bedrock[:,0],[ny+1,1]))
lowerSurf = numpy.maximum(-H*rho_i/rho_w,bedrock[:,:])
surface[:,:] = lowerSurf[:,:]+H
bedrock[:,:] = lowerSurf[:,:]
surf_out   = open(surfFile,'w')
lsurf_out  = open(bedFile,'w')
bed_out    = open(bedrockFile,'w')
for y in range(ny+1):
    for x in range(nx+1):
        surf_out.write(str(xAxis[x])+" "+str(yAxis[y])+" "+str(surface[x,y])+"\n")
        lsurf_out.write(str(xAxis[x])+" "+str(yAxis[y])+" "+str(lowerSurf[x,y])+"\n")
        bed_out.write(str(xAxis[x])+" "+str(yAxis[y])+" "+str(bedrock[x,y])+"\n")
surf_out.close()
lsurf_out.close()
bed_out.close()


# create sif section for reading the bedrock
print "\ncreating sif section for bedrock reading..."
sif = open(sif_bed,'w')
sif.write('\nSolver 1')
sif.write('\n  Exec Solver = Before Simulation')
sif.write('\n  Equation = \"Read Bedrock\"')
sif.write('\n  Procedure = \"ElmerIceSolvers\" \"Grid2DInterpolator\"')
sif.write('\n  Variable 1 = String "bedrock"')
sif.write('\n  Variable 1 data file = File \"'+bedrockFile+'\"')
sif.write('\n  Variable 1 x0 = Real 0.0')
sif.write('\n  Variable 1 y0 = Real 0.0')
sif.write('\n  Variable 1 lx = Real '+str(xlim))
sif.write('\n  Variable 1 ly = Real '+str(ylim))
sif.write('\n  Variable 1 Nx = Integer '+str(nx+1))
sif.write('\n  Variable 1 Ny = Integer '+str(ny+1))
sif.write('\nEnd\n')
sif.write('\n!Solver 2')
sif.write('\n!  Equation = \"Navier-Stokes\"')
sif.write('\n!  Exported Variable 1 = -dofs 1 bedrock')
sif.write('\n!End\n')
sif.close()


# create 2D mesh
print "\ncreating 2D mesh..."
f_in = open("FISOC_mesh1_2D_template.grd","r")
f_in_content = f_in.readlines()
f_in.close()

f_out = open(gridName2D+".grd","w")
for line in f_in_content:
    line=line.replace("XXX_nx",str(nx))
    line=line.replace("XXX_ny",str(ny))
    line=line.replace("XXX_xlim",str(xlim))
    line=line.replace("XXX_ylim",str(ylim))
    f_out.write(line)
f_out.close()

cmd = "ElmerGrid 1 2 "+gridName2D+".grd > ElmerGrid.log"
print "\nstarting ElmerGrid with following command :\n",cmd
ret = os.system(cmd)
if ret == 0.:
    print "ElmerGrid complete, see ElmerGrid.log"
else:
    sys.exit("ElmerGrid appeared to fail, see ElmerGrid.log")

# create 3D mesh (extrude 2D mesh using DEM)
print "\ncreating 3D mesh..."
cmd = "ExtrudeMesh "+gridName2D+" "+gridName3D+" "+str(nz)+" 10. 1 0 0 0 1 "+DEM_dir+" 1000. 1.5 -999999. > ExtrudeMesh.log"
print "\nstarting ExtrudeMesh with following command :\n",cmd
ret = os.system(cmd)
if ret == 0.:
    print "ExtrudeMesh complete, see ExtrudeMesh.log"
else:
    sys.exit("ExtrudeMesh appeared to fail, see ExtrudeMesh.log")

# convert to ElmerPost format to check 3D mesh
cmd = "ElmerGrid 2 3 "+gridName3D+" >> ElmerGrid.log"
print "\nstarting ElmerGrid with following command :\n",cmd
ret = os.system(cmd)
if ret == 0.:
    print "ElmerGrid complete, see ElmerGrid.log\n"
else:
    sys.exit("ElmerGrid appeared to fail, see ElmerGrid.log\n")


# note on resulting boundary labels:
# 1 - lower boundary (bedrock)
# 2 - upper surface
# 3 - y = 0 side wall
# 4 - x = 500km calving front
# 5 - y = 20km side wall
# 6 - x = 0 inland boundary

# elmerf90 DummySolver.f90 -o DummySolver.so
# http://elmerice.elmerfem.org/wiki/doku.php?id=userfunctions:buoyancy
