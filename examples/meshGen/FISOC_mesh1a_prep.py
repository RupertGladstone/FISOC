#!/usr/bin/python
import os, sys, numpy

# create DEM
print "\ncreating DEM..."
# first some parameters used in the DEM and the grid creation
# Note nx,ny,nz indicate number of cells or elements.  Number of nodes is +1.

# xg - x atgrounding line 
# xc - x at calving front 
# x0 - x at start of domain
# b0 - bottom of ice at start of domain
# s0 - surface at start of domain
# bg,sg bottom and surface at grounding line
# bc,sc bottom and surfsace at calving front

numElements_fl = 10000

ny   = 6
#ylim =  20000. # in metres
ylim = 1.0 # for creating unit height mesh for use with structured mesh mapper

nx   = 100

squish = 1.

#xlim = 50000. # in metres
#bedShape = "linear" 
#DEM_dir    = "/home/elmeruser/Source/FISOC/examples/meshGen/FISOC_mesh1a_DEM"
#DEM_dir_fl    = "/home/elmeruser/Source/FISOC/examples/meshGen/FISOC_mesh1a_DEM_fl"
#gridName1D_fl = "FISOC_mesh1a_1D_flowline"
#gridName2D_fl = "FISOC_mesh1a_2D_flowline"
xlim = 1800000./squish # in metres
bedShape = "MM_overdeepened"
DEM_dir    = "/home/elmeruser/Source/FISOC/examples/meshGen/FISOC_mesh1b_DEM"
DEM_dir_fl    = "/home/elmeruser/Source/FISOC/examples/meshGen/FISOC_mesh1b_DEM_fl"
gridName1D_fl = "FISOC_mesh1b_1D_flowline"
gridName2D_fl = "FISOC_mesh1b_2D_flowline"

# FISOC meshes:
# 1a linear 
# 1b MISMIP flowline overdeepened

# defines slope and sea level for linear bed
nz    = 8
z0    = 100.
zlast = -900.

# initial thickness
H = 200.

#nx   = 40
#xlim = 1700000. # in metres
#bedShape = "MM_overdeepened" 

# used in determining floatation
rho_i = 910.0   # ice density in kg/m^3
rho_w = 1025.0  # ocean water density in kg/m^3

surfFile   = DEM_dir+"/surf.xyz"
bedFile    = DEM_dir+"/bed.xyz" # called bed because extrudemesh requires it, but actualy lower surface
bedrockFile= DEM_dir+"/bedrock.dat" # this really is the bedrock!
gridName2D    = "FISOC_mesh1a_2D"
gridName3D    = "FISOC_mesh1a_3D"
sif_bed       = "/home/elmeruser/Source/PYELA/sifLib/PYELA_sp_readBedrock.sif"

surfFile_fl   = DEM_dir_fl+"/upper_surf.dat"
bedFile_fl    = DEM_dir_fl+"/lower_surf.dat" 
bedrockFile_fl= DEM_dir_fl+"/bedrock.dat"
sif_bed_fl    = "/home/elmeruser/Source/PYELA/sifLib/PYELA_sp_readGeometry_fl.sif"

x0 = 0.0
xc = xlim
b0 = z0
s0 = z0 + H

# depth at grounding line:
# zg(xg) = z0 + (zlast-z0)*xg/xlim
# floatation: 
# rho_i*H = rho_w*(-zg) => zg = -rho_i*H/rho_w
#
zg = -rho_i*H/rho_w
bg = zg
sg = zg + H
bc = bg
sc = sg
print "depth at grounding line "+str(zg)
xg = (zg - z0)*xlim/((zlast-z0))
print "grounding line distance "+str(xg)

# valid bedshapes:
# linear - downsloping towards ocean
# MM_overdeepened - MISMIP polynomial overdeepened bed (designed for 1800km domain but seems fine with 1700km)
# b(x) = - [ 729 - 2184.8(x/750km)^2 + 1031.72(x/750km)^4 - 151.72(x/750km)^6 ]m

# MATC syntax...
# $function bedrock(x) {_bedrock = ( 729.0 - 2184.8 * (x/750000.0)^2 + 1031.72 * (x/750000.0)^4 - 151.72 * (x/750000.0)^6 ) }
# $function bedrock(x) {_bedrock = 100.0 + (-900.0-100.0) * x / 500000.0}


#if DEM_dir not in os.listdir('.'):
if not os.path.isdir(DEM_dir):
    ret = os.mkdir(DEM_dir)
# the +1 is because nx and ny refer to the number of elements, 
# but here we are writing values at nodes
xAxis = numpy.linspace(0,xlim*squish,nx+1)
yAxis = numpy.linspace(0,ylim*squish,ny+1)
bedrock_fl = numpy.zeros([nx+1])
surface_fl = numpy.zeros([nx+1])
bedrock  = numpy.zeros([nx+1, ny+1])
surface  = numpy.zeros([nx+1, ny+1])
if bedShape == "linear":
    bedrock_fl[:] = numpy.linspace(z0,zlast,nx+1)
    bedrock[:,0]  = numpy.linspace(z0,zlast,nx+1)
elif bedShape == "MM_overdeepened":
    bedrock_fl[:] =    (729. - 2184.8*(xAxis[:]/750000.)**2 + 1031.72*(xAxis[:]/750000.)**4 - 151.72*(xAxis[:]/750000.)**6 )
    bedrock[:,0]  =    (729. - 2184.8*(xAxis[:]/750000.)**2 + 1031.72*(xAxis[:]/750000.)**4 - 151.72*(xAxis[:]/750000.)**6 )
else:
    sys.exit("ERROR: bedShape not recognised")
bedrock[:,:] = numpy.transpose(numpy.tile(bedrock[:,0],[ny+1,1]))
lowerSurf    = numpy.maximum(-H*rho_i/rho_w,bedrock[:,:])
lowerSurf_fl = numpy.maximum(-H*rho_i/rho_w,bedrock_fl[:])
surface[:,:] = lowerSurf[:,:]+H
surface_fl[:]= lowerSurf_fl[:]+H

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

surf_out   = open(surfFile_fl,'w')
lsurf_out  = open(bedFile_fl,'w')
bed_out    = open(bedrockFile_fl,'w')
for x in range(nx+1):
    surf_out.write(str(xAxis[x])+" 0.0 "+str(surface_fl[x])+"\n")
    lsurf_out.write(str(xAxis[x])+" 0.0 "+str(lowerSurf_fl[x])+"\n")
    bed_out.write(str(xAxis[x])+" 0.0 "+str(bedrock_fl[x])+"\n")
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
sif.write('\nSolver 2')
sif.write('\n  Equation = \"GroundedMaskInit\"')
sif.write('\n  Exported Variable 1 = -dofs 1 bedrock')
sif.write('\nEnd\n')
sif.close()

# Create sif section for reading the geometry (bed, and lower and upper surfaces).
# This is for internal extrusion and mesh mapper.
print "\ncreating sif section for flowline geometry reading..."
sif = open(sif_bed_fl,'w')
sif.write('\nSolver 1')
sif.write('\n  Exec Solver = Before Simulation')
sif.write('\n  Equation = \"Read Geometry\"')
sif.write('\n  Procedure = \"ElmerIceSolvers\" \"Grid2DInterpolator\"')
sif.write('\n  Variable 1 = String "bedrock"')
sif.write('\n  Variable 1 data file = File \"'+bedrockFile_fl+'\"')
sif.write('\n  Variable 1 x0 = Real 0.0')
sif.write('\n  Variable 1 y0 = Real 0.0')
sif.write('\n  Variable 1 lx = Real '+str(xlim))
sif.write('\n  Variable 1 ly = Real 0.0')
sif.write('\n  Variable 1 Nx = Integer '+str(nx+1))
sif.write('\n  Variable 1 Ny = Integer 1')
sif.write('\n  Variable 2 = String "lsurfInit"')
sif.write('\n  Variable 2 data file = File \"'+bedFile_fl+'\"')
sif.write('\n  Variable 2 x0 = Real 0.0')
sif.write('\n  Variable 2 y0 = Real 0.0')
sif.write('\n  Variable 2 lx = Real '+str(xlim))
sif.write('\n  Variable 2 ly = Real 0.0')
sif.write('\n  Variable 2 Nx = Integer '+str(nx+1))
sif.write('\n  Variable 2 Ny = Integer 1')
sif.write('\n  Variable 3 = String "usurfInit"')
sif.write('\n  Variable 3 data file = File \"'+surfFile_fl+'\"')
sif.write('\n  Variable 3 x0 = Real 0.0')
sif.write('\n  Variable 3 y0 = Real 0.0')
sif.write('\n  Variable 3 lx = Real '+str(xlim))
sif.write('\n  Variable 3 ly = Real 0.0')
sif.write('\n  Variable 3 Nx = Integer '+str(nx+1))
sif.write('\n  Variable 3 Ny = Integer 1')
sif.write('\nEnd\n')
sif.write('\nSolver 2')
sif.write('\n  Equation = \"GroundedMaskInit\"')
sif.write('\n  Exported Variable 1 = -dofs 1 bedrock')
sif.write('\n  Exported Variable 2 = -dofs 1 lsurfInit')
sif.write('\n  Exported Variable 3 = -dofs 1 usurfInit')
sif.write('\nEnd\n')
sif.close()


## create sif section for defining the bedrock with MATC
#if bedShape == "linear":
#    print "\ncreating sif section for flowline bedrock function..."
#    print "\n (fiel is "+sif_bed_fl+")"
#    sif = open(sif_bed_fl,'w')
#    sif.write('$function bedrock(x) { ')
#    sif.write('_bedrock = '+str(z0)+' + ('+str(zlast)+'-'+str(z0)+') * x / '+str(xlim)+' ')
#    sif.write('}')
#    sif.write('\nBody 2')
#    sif.write('\n  Name = "lower_surface"')
#    sif.write('\n  Material = 2')
#    sif.write('\nEnd')
#    sif.write('\nMaterial 2')
#    sif.write('\n  Min Zs Bottom = Variable Coordinate 1')
#    sif.write('\n    Real MATC "bedrock(tx)"')
#    sif.write('\nEnd')
#    sif.close()
#else:
#    print "WARNING: can only write linear MATC function so far..."

## create sif section for reading the bedrock
#print "\ncreating sif section for flowline bedrock reading..."
#sif = open(sif_bed_fl,'w')
#sif.write('\nSolver 1')
#sif.write('\n  Exec Solver = Before Simulation')
#sif.write('\n  Equation = \"Read Bedrock\"')
#sif.write('\n  Procedure = \"ElmerIceSolvers\" \"Grid2DInterpolator\"')
#sif.write('\n  Variable 1 = String "bedrock"')
#sif.write('\n  Variable 1 data file = File \"'+bedrockFile+'\"')
#sif.write('\n  Variable 1 x0 = Real 0.0')
#sif.write('\n  Variable 1 y0 = Real 0.0')
#sif.write('\n  Variable 1 lx = Real '+str(xlim))
#sif.write('\n  Variable 1 ly = Real 0.0')
#sif.write('\n  Variable 1 Nx = Integer '+str(nx+1))
#sif.write('\n  Variable 1 Ny = Integer 1')
#sif.write('\nEnd\n')
#sif.write('\nSolver 2')
#sif.write('\n  Equation = \"GroundedMaskInit\"')
#sif.write('\n  Exported Variable 1 = -dofs 1 bedrock')
#sif.write('\nEnd\n')
#sif.close()


# create flowline mesh
print "\ncreating 1D mesh..."
f_in = open("FISOC_mesh1_1D_template.grd","r")
f_in_content = f_in.readlines()
f_in.close()

f_out = open(gridName1D_fl+".grd","w")
for line in f_in_content:
    line=line.replace("XXX_nx",str(nx))
    line=line.replace("XXX_xlim",str(xlim))
    f_out.write(line)
f_out.close()

cmd = "ElmerGrid 1 2 "+gridName1D_fl+".grd > ElmerGrid.log"
print "\nstarting ElmerGrid with following command :\n",cmd
ret = os.system(cmd)
if ret == 0.:
    print "ElmerGrid complete, see ElmerGrid.log"
else:
    sys.exit("ElmerGrid appeared to fail, see ElmerGrid.log")

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

# create 2D flowline mesh
print "\ncreating 2D flowline mesh..."
f_in = open("flowlineTemplate.grd","r")
f_in_content = f_in.readlines()
f_in.close()

f_out = open(gridName2D_fl+".grd","w")
for line in f_in_content:
    line=line.replace("XXX_H" ,str(H))
    line=line.replace("XXX_x0",str(x0))
    line=line.replace("XXX_xg",str(xg))
    line=line.replace("XXX_xc",str(xc))
    line=line.replace("XXX_b0",str(b0))
    line=line.replace("XXX_bg",str(bg))
    line=line.replace("XXX_bc",str(bc))
    line=line.replace("XXX_s0",str(s0))
    line=line.replace("XXX_sg",str(sg))
    line=line.replace("XXX_sc",str(sc))
    line=line.replace("XXX_numElements",str(numElements_fl))
    f_out.write(line)
f_out.close()

cmd = "ElmerGrid 1 2 "+gridName2D_fl+".grd > ElmerGrid.log"
print "\nstarting ElmerGrid with following command :\n",cmd
ret = os.system(cmd)
if ret == 0.:
    print "ElmerGrid complete, see ElmerGrid.log"
else:
    sys.exit("ElmerGrid appeared to fail, see ElmerGrid.log")


# convert to ElmerPost format to check 3D mesh
cmd = "ElmerGrid 2 3 "+gridName2D_fl+" >> ElmerGrid.log"
print "\nstarting ElmerGrid with following command :\n",cmd
ret = os.system(cmd)
if ret == 0.:
    print "ElmerGrid complete, see ElmerGrid.log\n"
else:
    sys.exit("ElmerGrid appeared to fail, see ElmerGrid.log\n")

## create 3D mesh (extrude 2D mesh using DEM)
#print "\ncreating 3D mesh..."
#cmd = "ExtrudeMesh "+gridName2D+" "+gridName3D+" "+str(nz)+" 10. 1 0 0 0 1 "+DEM_dir+" 1000. 1.5 -999999. > ExtrudeMesh.log"
#print "\nstarting ExtrudeMesh with following command :\n",cmd
#ret = os.system(cmd)
#if ret == 0.:
#    print "ExtrudeMesh complete, see ExtrudeMesh.log"
#else:
#    sys.exit("ExtrudeMesh appeared to fail, see ExtrudeMesh.log")

# create 2D flowline mesh (extrude 1D flowline mesh using DEM)
#print "\ncreating 2D flowline mesh..."
#cmd = "ExtrudeMesh "+gridName1D_fl+" "+gridName2D_fl+" "+str(nz)+" 10. 1 0 0 0 1 "+DEM_dir+" 1000. 1.5 -999999. > ExtrudeMesh.log"
#print "\nstarting ExtrudeMesh with following command :\n",cmd
#ret = os.system(cmd)
#if ret == 0.:
#    print "ExtrudeMesh complete, see ExtrudeMesh.log"
#else:
#    sys.exit("ExtrudeMesh appeared to fail, see ExtrudeMesh.log")

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
