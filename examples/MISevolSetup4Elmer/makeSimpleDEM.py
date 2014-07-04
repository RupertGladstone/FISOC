#!/usr/bin/python
import os, sys, numpy

# this is where the digital elevation maps for the upper and lower ice surfaces and 
# the bedrock (which differs from lower surface under ice shelves) will be stored.
DEM_dir     = "/home/elmeruser/Source/FISOC/examples/MISevolSetup4Elmer/longThinMIS_DEM"
uSurfFile   = DEM_dir+"/usurf.dat"
lSurfFile   = DEM_dir+"/lsurf.dat"
bedrockFile = DEM_dir+"/bedrock.dat"
if not os.path.isdir(DEM_dir):
    ret = os.mkdir(DEM_dir)

# some hard coded parameters governing the DEMs
bedShape = "MM_overdeepened" # MM is short for MISMIP, so this is Schoof's overdeepened bed
#bedShape = "linear" 
H = 200. # initial thickness in metres
z0    = 100. # z0 and zlast are highest and lowest bed elevations for sue with linear bedShape.
zlast = -900.

# used in determining floatation
rho_i = 910.0   # ice density in kg/m^3
rho_w = 1025.0  # ocean water density in kg/m^3

# read some grid parameters from the elmergrid input file
f = open('longThinMIS.grd','r')
lines = f.readlines()
f.close()
noData=-999;nx=noData;ny=noData;xLen=noData;yLen=noData;dx=noData;dy=noData
for line in lines:
    if '=' in line:
        key, value = line.split('=')
        if key.strip() == 'Element Divisions 1':
            nx = int(value)
        if key.strip() == 'Element Divisions 2':
            ny = int(value)
        if key.strip() == 'Subcell Limits 1':
            xLen = int(float(value.split()[1]))
        if key.strip() == 'Subcell Limits 2':
            yLen = int(float(value.split()[1]))

if (nx!=noData & xLen!=noData):
    dx = xLen / nx
if (ny!=noData & yLen!=noData):
    dy = yLen / ny

# extend one block in the y direction because the grid interpolator (an elmerice solver) isnt really 
# set up well for grids that match exactly the domain size...
ny   = ny + 1
yLen = yLen + dy
print ny, yLen

#print nx,ny,xLen,yLen,dx,dy
#sys.exit(0)

# generate the DEMs
xAxis     = numpy.linspace(0,xLen,nx+1)
yAxis     = numpy.linspace(0,yLen,ny+1)
bedrock   = numpy.zeros([nx+1, ny+1])
surface   = numpy.zeros([nx+1, ny+1])
lowerSurf = numpy.zeros([nx+1, ny+1])
if bedShape == "linear":
    bedrock[:,0]  = numpy.linspace(z0,zlast,nx+1)
elif bedShape == "MM_overdeepened":
    bedrock[:,0]  =    (729. - 2184.8*(xAxis[:]/750000.)**2 + 1031.72*(xAxis[:]/750000.)**4 - 151.72*(xAxis[:]/750000.)**6 )
else:
    sys.exit("ERROR: bedShape not recognised")

bedrock[:,:] = numpy.transpose(numpy.tile(bedrock[:,0],[ny+1,1]))
lowerSurf    = numpy.maximum(-H*rho_i/rho_w,bedrock[:,:])
surface[:,:] = lowerSurf[:,:]+H


# write out the DEMs
surf_out   = open(uSurfFile,'w')
lsurf_out  = open(lSurfFile,'w')
bed_out    = open(bedrockFile,'w')
for y in range(ny+1):
    for x in range(nx+1):
        surf_out.write(str(xAxis[x])+" "+str(yAxis[y])+" "+str(surface[x,y])+"\n")
        lsurf_out.write(str(xAxis[x])+" "+str(yAxis[y])+" "+str(lowerSurf[x,y])+"\n")
        bed_out.write(str(xAxis[x])+" "+str(yAxis[y])+" "+str(bedrock[x,y])+"\n")
surf_out.close()
lsurf_out.close()
bed_out.close()


# after making the DEM, make some adjustments to the sif file 
# specific to DEM parameters, and make a copy of the sif file 
# that can be run purely for the purpose of outputting the full 
# 3D (extruded) mesh.

sifName = 'longThinMIS'

replacements = {
    'XXX_xLen'           : str(xLen),
    'XXX_yLen'           : str(yLen),
    'XXX_Nx'             : str(nx+1),
    'XXX_Ny'             : str(ny+1),
    'XXX_execNever'      : '',
    'XXX_execReadInputs' : 'Exec Solver = Before Simulation',
    'XXX_nt'             : '5000',
    'XXX_outFreq'        : '10',
    'XXX_restartFile'    :  '',
    'XXX_restartPosition':  '',
    'XXX_restartB4Init'  :  '',
    }

restartReplacements = {
    'XXX_restartFile'    : '  Restart File = $restartName',
    'XXX_restartPosition': '  Restart position = 0',
    'XXX_restartB4Init'  : '  Restart Before Initial Conditions = Logical False',
    'XXX_execReadInputs' : 'Exec Solver = Never',
    }

meshOnlyReplacements = {
    'XXX_execNever'      : 'exec solver = never',
    'XXX_nt'             : '1',
    'XXX_outFreq'        : '1',
}

restartReplacements = dict(replacements.items() + restartReplacements.items())
meshOnlyReplacements = dict(replacements.items() + meshOnlyReplacements.items())

# read sif template then make replacements to create specific
# sifs for use
f = open(sifName+'.sif.template','r')
lines = f.readlines()
f.close()
f    = open(sifName+'.sif','w')
f_mo = open(sifName+'_meshOnly.sif','w')
f_re = open(sifName+'_restart.sif','w')
for line in lines:
    lineMO = line
    lineRE = line
    for key in replacements.keys():
        line   = line.replace(key,replacements[key])
    for key in restartReplacements.keys():
        lineRE = lineRE.replace(key,restartReplacements[key])
    for key in meshOnlyReplacements.keys():
        lineMO = lineMO.replace(key,meshOnlyReplacements[key])
    f.write(line)
    f_re.write(lineMO)
    f_mo.write(lineMO)
f.close()
f_mo.close()
f_re.close()
