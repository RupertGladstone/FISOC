#!/usr/bin/python
import os, sys, numpy

# create DEM
print "\ncreating DEM..."
# first some parameters used in the DEM and the grid creation
nx   = 20
ny   = 4
xlim = 500000. # in metres
ylim =  20000. # in metres
DEM_dir   = "FISOC_mesh1a_DEM"
surfFile  = DEM_dir+"/surf.xyz"
bedFile   = DEM_dir+"/bed.xyz"
gridName2D = "FISOC_mesh1a_2D"
gridName3D = "FISOC_mesh1a_3D"

if DEM_dir not in os.listdir('.'):
    ret = os.mkdir(DEM_dir)
# the +1 is because nx and ny refer to the number of elements, 
# but here we are writing values at nodes
xAxis = numpy.linspace(0,xlim,nx+1)
yAxis = numpy.linspace(0,ylim,ny+1)
bedrock  = numpy.zeros([nx+1, ny+1])
surface  = numpy.zeros([nx+1, ny+1])
surface[:,:] = 20000.
surf_out = open(surfFile,'w')
bed_out  = open(bedFile,'w')
for x in range(nx+1):
    for y in range(ny+1):
        surf_out.write(str(xAxis[x])+" "+str(yAxis[y])+" "+str(surface[x,y])+"\n")
        bed_out.write(str(xAxis[x])+" "+str(yAxis[y])+" "+str(bedrock[x,y])+"\n")
surf_out.close()
bed_out.close()

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
cmd = "ExtrudeMesh "+gridName2D+" "+gridName3D+" 6 10. 1 0 0 0 1 "+DEM_dir+" 1000. 1.5 -999999. > ExtrudeMesh.log"
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

