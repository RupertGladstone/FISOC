#!/usr/bin/python

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np



ROMS_grp  = Dataset("../examples/Ex5_BensBox_gl/ocean_his.nc", "r")

#print ROMS_grp.variables

melt = ROMS_grp.variables["m"]
draft = ROMS_grp.variables["draft"]
time = ROMS_grp.variables["ocean_time"]

#print draft[1,:,:]


locI = 18; locJ = 4
print melt.shape
for t in range(time.size-1):
#    print sum(np.where(melt[t,:,:] > 1000.,1,0))
    try:
        draftDiff = draft[t+1,locI,locJ]-draft[t,locI,locJ]
        meltHere  = melt[t,locI,locJ]
    except:
        print t
        draftDiff = np.zeros(2)

#    if ( draftDiff != 0.0 ):
    if ( draftDiff.sum() != 0.0 ):
        print t,draft[t,locI,locJ],draftDiff.sum(),melt[t,locI,locJ]
#        print draftDiff
#print t
