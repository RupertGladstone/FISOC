#!/usr/bin/python

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl


ROMS_grp  = Dataset("../examples/Ex4_BensBox/ocean_his.nc", "r")

#print ROMS_grp.variables

draft = ROMS_grp.variables["draft"]
time = ROMS_grp.variables["ocean_time"]

#print draft[1,:,:]

for t in range(time.size-1):
    try:
        draftDiff = draft[t+1,3,3]-draft[t,3,3]        
    except:
        print t

    if ( draftDiff != 0.0 ):
#    if ( draftDiff.sum() != 0.0 ):
        print t,draft[t,3,3],draft[t+1,3,3],draftDiff.sum()
#        print draftDiff
print t
