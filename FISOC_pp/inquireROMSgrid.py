#!/usr/bin/python

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl

ROMS_grp  = Dataset("../ocean_his.nc", "r")

#print ROMS_grp.variables

ROMS_bmb = ROMS_grp.variables["m"]
xl = ROMS_grp.variables["xl"]
el = ROMS_grp.variables["el"]

print "ROMS domain length in XI direction ::: ",xl[:]/1000.," km" 
print "ROMS domain length in ETA direction ::: ",el[:]/1000.," km" 

ROMS_grp.close()


