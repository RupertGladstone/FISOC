from netCDF4 import Dataset
import numpy

# This operates on a netcdf file containing ROMS output after it has been processed by
# ROMS2para.m to set up the coords.
# The time series of ocean volume is written to text file.

inFile = "/media/sf_VBshare/ocean_his_bil.nc"

ib = 1 # ignore boundary: number of grid cells to ignore at edge of domain

outFile = inFile.replace(".nc","_vol.asc")
outHandle = open(outFile,"w+")

# access netcdf variables
dataset = Dataset(inFile)
bathy    = dataset.variables['h']
iceDraft = dataset.variables['draft']
xi_rho   = dataset.variables['xi_rho']
eta_rho  = dataset.variables['eta_rho']

# process some meta data
dx = xi_rho[1]  - xi_rho[0]  # assume grid cells are equally spaced so we can use only one value for dx
dy = eta_rho[1] - eta_rho[0] # assume grid cells are equally spaced so we can use only one value for dy
area = dx * dy
nTimes = iceDraft.shape[0] # use first array index to get number of timesteps, because iceDraft is time evolving.

# remove the boundary (ROMS has halo cells not matched by Elmer)
if ib>0:
    bathy    = bathy[ib:-ib,ib:-ib]
    iceDraft = iceDraft[:,ib:-ib,ib:-ib]
    
# calculate volume at each timestep
vol = numpy.zeros(nTimes)
for tt in range(nTimes):
    wct = bathy[:,:] + iceDraft[tt,:,:] # water column thickness. bathy is positive down, iceDraft is positive up, hence addition not subtraction.
    vol = numpy.sum(wct) * area
    print tt, vol
    outHandle.write('%.10e'%vol)
    outHandle.write("\n")

dataset.close()
outHandle.close()

#loop over time
#extract inner (remove border)
#take difference between top and bottom coords (or draft and h maybe?)
#numpy.sum()
