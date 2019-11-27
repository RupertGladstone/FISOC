from netCDF4 import Dataset
import numpy

# This operates on a netcdf file containing ROMS output after it has been processed by
# ROMS2para.m to set up the coords.
# The time series of ocean volume is written to text file.

#inFile = "/media/sf_VBshare/FISOC_Ex4_cp/ocean_his_select.nc"     # FX4_cp (nearest)
#inFile = "/media/sf_VBshare/FISOC_Ex4_interp/ocean_his_select.nc" # ER_FX4_interp (bilin)
#inFile = "/media/sf_VBshare/FISOC_Ex4_RO/ocean_his_select.nc"     # FX4_RO ROMS only
inFile = "/media/sf_VBshare/FISOC_Ex5_bil2/ocean_his_select.nc" # ER_FX4_interp (bilin)

ib = 1 # ignore boundary: number of grid cells to ignore at edge of domain

outFile = inFile.replace(".nc","_gra.asc") # gra is short for grounded area
outHandle = open(outFile,"w+")

# access netcdf variables
dataset    = Dataset(inFile)
wetdryMask = dataset.variables['wetdry_mask_rho']
xi_rho     = dataset.variables['xi_rho']
eta_rho    = dataset.variables['eta_rho']

# process some meta data
dx     = xi_rho[1]  - xi_rho[0]  # assume grid cells are equally spaced so we can use only one value for dx
dy     = eta_rho[1] - eta_rho[0] # assume grid cells are equally spaced so we can use only one value for dy
area   = dx * dy                 # area of each cell (we assume cells all have the same area)
nTimes = wetdryMask.shape[0]     # use first array index to get number of timesteps, because wetdryMask is time evolving.

# remove the boundary (ROMS has halo cells not matched by Elmer)
if ib>0:
    wetdryMask = wetdryMask[:,ib:-ib,ib:-ib]
    
print wetdryMask.shape, type(wetdryMask)

# invert the mask (we want the ones to be grounded cells and zeros floating)
wetdryMask = 1.0 - wetdryMask

# calculate grounded area at each timestep
gra = numpy.zeros(nTimes)
for tt in range(nTimes):
    gra = numpy.sum(wetdryMask[tt,:,:]) * area
    print tt, gra, numpy.sum(wetdryMask)
    outHandle.write('%.10e'%gra)
    outHandle.write("\n")

dataset.close()
outHandle.close()
