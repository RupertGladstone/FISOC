#!/usr/bin/python

from ROMS2Paraview_params import *
from netCDF4 import Dataset
import numpy

outCounter = 0
for nn in range(1,numInFiles+1):
    NC_inName = NC_inNameTemplate.replace("XXXX",str(nn).zfill(4))
    print "Reading netcdf file: ",NC_inName
    dataset = Dataset(NC_inName)
    print "Netcdf data set has format: ",dataset.file_format
    # print dataset.variables.keys() 

    timeLength = dataset.dimensions[timeVarName].size 
    print "Length of time dimension: ",timeLength

    #print dataset.dimensions.keys() 
    #s_w = dataset.dimensions['s_w']
    #ww = dataset.variables['w']
    #x_rho = dataset.variables['x_rho']
    #u_eastward = dataset.variables['u_eastward']
    melt = dataset.variables['m']
    
    # x and y dims on rho points
    xrlen = dataset.dimensions['xi_rho'].size
    yrlen = dataset.dimensions['eta_rho'].size
    print "x and y dimensions ", xrlen,yrlen

    # read coordinate information that doesn't evolve within this file
    xArr = dataset.variables['x_rho'][:,:]
    yArr = dataset.variables['y_rho'][:,:]
    #    print numpy.min(xArr), numpy.max(xArr)
    #    print numpy.min(yArr), numpy.max(yArr)
    #print xArr[1,:].shape
    #    print yArr[:,1]
    
    # loop over timesteps within this file
    skipCounter = skipLen
    for tt in range(timeLength):
        
        # write out netcdf file for this timestep, with correct coord info for Paraview
        if (skipCounter == skipLen):
            skipCounter = 1
            outCounter = outCounter + 1

            # create netcdf file (overwrites existing)
            NC_outName = NC_outNameTemplate.replace("XXXX",str(outCounter).zfill(4))
            print "Writing to file: ", NC_outName
            datasetOut = Dataset(NC_outName, 'w', format='NETCDF4')

            # extract time varying data from intput file
            meltArr = dataset.variables['m'][tt,:,:]

            # write out dimensions
            datasetOut.createDimension("x_rho", xrlen)
            datasetOut.createDimension("y_rho", yrlen)

            # write out coordinate variables
            x_rho_var = datasetOut.createVariable("x_rho", "float32",dimensions=("x_rho"))
            y_rho_var = datasetOut.createVariable("y_rho", "float32",dimensions=("y_rho"))
            datasetOut.variables['x_rho'][:] = xArr[1,:]
            datasetOut.variables['y_rho'][:] = yArr[:,1]
            
            ## write out coordinate variables
            #x_rho_var = datasetOut.createVariable("x_rho", "float32",dimensions=("y_rho","x_rho"))
            #y_rho_var = datasetOut.createVariable("y_rho", "float32",dimensions=("y_rho","x_rho"))
            #datasetOut.variables['x_rho'][:] = xArr[:,:]
            #datasetOut.variables['y_rho'][:] = yArr[:,:]
            
            # write out the fields we need as netcdf variables
            y_rho_var = datasetOut.createVariable("melt", "float32",dimensions=("y_rho", "x_rho"))
            datasetOut.variables['melt'][:] = melt[tt,:,:]

            datasetOut.close()

        else:
            skipCounter = skipCounter + 1
