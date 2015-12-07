#!/usr/bin/python

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl

# Compare ROMS direct outputs against the same variables outputted from FISOC 
# (there should be no difference if the FISOC reformatting of vars to be ESMF 
# friendly is working)

timestep = 1000

ROMS_grp  = Dataset("../ocean_his.nc", "r")

FISOC_ncdfName = "../FISOC_OM_exp_t"+str(timestep/10)+".nc"
FISOC_grp = Dataset(FISOC_ncdfName, "r")

#rootgrp = Dataset("FISOC_OM_exp_t138.nc", "r", format="NETCDF4")
#rootgrp = Dataset("FISOC_OM_exp_t128.nc", "r", format="NETCDF4")
#print rootgrp.data_model
#print rootgrp.groups
#print rootgrp
#print rootgrp.dimensions

#for children in walktree(rootgrp):
#    for child in children:
#        print child

print ROMS_grp.variables

FISOC_dBdt_l0 = FISOC_grp.variables["OM_dBdt_l0"]
FISOC_temperature_l0 = FISOC_grp.variables["OM_temperature_l0"]

ROMS_dBdt_l0 = ROMS_grp.variables["m"]
ROMS_dBdt_l0 = ROMS_dBdt_l0[timestep,:,:]

ROMS_z_l0 = ROMS_grp.variables["draft"]
ROMS_z_l0 = ROMS_z_l0[timestep,:,:]

ROMS_temperature_l0 = ROMS_grp.variables["Tb"]
ROMS_temperature_l0 = ROMS_temperature_l0[timestep,:,:]

ROMS_temperature_lz = ROMS_grp.variables["temp"]
print ROMS_temperature_lz.shape
ROMS_temperature_lz10 = ROMS_temperature_lz[timestep,10,:,:]
ROMS_temperature_lz0 = ROMS_temperature_lz[timestep,0,:,:]


plt.subplot(1, 7, 1)
plt.imshow(ROMS_dBdt_l0,vmin=0.0,vmax=0.000001)
plt.colorbar()
plt.axis('off')  

plt.subplot(1, 7, 2)
plt.imshow(FISOC_dBdt_l0,vmin=0.0,vmax=0.000001)
plt.colorbar()
plt.axis('off')  

plt.subplot(1, 7, 3)
plt.imshow(ROMS_temperature_lz10,vmin=-2.3,vmax=-1.8)
plt.colorbar()
plt.axis('off')  

plt.subplot(1, 7, 4)
plt.imshow(ROMS_temperature_lz0,vmin=-2.3,vmax=-1.8)
plt.colorbar()
plt.axis('off')  

plt.subplot(1, 7, 5)
plt.imshow(ROMS_temperature_l0,vmin=-2.3,vmax=-1.8)
plt.colorbar()
plt.axis('off')  

plt.subplot(1, 7, 6)
plt.imshow(FISOC_temperature_l0,vmin=-2.3,vmax=-1.8)
plt.colorbar()
plt.axis('off')  

plt.subplot(1, 7, 7)
plt.imshow(ROMS_z_l0)
plt.colorbar()
plt.axis('off')  

plt.savefig("FISOC_testing.pdf")
plt.show()


ROMS_grp.close()
FISOC_grp.close()

def walktree(top):
    values = top.groups.values()
    yield values
    for value in top.groups.values():
        for children in walktree(value):
            yield children

