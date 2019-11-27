
from paraview.simple import *
from vtk.util.numpy_support import vtk_to_numpy
import glob
import numpy as np
#import matplotlib.pyplot as plt 

GeomID_us = 103
GeomID_ls = 102

vtu_path = "/media/sf_VBshare/FISOC_Ex5_bil2/"

#-----------------------------------------------------------------------------------
def access_pvtu(pvtu_file):
    "Uses paraview.simple functionality to get data and coords from pvtu"
    
    # using paraview.simple (functionality overlap with vtk) to extract data and coords
    vtuHandle = XMLPartitionedUnstructuredGridReader(FileName=pvtu_file)
    vtuHandle.UpdatePipeline()
    
    # create a new 'Threshold' to access lower surface.   This is not completely robust.
    # (how to use find data in paraview.simple?)
    threshold_ls = Threshold(Input=vtuHandle)
    threshold_ls.Scalars = ['POINTS', 'height']
    threshold_ls.ThresholdRange = [0.0, 0.01]
    threshold_ls.UpdatePipeline()
    ls_data = servermanager.Fetch(threshold_ls)
    ls_PointData = ls_data.GetPointData()
    print pvtu_file
    ls_CellData = ls_data.GetCellData()
    ls_coords = vtk_to_numpy(ls_data.GetPoints().GetData())
    ls_coords_xy = ls_coords[:,0:2]
    
    # similar for upper surface
    threshold_us = Threshold(Input=vtuHandle)
    threshold_us.Scalars = ['POINTS', 'depth']
    threshold_us.ThresholdRange = [0.0, 0.01]
    threshold_us.UpdatePipeline()
    us_data = servermanager.Fetch(threshold_us)
    us_PointData = us_data.GetPointData()
    us_CellData = us_data.GetCellData()
    us_coords = vtk_to_numpy(us_data.GetPoints().GetData())
    us_coords_xy = us_coords[:,0:2]

    Delete(threshold_ls)
    Delete(threshold_us)
    Delete(vtuHandle)
    del threshold_ls
    del threshold_us
    del vtuHandle

    return ls_data, ls_CellData, ls_PointData, ls_coords_xy, us_data, us_CellData, us_PointData, us_coords_xy


if __name__ == "__main__":
    
    filesElmer=glob.glob(vtu_path+'*.pvtu')
    filesElmer.sort()
    
    graOutFile = open("graOverTime.asc","w+")
    graArr = np.zeros(len(filesElmer))

    for file_name in filesElmer:
        ls_data, ls_CellData, ls_PointData, ls_coords_xy, us_data, us_CellData, us_PointData, us_coords_xy = access_pvtu(file_name)
        GeometryIDS  = vtk_to_numpy(ls_CellData.GetArray(0))
        groundedMask = vtk_to_numpy(ls_PointData.GetArray('groundedmask'))
        gra=0.0
        for i in range(GeometryIDS.shape[0]):           # loop over elements
            celda1=ls_data.GetCell(i)                   # "Cell" is paraview name for Elmer "element"
            ids=celda1.GetPointIds()
            grounded = True
            for j in np.arange(ids.GetNumberOfIds()):   # loop over nodes for this element
                if (groundedMask[ids.GetId(j)] < 0.0):  # a node with groundedmask < 1 is floating
                    grounded = False
            if grounded:
                gra = gra + celda1.ComputeArea()        # add current element area only if grounded
 
        graOutFile.write('%.10e'%gra)
        graOutFile.write("\n")
    graOutFile.close()
        

