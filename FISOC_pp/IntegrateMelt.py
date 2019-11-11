
from paraview.simple import *
from vtk.util.numpy_support import vtk_to_numpy
import glob
import numpy as np
#import matplotlib.pyplot as plt 

GeomID_us = 103
GeomID_ls = 102

vtu_path = "/media/sf_VBshare/tmp/"

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

    meltOutFile = open("meltOverTime.asc","w+")

    meltArr = np.zeros(len(filesElmer))
    counter = -1
    for file_name in filesElmer:

        counter = counter + 1
        ls_data, ls_CellData, ls_PointData, ls_coords_xy, us_data, us_CellData, us_PointData, us_coords_xy = access_pvtu(file_name)
        GeometryIDS=vtk_to_numpy(ls_CellData.GetArray(0))
        
        bmb = vtk_to_numpy(ls_PointData.GetArray('meltrate'))
        #GL = vtk_to_numpy(ls_PointData.GetArray('groundedmask'))
        
        BMBInteg=0
        indexGEO=np.where(GeometryIDS==GeomID_ls) #bottom layer
        for i in indexGEO[0]:
            celda1=ls_data.GetCell(i)
            ids=celda1.GetPointIds()
            if ids.GetNumberOfIds()==3:
                mean=0
                #            meanGL=0
                # compute the mean bmb of three points
                for j in np.arange(3):
                    #                meanGL=meanGL+GL[ids.GetId(j)]
                    mean=mean+bmb[ids.GetId(j)]/3
                #        if meanGL==0:
            BMBLocal=mean*celda1.ComputeArea() # BMB*area
            BMBInteg=BMBInteg+BMBLocal
        print "Total BMB flux ",BMBInteg 
        meltArr[counter] = BMBInteg

        #Grounded Area
        AreaInteg=0
        for i in indexGEO[0]:
            celda1=ls_data.GetCell(i)
            ids=celda1.GetPointIds()
            if ids.GetNumberOfIds()==3:
#                mean=0
#                for j in np.arange(3):
#                    mean=mean+GL[ids.GetId(j)]/3
#                if mean>0:
                AreaLocal=celda1.ComputeArea()
#                else:
#                    AreaLocal=0.             
            AreaInteg=AreaInteg+AreaLocal
            
        #Floating ice area
        AreaWholeInteg=0
        for i in indexGEO[0]:
            celda1=ls_data.GetCell(i)
            AreaLocal=celda1.ComputeArea()
            AreaWholeInteg=AreaWholeInteg+AreaLocal

        meltOutFile.write('%.10e'%BMBInteg)
        meltOutFile.write("\n")
    meltOutFile.close()
        

  
## line 1 points 
#x1 = [1,2,3] 
#y1 = [2,4,1] 
## plotting the line 1 points  
#plt.plot(x1, y1, label = "line 1") 
  
## line 2 points 
#x2 = [1,2,3] 
#y2 = [4,1,3] 
## plotting the line 2 points  
#plt.plot(x2, y2, label = "line 2") 
  
## naming the x axis 
#plt.xlabel('x - axis') 
## naming the y axis 
#plt.ylabel('y - axis') 
## giving a title to my graph 
#plt.title('Two lines on same graph!') 
  
## show a legend on the plot 
#plt.legend() 
  
## function to show the plot 
#plt.show() 
