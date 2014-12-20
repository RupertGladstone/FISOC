#!/usr/bin/python

import sys, numpy
from netCDF4 import Dataset
import time

allElementTypes = {
    101 : (1,"nodal element"),
    202 : (2,"line segment"),
    203 : (2,"line segment quadratic"),
    204 : (2,"line segment cubic"),
    303 : (3,"triangle linear"),
    306 : (3,"triangle quadratic"),
    310 : (3,"triangle cubic"),
    404 : (4,"quadrilateral bilinear"),
    408 : (4,"quadrilateral quadratic"),
    409 : (4,"quadrilateral quadratic"),
    412 : (4,"quadrilateral cubic"),
    504 : (5,"tetrahedron linear"),
    510 : (5,"tetrahedron quadratic"),
    605 : (6,"pyramid linear"),
    613 : (6,"pyramid quadratic"),
    706 : (7,"wedges linear"),
    715 : (7,"wedges quadratic"),
    808 : (8,"hexahedron trilinear"),
    820 : (8,"hexahedron quadratic"),
    827 : (8,"hexahedron quadratice")
    }

validElementTypes = (303,404)

# we may need to write out the mapping for node:partition in case we cant store this info in ESMF format.
# in any case this code needs to be extended at some point to cope with partitioned meshes.

# Check element connectivity (node ordering) in elmer is consistent with esmf, i.e. counter clockwise.
# Node ordering looks fine, see ElmerSolver manual, appendix D

#--------------------------------------------------------------------------------------------
def readElmerMesh(meshFileName): 
    "Read an Elmer formatted mesh into python numpy arrays"
#--------------------------------------------------------------------------------------------

    headerFileName   = meshFileName+'/mesh.header'
    boundaryFileName = meshFileName+'/mesh.boundary'
    elementsFileName = meshFileName+'/mesh.elements'
    nodesFileName    = meshFileName+'/mesh.nodes'

    (header, elementTypes) = readElmerHeader(headerFileName)

#    boundary = readElmerBoundary(boundaryFileName)
    elements = readElmerElements(elementsFileName)        # get element -  node mapping
    validNodeIds = elmerGetValidNodeIds(elements)         # we only want the nodes corresponding to the elements we have kept
    nodes    = readElmerNodes(nodesFileName,validNodeIds) # get node coords

    return (header, elements, nodes, elementTypes)

#--------------------------------------------------------------------------------------------


def elmerGetValidNodeIds(elements):
#--------------------------------------------------------------------------------------------
    uniqueNodeIds = set()
    for elementId in elements.keys():
        uniqueNodeIds.update(elements[elementId])

    return uniqueNodeIds

#--------------------------------------------------------------------------------------------



def readElmerNodes(nodesFileName,validNodeIds):
#--------------------------------------------------------------------------------------------
    nodes = {}

    fMeshNodes = open(nodesFileName,'r')
    lines      = fMeshNodes.readlines()
    fMeshNodes.close()

    for line in lines:
        lineList    = line.split()
        nodeId      = int(lineList[0])
        partitionId = int(lineList[1])
        if nodeId in validNodeIds:
            nodes[nodeId] = map(float,lineList[2:5]) # 3 coords (z will be zero for 3d meshes) for this node id

    return nodes
#--------------------------------------------------------------------------------------------


def readElmerElements(elementsFileName):
#--------------------------------------------------------------------------------------------
    elements = {}

    fMeshElements = open(elementsFileName,'r')
    lines         = fMeshElements.readlines()
    fMeshElements.close()

    for line in lines:
        lineList    = line.split()
        elementId   = int(lineList[0])
        bodyId      = int(lineList[1])
        elementType = int(lineList[2])
        if (elementType in validElementTypes): 
            nodesPerElement = allElementTypes[elementType][0]
            elements[elementId] = map(int,lineList[3:3+nodesPerElement]) # list of node ids for this element id

    return elements
#--------------------------------------------------------------------------------------------


def readElmerHeader(headerFileName):
#--------------------------------------------------------------------------------------------

    elementTypes = set()
    header = {}
    header["numElements"] = {}

    fMeshHeader  = open(headerFileName,'r')
    lines        = fMeshHeader.readlines()
    fMeshHeader.close()

    for line in lines:
        lineList = line.split()
        if (len(lineList)==1):
            header["numTypes"]  = int(lineList[0])
        elif (len(lineList)==2):
            elementType         = int(lineList[0])
            numElements         = int(lineList[1])
            if (elementType not in allElementTypes.keys()):
                msg = "I don't recognise this element type: "+elementType
                raise ValueError(msg)
            if elementType in validElementTypes:
                elementTypes.update([elementType])
                header["numElements"][elementType] = numElements
        elif (len(lineList)==3):
            numNodes            = int(lineList[0])
            numBodyElements     = int(lineList[1])
            numBoundaryElements = int(lineList[2])
        else:
            msg = "I can't interpret this line in the mesh.header file: "+line
            raise ValueError(msg)


#header["numQuadrilaterals"] 
#header["numTriangles"] 

#    boundaryElementType = int(lines[2].split()[0])
#    bodyElementType     = int(lines[3].split()[0])

    return header, elementTypes
#--------------------------------------------------------------------------------------------


def convertMesh(elmerMesh):
#--------------------------------------------------------------------------------------------
    (eHeader, eElements, eNodes, eElementTypes) = elmerMesh # "e" for elmer

    eElementTypes = list(eElementTypes)

    maxNumNodes = 0
    for elementType in eElementTypes:
        maxNumNodes = max(maxNumNodes,allElementTypes[elementType][0])
    ESMF_maxNodePElement = maxNumNodes

    ESMF_nodeCount = len(eNodes)
    ESMF_elementCount = len(eElements)

    ESMF_coordDim = 2 # should this perhaps not be hard coded?  But I think we only need to used the bottom ice surface for exchange of variables in ESMF, at least to start with.

    # create the data using proper arrays.  Now the node and element Ids from the lists become array indices.

    ESMF_nodeCoords = numpy.zeros([len(eNodes), 2])
    for nodeId in eNodes.keys():
        ESMF_nodeCoords[nodeId-1,:] = eNodes[nodeId][0:2] # num dimensions hard coded to 2.  thus nodes will have 2 coords each

    ESMF_elementConn    = numpy.zeros([len(eElements), ESMF_maxNodePElement])
    ESMF_elementConn[:] = -1.0
    ESMF_numElementConn = numpy.zeros([len(eElements)])
    for elementId in eElements.keys():
        ESMF_numElementConn[elementId-1] = len(eElements[elementId])
        ESMF_elementConn[elementId-1,:]  = numpy.array(eElements[elementId])

    return ESMF_maxNodePElement, ESMF_nodeCount, ESMF_elementCount, ESMF_coordDim, ESMF_nodeCoords, ESMF_elementConn, ESMF_numElementConn


#--------------------------------------------------------------------------------------------


def writeESMFmesh(ESMFmesh,ESMFmeshFileName):
#--------------------------------------------------------------------------------------------

    (ESMF_maxNodePElement, ESMF_nodeCount, ESMF_elementCount, ESMF_coordDim, ESMF_nodeCoords, ESMF_elementConn, ESMF_numElementConn) = ESMFmesh

    rootgrp = Dataset(ESMFmeshFileName, 'w', format='NETCDF4')

    nodeCount       = rootgrp.createDimension('nodeCount', ESMF_nodeCount)
    elementCount    = rootgrp.createDimension('elementCount', ESMF_elementCount)
    maxNodePElement = rootgrp.createDimension('maxNodePElement', ESMF_maxNodePElement)
    coordDim        = rootgrp.createDimension('coordDim',ESMF_coordDim)

    rootgrp.history = 'Created ' + time.ctime(time.time())
    rootgrp.gridType= 'unstructured'
#    rootgrp.version = '0.9'

    nodeCoords     = rootgrp.createVariable('nodeCoords','f8',('nodeCount','coordDim'))
    elementConn    = rootgrp.createVariable('elementConn','i4',('elementCount','maxNodePElement'),fill_value=-1)
    numElementConn = rootgrp.createVariable('numElementConn','i4',('elementCount'))
#    centerCoords   = rootgrp.createVariable('centerCoords','f8',('elementCount','coordDim')) # optional!
#    elementArea    = rootgrp.createVariable('elementArea','f8',('elementCount')) # optional!
#    elementMask    = rootgrp.createVariable('elementMask','i4',('elementCount'),fill_value=-9999) # optional!

    nodeCoords.units         = 'degrees' 
    elementConn.long_name    = 'Node Indices that define the element connectivity'
    numElementConn.long_name = 'Number of nodes per element' 
#    centerCoords.units       = 'degrees' 
#    elementArea.units        = 'radians^2' 
#    elementArea.long_name    = 'area weights' 
 
    nodeCoords[:] =  ESMF_nodeCoords   
    numElementConn[:] = ESMF_numElementConn
    elementConn[:] = ESMF_elementConn

    rootgrp.close()

    rc = -1
    return rc

#--------------------------------------------------------------------------------------------

if __name__ == "__main__":

    if (len(sys.argv) < 2):
        raise ValueError("Missing command line arguments.  Please specify at least the Elmer mesh directory.")
    elif (len(sys.argv) == 2):
        ESMFmeshFileName = sys.argv[1]+'.nc'
        ESMFmeshFileName = ESMFmeshFileName.replace('/','')
    elif (len(sys.argv) == 3):
        ESMFmeshFileName = sys.argv[2]
    elif (len(sys.argv) > 3):
        raise ValueError("Too many command line arguments.  Please specify at most the Elmer mesh directory and the ESMF format mesh file to be written.")

    ElmerMeshDir = sys.argv[1]

    print 'I\'ll try to read Elmer mesh from ', ElmerMeshDir
    print 'I\'ll try to write converted ESMF mesh to ', ESMFmeshFileName

    elmerMesh = readElmerMesh(ElmerMeshDir)
    ESMFmesh = convertMesh(elmerMesh)
    rc = writeESMFmesh(ESMFmesh,ESMFmeshFileName)

    print "conversion complete"
