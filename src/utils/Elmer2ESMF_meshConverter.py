#!/usr/bin/python

import sys, numpy
from netCDF4 import Dataset
import time

elementTypes = {
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

    header   = readElmerHeader(headerFileName)

    print len(header), len(header["numElements"])
#    boundary = readElmerBoundary(boundaryFileName)
    elements = readElmerElements(elementsFileName)  # get element -  node mapping
#    print elements[1]
    nodes    = readElmerNodes(nodesFileName)        # get node coords
#    print nodes[2]

    return (header, elements, nodes)

#--------------------------------------------------------------------------------------------


def readElmerNodes(nodesFileName):
#--------------------------------------------------------------------------------------------
    nodes = {}

    fMeshNodes = open(nodesFileName,'r')
    lines         = fMeshNodes.readlines()
    fMeshNodes.close()

    for line in lines:
        lineList    = line.split()
        nodeId      = int(lineList[0])
        partitionId = int(lineList[1])
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
            nodesPerElement = elementTypes[elementType][0]
            elements[elementId] = map(int,lineList[3:3+nodesPerElement]) # list of node ids for this element id

    return elements
#--------------------------------------------------------------------------------------------


def readElmerHeader(headerFileName):
#--------------------------------------------------------------------------------------------

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
            if (elementType not in elementTypes.keys()):
                msg = "I don't recognise this element type: "+elementType
                raise ValueError(msg)
            if elementType in validElementTypes:
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

    return header
#--------------------------------------------------------------------------------------------

if __name__ == "__main__":

    if (len(sys.argv) < 2):
        raise ValueError("Missing command line arguments.  Please specify at least the Elmer mesh directory.")
    elif (len(sys.argv) == 2):
        ESMFmeshFileName = sys.argv[1]+'.nc'
    elif (len(sys.argv) == 3):
        ESMFmeshFileName = sys.argv[2]
    elif (len(sys.argv) > 3):
        raise ValueError("Too many command line arguments.  Please specify at most the Elmer mesh directory and the ESMF format mesh file to be written.")

    ElmerMeshDir = sys.argv[1]

    print ElmerMeshDir, ESMFmeshFileName

    mesh = readElmerMesh(ElmerMeshDir)
