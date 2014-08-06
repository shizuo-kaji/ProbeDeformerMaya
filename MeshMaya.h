/**
 * @file MeshMaya.h
 * @brief Mesh manipulation functions for Maya
 * @section LICENSE The MIT License
 * @section requirements:  Eigen library,   Maya
 * @version 0.10
 * @date  Aug. 2014
 * @author Shizuo KAJI
 */

#pragma once


//#include "mayaHeaders.h"
//#include "tetrise.h"

using namespace Tetrise;

// get mesh data
int getMeshData(MDataBlock& data, MObject& input, MObject& inputGeom, unsigned int mIndex,
                short tetMode, const std::vector<Vector3d>& pts, std::vector<int>& tetList,
                std::vector<int>& faceList, std::vector<edge>& edgeList,
                std::vector<vertex>& vertexList, std::vector<Matrix4d>& tetMat){
    // returns total number of pts including ghost ones
    // read mesh data
    int numPts = (int) pts.size();
    MStatus status;
    MArrayDataHandle hInput = data.outputArrayValue( input, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    status = hInput.jumpToElement( mIndex );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MObject oInputGeom = hInput.outputValue().child( inputGeom ).asMesh();
    MFnMesh inputMesh(oInputGeom);
    // face list
    MIntArray count, triangles;
    inputMesh.getTriangles( count, triangles );
    int numFaces = triangles.length()/3;
    faceList.resize(3*numFaces);
    for(int i=0;i<3*numFaces;i++){
        faceList[i]=triangles[i];
    }
    // vertex list
    MItMeshPolygon iter(oInputGeom);
    MIntArray faceVertices;
    vertexList.resize(numPts);
    for(int i=0;i<numPts;i++){
        vertexList[i].index = i;
        vertexList[i].connectedTriangles.clear();
    }
    for( ; ! iter.isDone(); iter.next()){
        status = iter.getVertices(faceVertices);
        int count = (int) faceVertices.length();
        for(int j=0;j<count;j++){
            vertexList[faceVertices[j]].connectedTriangles.push_back(faceVertices[(j+1)%count]);
            vertexList[faceVertices[j]].connectedTriangles.push_back(faceVertices[(j+count-1)%count]);
        }
    }
    makeEdgeList(faceList, edgeList);
    int dim=makeTetList(tetMode, numPts, faceList, edgeList, vertexList, tetList);
    makeTetMatrix(tetMode, pts, tetList, faceList, edgeList, vertexList, tetMat);
    return dim;
}

// visualise vertex assigned values
void visualise(MDataBlock& data, MObject& outputGeom, std::vector<double>& ptsColour){
    // load target mesh output
    MStatus status;
    MArrayDataHandle outputArray = data.outputArrayValue(outputGeom , &status );
    MDataHandle hOutput = outputArray.inputValue(&status);
    MFnMesh outMesh(hOutput.data());
    MColorArray Colour;
    MIntArray Index;
    // set vertex colour
    for(int i=0;i<ptsColour.size();i++){
        Colour.append( ptsColour[i] ,0,0);
        Index.append(i);
    }
    outMesh.setVertexColors(Colour, Index);
}

// read array of matrix attributes and convert them to Eigen matrices
void readMatrixArray(MArrayDataHandle& handle, std::vector<Matrix4d>& m){
    int numPrb=handle.elementCount();
    m.resize(numPrb);
    MMatrix mat;
    for(int i=0;i<numPrb;i++){
        handle.jumpToArrayElement(i);
        mat=handle.inputValue().asMatrix();
        m[i] << mat(0,0), mat(0,1), mat(0,2), mat(0,3),
        mat(1,0), mat(1,1), mat(1,2), mat(1,3),
        mat(2,0), mat(2,1), mat(2,2), mat(2,3),
        mat(3,0), mat(3,1), mat(3,2), mat(3,3);
    }
}

