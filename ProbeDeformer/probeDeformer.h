#pragma once

#pragma comment(linker, "/export:initializePlugin /export:uninitializePlugin")

#include <maya/MFnPlugin.h>

#include <numeric>
#include <unsupported/Eigen/MatrixFunctions>

#include "../affinelib.h"

using namespace Eigen;

class probeDeformerNode : public MPxDeformerNode
{
public:
    probeDeformerNode()  {};
    virtual MStatus deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex );
	virtual MStatus accessoryNodeSetup( MDagModifier& cmd );
    static  void*   creator();
    static  MStatus initialize();
 
    static MTypeId      id;
    static MString      nodeName;
    static MObject      aInitMatrix;
    static MObject      aMatrix;
    static MObject      aBlendMode;
    static MObject      aWorldMode;
	static MObject		aWeightMode;
	static MObject		aWeightCurveR;
	static MObject		aWeightCurveS;
	static MObject		aWeightCurveL;
	static MObject		aMaxDist;
	static MObject		aRotationConsistency;
	static MObject		aFrechetSum;
    static MObject      aNormExponent;

private:
    void readMatrixArray(MArrayDataHandle& handle, std::vector<Matrix4d>& m);

	std::vector<Vector3d> prevNs;
	std::vector<double> prevThetas;
};