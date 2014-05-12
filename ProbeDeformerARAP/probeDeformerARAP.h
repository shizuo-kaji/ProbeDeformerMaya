//   Probe Deformer Maya Plugin
//   by Shizuo KAJI,     Nov. 2013
//
// requirements:  Maya,  Eigen library,   matrixlib

#pragma once

#pragma comment(linker, "/export:initializePlugin /export:uninitializePlugin")

#include <maya/MFnPlugin.h>

#include <numeric>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/MatrixFunctions>

#include "affinelib.h"

typedef Eigen::SparseMatrix<float> SpMat;
typedef Eigen::Triplet<double> T;

using namespace Eigen;

class probeDeformerARAPNode : public MPxDeformerNode
{
public:
    probeDeformerARAPNode(): numTet(0), numPts(0), numPrb(0), transWeight(0.0f)  {};
    virtual MStatus deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex );
	virtual MStatus accessoryNodeSetup( MDagModifier& cmd );
    static  void*   creator();
    static  MStatus initialize();
 
    static MTypeId      id;
    static MString      nodeName;
    static MObject      aInitMatrix;
    static MObject      aMatrix;
    static MObject      aBlendMode;
    static MObject      aTransMode;
	static MObject		aWeightMode;
	static MObject		aWeightCurveR;
	static MObject		aWeightCurveS;
	static MObject		aWeightCurveL;
	static MObject		aMaxDist;
    static MObject      aTransWeight;
	static MObject		aRotationConsistency;
	static MObject		aFrechetSum;

private:
    void readMatrixArray(MArrayDataHandle& handle, std::vector<Matrix4f>& m);
    void tetMatrixC(const MPointArray& p, const MIntArray& triangles, std::vector<Matrix4f>& m, std::vector<Vector3f>& tetCenter);
	void arapHI(const std::vector<Matrix4f>& PI, const MIntArray& triangles);
	void arapG(const std::vector< Matrix4f>& At, const std::vector<Matrix4f>& PI,
                  const MIntArray& triangles, const std::vector<Matrix4f>& aff, MatrixXf& G);
	std::vector<Vector3f> prevNs;
	std::vector<float> prevThetas;
	std::vector<Matrix4f> PI;
	std::vector<Vector3f> tetCenter;     // barycenter of tetrahedra
    std::vector<Vector3f> probeCenter;
    float transWeight;
    SparseLU<SpMat> solver;
	MIntArray triangles;
    MPointArray pts;
	unsigned int numPts;
    unsigned int numTet;
    int numPrb;
    std::vector<int> constraintTet;
    std::vector<RowVector4f> constraintVector;
    std::vector<float> sidist;
    std::vector< std::vector<float> > idist;
};