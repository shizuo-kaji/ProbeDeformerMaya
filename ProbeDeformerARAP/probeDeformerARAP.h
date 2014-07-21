#pragma once

#pragma comment(linker, "/export:initializePlugin /export:uninitializePlugin")

#include <maya/MFnPlugin.h>

#include <numeric>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/MatrixFunctions>

#include "affinelib.h"

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

using namespace Eigen;

class probeDeformerARAPNode : public MPxDeformerNode
{
public:
    probeDeformerARAPNode(): numTet(0), numPts(0), numPrb(0), transWeight(0.0)  {};
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
    static MObject      aTransWeight;
    static MObject      aConstraintWeight;
	static MObject		aRotationConsistency;
	static MObject		aFrechetSum;
    static MObject      aNormExponent;
    
private:
    void readMatrixArray(MArrayDataHandle& handle, std::vector<Matrix4d>& m);
    void tetMatrixC(const MPointArray& p, const MIntArray& triangles, std::vector<Matrix4d>& m, std::vector<Vector3d>& tetCenter);
	void arapHI(const std::vector<Matrix4d>& PI, const MIntArray& triangles);
	void arapG(const std::vector< Matrix4d>& At, const std::vector<Matrix4d>& PI,
                  const MIntArray& triangles, const std::vector<Matrix4d>& aff, MatrixXd& G);
	std::vector<Vector3d> prevNs;
	std::vector<double> prevThetas;
	std::vector<Matrix4d> PI;
	std::vector<Vector3d> tetCenter;     // barycenter of tetrahedra
    std::vector<Vector3d> probeCenter;
    double transWeight;
    double constraintWeight;
    double normExponent;
    bool worldMode;
//    SimplicialLDLT<SpMat> solver;
    SparseLU<SpMat> solver;
    SpMat F;                // ARAP constraint matrix
	MIntArray triangles;
    MPointArray pts;
	int numPts;
    int numTet;
    int numPrb;
    std::vector<int> constraintTet;
    std::vector<RowVector4d> constraintVector;
    std::vector<double> sidist;
    std::vector< std::vector<double> > idist;
};