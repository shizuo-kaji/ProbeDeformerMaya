#pragma once

#pragma comment(linker, "/export:initializePlugin /export:uninitializePlugin")

#include <maya/MFnPlugin.h>

#include <numeric>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Sparse>

#include "affinelib.h"
#include "tetrise.h"
#include "MeshMaya.h"
#include "ARAP.h"

using namespace Eigen;

typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;


class probeDeformerNode : public MPxDeformerNode
{
public:
    probeDeformerNode(): numPrb(0), numPts(0) {};
    virtual MStatus deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex );
	virtual MStatus accessoryNodeSetup( MDagModifier& cmd );
    static  void*   creator();
    static  MStatus initialize();
    void    postConstructor();
 
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
	static MObject		aEffectRadius;
    static MObject      aNormaliseWeight;
	static MObject		aRotationConsistency;
	static MObject		aFrechetSum;
    static MObject      aNormExponent;
    static MObject      aVisualisationMode;
    static MObject      aProbeWeight;
    static MObject      aComputeWeight;
    static MObject      aVisualisationMultiplier;
    static MObject      aAreaWeighted;

private:
	std::vector<Matrix4d> logSE;
	std::vector<Matrix3d> logR;
    std::vector< std::vector<double> > wr,ws,wl;
    int numPrb, numPts;
};
