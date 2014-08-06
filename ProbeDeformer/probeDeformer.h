#pragma once

#pragma comment(linker, "/export:initializePlugin /export:uninitializePlugin")

#include <maya/MFnPlugin.h>

#include <numeric>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Sparse>

#include "affinelib.h"
#include "tetrise.h"
#include "MeshMaya.h"

using namespace Eigen;

typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;


class probeDeformerNode : public MPxDeformerNode
{
public:
    probeDeformerNode(): numPrb(0) {};
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
	static MObject		aMaxDist;
	static MObject		aRotationConsistency;
	static MObject		aFrechetSum;
    static MObject      aNormExponent;
    static MObject      aVisualisationMode;
    static MObject      aProbeWeight;
    static MObject      aComputeWeight;
    static MObject      aVisualisationMultiplier;

private:
    void harmonicWeight(const std::vector<double>& probeWeight, const std::vector<int>& faceList,
                        const std::vector<Vector3d>& pts, const std::vector< std::vector<double> >& dist);

	std::vector<Vector3d> prevNs;
	std::vector<double> prevThetas;
    std::vector< std::vector<double> > wr,ws,wl;
    int numPrb;
};