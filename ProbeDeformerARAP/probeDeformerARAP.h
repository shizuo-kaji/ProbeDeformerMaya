#pragma once

#pragma comment(linker, "/export:initializePlugin /export:uninitializePlugin")

#include <maya/MFnPlugin.h>

#include <numeric>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>

#include "affinelib.h"
#include "tetrise.h"

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

using namespace Eigen;

//deformer
class probeDeformerARAPNode : public MPxDeformerNode
{
public:
    probeDeformerARAPNode(): numPrb(0), tetMode(-1), isError(0)  {};
    virtual MStatus deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex );
	virtual MStatus accessoryNodeSetup( MDagModifier& cmd );
    void    postConstructor();
    static  void*   creator();
    static  MStatus initialize();
 
    static MTypeId      id;
    static MString      nodeName;
    static MObject      aInitMatrix;
    static MObject      aMatrix;
    static MObject      aBlendMode;
    static MObject      aTetMode;
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
    static MObject      aIteration;
    static MObject      aConstraintRadius;
    static MObject      aConstraintMode;
    static MObject      aVisualisationMode;
    static MObject      aVisualisationMultiplier;
    static MObject      aSupervisedMesh;
    static MObject      aStiffness;
    
private:
    void readMatrixArray(MArrayDataHandle& handle, std::vector<Matrix4d>& m);
	void arapHI(const std::vector<Matrix4d>& PI, const std::vector<int>& tetList);
	void arapG(const std::vector< Matrix4d>& At, const std::vector<Matrix4d>& PI,
                  const std::vector<int>& tetList, const std::vector<Matrix4d>& Aff, MatrixXd& G);
    void visualise(MDataBlock& data, std::vector<double>& ptsColour);
	std::vector<Vector3d> prevNs;
	std::vector<double> prevThetas;
	std::vector<Matrix4d> PI;
    std::vector<Vector3d> probeCenter;
    std::vector<Vector3d> tetCenter;
    std::vector<vertex> vertexList;
    std::vector<int> tetList;
    std::vector<edge> edgeList;
    std::vector<int> faceList;
    std::vector<Vector3d> pts;
    std::vector<double> tetWeight;
    double transWeight;
    double constraintWeight;
    double normExponent, constraintRadius;
    bool worldMode;
    short tetMode, constraintMode, stiffnessMode, isError;
    SimplicialLDLT<SpMat> solver;
//    SimplicialCholesky<SpMat> solver;
//    SparseLU<SpMat> solver;
    SpMat F;                // ARAP constraint matrix
	MIntArray triangles;
    int numPrb;
    int dim;
    std::vector< std::map<int,double> > constraint;
    std::vector< std::vector<double> > dist;
};