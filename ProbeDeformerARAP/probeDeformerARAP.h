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
    probeDeformerARAPNode(): numPrb(0), isError(0)  {};
    virtual MStatus deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex );
	virtual MStatus accessoryNodeSetup( MDagModifier& cmd );
    void    postConstructor();
    static  void*   creator();
    static  MStatus initialize();
 
    static MTypeId      id;
    static MString      nodeName;
    static MObject      aARAP;   // this attr will be dirtied when ARAP recomputation is needed
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
    static MObject      aProbeWeight;
    static MObject      aProbeConstraintRadius;
    
private:
    void readMatrixArray(MArrayDataHandle& handle, std::vector<Matrix4d>& m);
	void arapHI(const std::vector<Matrix4d>& PI, const std::vector<int>& tetList);
	void arapG(const std::vector< Matrix4d>& At, const std::vector<Matrix4d>& PI,
                  const std::vector<int>& tetList, const std::vector<Matrix4d>& Aff, MatrixXd& G);
    void visualise(MDataBlock& data, std::vector<double>& ptsColour);
    // variables
	std::vector<Vector3d> prevNs;   // for rotation consistency
	std::vector<double> prevThetas; // for rotation consistency
	std::vector<Matrix4d> PI;   // inverse of initial tet matrix
    std::vector<Vector3d> probeCenter; // initial location of probes
    std::vector<Vector3d> tetCenter; // center of tets
    std::vector<vertex> vertexList;   // mesh data
    std::vector<int> tetList;   // mesh data
    std::vector<edge> edgeList;   // mesh data
    std::vector<int> faceList;   // mesh data
    std::vector<Vector3d> pts;   // coordinates for mesh points
    std::vector<double> tetWeight; // tetWeight[j] is the stiffness of j-th tet
    double transWeight; // how much translation part affects in ARAP computation
    short isError;  // to catch error
    int numPrb;  // number of probes
    int dim; // total number of pts including the "ghost" added for forming tetrahedra
    std::vector< std::map<int,double> > constraint;  // if constraint[i] contains {j:x}, it means i-th probe constraints j-th point with weight x
    std::vector< std::vector<double> > dist;    // dist[j][i] is the distance from j-th tet to i-th probe
    SimplicialLDLT<SpMat> solver;   // ARAP solver
//    SimplicialCholesky<SpMat> solver;
//    SparseLU<SpMat> solver;
    SpMat F;                // ARAP constraint matrix
};