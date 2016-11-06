#pragma once

#pragma comment(linker, "/export:initializePlugin /export:uninitializePlugin")

#include <maya/MFnPlugin.h>

#include <numeric>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>

#include "affinelib.h"
#include "tetrise.h"
#include "MeshMaya.h"
#include "laplacian.h"
#include "blendAff.h"

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
    static MObject      aComputeWeight; // this attr will be dirtied when weight recomputation is needed
    static MObject      aInitMatrix;
    static MObject      aMatrix;
    static MObject      aBlendMode;
    static MObject      aTetMode;
    static MObject      aWorldMode;
	static MObject		aWeightMode;
	static MObject		aWeightCurveR;
	static MObject		aWeightCurveS;
	static MObject		aWeightCurveL;
	static MObject		aEffectRadius;     // global effect radius of probes
    static MObject      aTransWeight; // how much translation part affects in ARAP computation
    static MObject      aConstraintWeight; // global constraint weight of probes
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
    static MObject      aNormaliseWeight;
    static MObject      aAreaWeighted;
    static MObject      aNeighbourWeighting;
    
private:
    // variables
    BlendAff B;
    Laplacian mesh;
    std::vector<Vector3d> tetCenter; // center of tets
    std::vector<vertex> vertexList;   // mesh data
    std::vector<edge> edgeList;   // mesh data
    std::vector<int> faceList;   // mesh data
    std::vector<Vector3d> pts, new_pts;   // coordinates for mesh points
    std::vector< std::vector<double> > wr, ws, wl; // wr[j][i] is the weight of ith probe on j-th tet
    std::vector<int> closestPts, closestTet; // closestPts[i] is the index of pt closest to i-th probe
    short isError;  // to catch error
    int numPrb;  // number of probes
    std::vector<T> constraint;  // [row,col,value): row probe constraints col point with strength value
    std::vector< std::vector<double> > dist;    // dist[i][j] is the distance from j-th tet to i-th probe
    std::vector< std::vector<double> > distPts; // distPts[i][j] is the distance from i-th probe to j-th pt
    std::vector<Matrix4d> A,Q, blendedSE;  //temporary
    std::vector<Matrix3d> blendedR, blendedS;
    std::vector<Vector4d> blendedL;
    std::vector<double> tetEnergy;
    
};
