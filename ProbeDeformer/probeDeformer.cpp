/**
 * @file probeDeformer.cpp
 * @brief Probe Deformer plugin for Maya
 * @section LICENSE The MIT License
 * @section  requirements:  Eigen library, Maya
 * @version 0.14
 * @date  3/Nov/2013
 * @author Shizuo KAJI
 */

#include "StdAfx.h"
#include "probeDeformer.h"

using namespace Eigen;
using namespace AffineLib;
using namespace Tetrise;

// parametrisation mode
#define BM_SRL 0
#define BM_SES 1
#define BM_LOG3 3
#define BM_LOG4 4
#define BM_QSL 5
#define BM_AFF 10
#define BM_OFF -1

// weight mode
#define WM_INV_DISTANCE 0
#define WM_CUTOFF_DISTANCE 1
#define WM_DRAW 2
#define WM_HARMONIC 3

// visualisation mode
#define VM_OFF 0
#define VM_EFFECT 2
#define VM_STIFFNESS 4

MTypeId probeDeformerNode::id( 0x00000103 );
MString probeDeformerNode::nodeName( "probeDeformer" );
MObject probeDeformerNode::aMatrix;
MObject probeDeformerNode::aInitMatrix;
MObject probeDeformerNode::aWorldMode;
MObject probeDeformerNode::aBlendMode;
MObject probeDeformerNode::aWeightMode;
MObject probeDeformerNode::aWeightCurveR;
MObject probeDeformerNode::aWeightCurveS;
MObject probeDeformerNode::aWeightCurveL;
MObject probeDeformerNode::aMaxDist;
MObject probeDeformerNode::aRotationConsistency;
MObject probeDeformerNode::aFrechetSum;
MObject probeDeformerNode::aNormExponent;
MObject probeDeformerNode::aProbeWeight;
MObject probeDeformerNode::aVisualisationMode;
MObject probeDeformerNode::aComputeWeight;
MObject probeDeformerNode::aVisualisationMultiplier;

void* probeDeformerNode::creator() { return new probeDeformerNode; }
 
MStatus probeDeformerNode::deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex ){
	MObject thisNode = thisMObject();
    MStatus status;
    MThreadUtils::syncNumOpenMPThreads();    // for OpenMP
    bool worldMode = data.inputValue( aWorldMode ).asBool();
    short blendMode = data.inputValue( aBlendMode ).asShort();
    short weightMode = data.inputValue( aWeightMode ).asShort();
    short visualisationMode = data.inputValue( aVisualisationMode ).asShort();
    double maxDist = data.inputValue( aMaxDist ).asDouble();
    double normExponent = data.inputValue( aNormExponent ).asDouble();
	bool rotationCosistency = data.inputValue( aRotationConsistency ).asBool();
	bool frechetSum = data.inputValue( aFrechetSum ).asBool();
    double visualisationMultiplier = data.inputValue(aVisualisationMultiplier).asDouble();
    MArrayDataHandle hMatrixArray = data.inputArrayValue(aMatrix);
    MArrayDataHandle hInitMatrixArray = data.inputArrayValue(aInitMatrix);
    bool isNumProbeChanged = (numPrb != hMatrixArray.elementCount());
    numPrb = hMatrixArray.elementCount();
    if(numPrb != hInitMatrixArray.elementCount() || numPrb==0 || blendMode == BM_OFF)
        return MS::kSuccess;
	if( ! rotationCosistency || numPrb != prevNs.size()){
		prevThetas.clear();
		prevThetas.resize(numPrb, 0.0);
		prevNs.clear();
		prevNs.resize(numPrb, Vector3d::Zero());
	}
// setting transformation matrix
    std::vector<Matrix4d> initMatrix(numPrb), matrix(numPrb);
    readMatrixArray(hInitMatrixArray, initMatrix);
    readMatrixArray(hMatrixArray, matrix);
    std::vector<Matrix3d> logR(numPrb),R(numPrb),logS(numPrb),S(numPrb),logGL(numPrb);
    std::vector<Matrix4d> logSE(numPrb),SE(numPrb),logAff(numPrb),Aff(numPrb);
    std::vector<Vector3d> L(numPrb);
    std::vector<Vector4d> quat(numPrb);
    std::vector<Vector3d> probeCenter(numPrb);
    for( int i=0;i<numPrb;i++){
        Aff[i]=initMatrix[i].inverse()*matrix[i];
        probeCenter[i] << transPart(initMatrix[i]);
    }
    if(blendMode == BM_SRL || blendMode == BM_SES || blendMode == BM_QSL){
        for(int i=0;i<numPrb;i++){
            parametriseGL(Aff[i].block(0,0,3,3), logS[i] ,R[i]);
            L[i] = transPart(Aff[i]);
            if(blendMode == BM_SRL){
                logR[i]=logSOc(R[i], prevThetas[i], prevNs[i]);
            }else if(blendMode == BM_SES){
                SE[i]=pad(R[i], L[i]);
                logSE[i]=logSEc(SE[i], prevThetas[i], prevNs[i]);
            }else if(blendMode == BM_QSL){
                Quaternion<double> Q(R[i].transpose());
                quat[i] << Q.x(), Q.y(), Q.z(), Q.w();
                S[i]=expSym(logS[i]);
            }
        }
    }else if(blendMode == BM_LOG3){
        for(int i=0;i<numPrb;i++){
            logGL[i] = Aff[i].block(0,0,3,3).log();
            L[i] = transPart(Aff[i]);
        }
    }else if(blendMode == BM_LOG4){
        for(int i=0;i<numPrb;i++){
            logAff[i] = Aff[i].log();
        }
    }
    
// transform target vertices
    // get positions
    MPointArray Mpts;
    itGeo.allPositions(Mpts);
    int numPts = Mpts.length();
    if(worldMode){
        for(int j=0; j<numPts; j++ )
            Mpts[j] *= localToWorldMatrix;
    }
    std::vector<Vector3d> pts(numPts);
    for(int i=0;i<numPts;i++){
        pts[i] << Mpts[i].x, Mpts[i].y, Mpts[i].z;
    }
    
    // load per vertex weights
    std::vector<double> ptsWeight(numPts);
    for (int i=0; !itGeo.isDone(); itGeo.next()){
        ptsWeight[i++] = weightValue(data, mIndex, itGeo.index());
    }
    

    // weight computation
    if(!data.isClean(aComputeWeight) || isNumProbeChanged){
        status = data.setClean(aComputeWeight);
        // load probe weights
        std::vector<double> probeWeight(numPrb), probeRadius(numPrb);
        MArrayDataHandle handle = data.inputArrayValue(aProbeWeight);
        if(handle.elementCount() != numPrb){
            MGlobal::displayInfo("# of Probes and probeWeight are different");
            return MS::kFailure;
        }
        for(int i=0;i<numPrb;i++){
            handle.jumpToArrayElement(i);
            probeWeight[i]=handle.inputValue().asDouble();
            probeRadius[i] = probeWeight[i] * maxDist;
        }
        //
        wr.resize(numPts),ws.resize(numPts),wl.resize(numPts);
        std::vector< std::vector<double> > dist(numPts);
        for(int j=0; j<numPts; j++ ){
            wr[j].resize(numPrb);ws[j].resize(numPrb);wl[j].resize(numPrb);
            dist[j].resize(numPrb);
            for( int i=0; i<numPrb; i++){
                dist[j][i] = (pts[j]-probeCenter[i]).norm();
            }
        }
        if(weightMode == WM_INV_DISTANCE){
            for(int j=0; j<numPts; j++ ){
                double sum=0;
                std::vector<double> idist(numPrb);
                for( int i=0; i<numPrb; i++){
                    idist[i] = probeWeight[i]/pow(dist[j][i],normExponent);
                    sum += idist[i];
                }
                assert( sum > 0);
                for( int i=0; i<numPrb; i++){
                    wr[j][i] = ws[j][i] = wl[j][i] = idist[i]/sum;
                }
            }
        }else if(weightMode == WM_CUTOFF_DISTANCE){
            for(int j=0; j<numPts; j++ ){
                for( int i=0; i<numPrb; i++){
                    wr[j][i] = ws[j][i] = wl[j][i] = (dist[j][i] > probeRadius[i])
                    ? 0 : pow((probeRadius[i]-dist[j][i])/probeRadius[i],normExponent);
                }
            }
        }else if(weightMode == WM_DRAW){
            MRampAttribute rWeightCurveR( thisNode, aWeightCurveR, &status );
            MRampAttribute rWeightCurveS( thisNode, aWeightCurveS, &status );
            MRampAttribute rWeightCurveL( thisNode, aWeightCurveL, &status );
            float val;
            for(int j=0; j<numPts; j++ ){
                for( int i=0; i<numPrb; i++){
                    rWeightCurveR.getValueAtPosition(dist[j][i]/probeRadius[i], val );
                    wr[j][i] = val;
                    rWeightCurveS.getValueAtPosition(dist[j][i]/probeRadius[i], val );
                    ws[j][i] = val;
                    rWeightCurveL.getValueAtPosition(dist[j][i]/probeRadius[i], val );
                    wl[j][i] = val;
                }
            }
        }else if(weightMode == WM_HARMONIC){
            // face list
            MArrayDataHandle hInput = data.outputArrayValue( input, &status );
            status = hInput.jumpToElement( mIndex );
            MObject oInputGeom = hInput.outputValue().child( inputGeom ).asMesh();
            MFnMesh inputGeom(oInputGeom);
            MIntArray count, triangles;
            inputGeom.getTriangles( count, triangles );
            std::vector<int> faceList;
            faceList.resize(triangles.length());
            for(int i=0;i<triangles.length();i++){
                faceList[i]=triangles[i];
            }
            harmonicWeight(probeWeight,faceList,pts,dist);
        }
        // normalise
        for(int j=0;j<numPts;j++){
            double sum = std::accumulate(wr[j].begin(), wr[j].end(), 0.0);
            if (sum > 1){
                for (int i = 0; i < numPrb; i++){
                    wr[j][i] /= sum;
                }
            }
            sum = std::accumulate(ws[j].begin(), ws[j].end(), 0.0);
            if (sum > 1){
                for (int i = 0; i < numPrb; i++){
                    ws[j][i] /= sum;
                }
            }
            sum = std::accumulate(wl[j].begin(), wl[j].end(), 0.0);
            if (sum > 1){
                for (int i = 0; i < numPrb; i++){
                    wl[j][i] /= sum;
                }
            }
        }
    }
    
#pragma omp parallel for
    for(int j=0; j<numPts; j++ ){
        std::vector<double> wrr(numPrb),wss(numPrb),wll(numPrb);
        for(int i=0;i<numPrb;i++){
            wrr[i]=ptsWeight[j]*wr[j][i];
            wss[i]=ptsWeight[j]*ws[j][i];
            wll[i]=ptsWeight[j]*wl[j][i];
        }
        // blend matrix
        Matrix4d mat;
        if(blendMode == BM_SRL){
            Matrix3d RR,SS=expSym(blendMat(logS, wrr));
            Vector3d l=blendMat(L, wll);
            if(frechetSum){
                RR = frechetSO(R, wrr);
            }else{
                RR = expSO(blendMat(logR, wrr));
            }
            mat = pad(SS*RR, l);
        }else if(blendMode == BM_SES){
            Matrix4d RR;
            Matrix3d SS=expSym(blendMat(logS, wss));
            if(frechetSum){
                RR = frechetSE(SE, wrr);
            }else{
                RR = expSE(blendMat(logSE, wrr));
            }
            mat = pad(SS,Vector3d::Zero()) * RR;
        }else if(blendMode == BM_LOG3){
            Matrix3d RR=blendMat(logGL, wrr).exp();
            Vector3d l=blendMat(L, wll);
            mat = pad(RR, l);
        }else if(blendMode == BM_LOG4){
            mat=blendMat(logAff, wrr).exp();
        }else if(blendMode == BM_QSL){
            Vector4d q=blendQuat(quat,wrr);
            Vector3d l=blendMat(L, wll);
            Matrix3d SS=blendMatLin(S,wss);
            Quaternion<double> Q(q);
            Matrix3d RR = Q.matrix().transpose();
            mat = pad(SS*RR, l);
        }else if(blendMode == BM_AFF){
            mat = blendMatLin(Aff,wrr);
        }
        // apply matrix
        RowVector4d p = pad(pts[j]) * mat;
        Mpts[j].x = p[0];
        Mpts[j].y = p[1];
        Mpts[j].z = p[2];
        if(worldMode)
            Mpts[j] *= localToWorldMatrix.inverse();
    }
    // set positions
    itGeo.setAllPositions(Mpts);
    
    
    // set vertex colour
    if(visualisationMode != VM_OFF){
        std::vector<double> ptsColour(numPts, 0.0);
        if(visualisationMode == VM_STIFFNESS){
            for(int i=0;i<numPts;i++){
                ptsColour[i] = 1.0 - ptsWeight[i];
            }
        }else if(visualisationMode == VM_EFFECT){
            for(int j=0;j<numPts;j++){
//                ptsColour[j] = std::accumulate(wr[j].begin(), wr[j].end(), 0.0);
                ptsColour[j] = visualisationMultiplier * wr[j][numPrb-1];
            }
        }
        visualise(data, outputGeom, ptsColour);
    }

    return MS::kSuccess;
}

/// harmonic weighting
void probeDeformerNode::harmonicWeight(const std::vector<double>& probeWeight, const std::vector<int>& faceList,
                                           const std::vector<Vector3d>& pts, const std::vector< std::vector<double> >& dist){
    std::vector<Matrix4d> P;
    int numPts=(int)pts.size();
    std::vector<vertex> dummyVertexList;
    std::vector<int> tetList;
    std::vector<edge> dummyEdgeList;
    makeTetList(TM_FACE, numPts, faceList, dummyEdgeList, dummyVertexList, tetList);
    makeTetMatrix(TM_FACE, pts, tetList, faceList, dummyEdgeList, dummyVertexList, P);
    int num = (int)tetList.size()/4;
    // find closest points to probes
    std::vector<int> closestPts(numPrb);
    for(int i=0;i<numPrb;i++){
        closestPts[i] = 0;
        double min_d = HUGE_VAL;
        for(int j=0;j<numPts;j++){
            if( dist[j][i] < min_d){
                min_d = dist[j][i];
                closestPts[i] = j;
            }
        }
    }
    // LHS
    std::vector<T> tripletListMat(0);
    tripletListMat.reserve(num*16+2*numPrb);
    Matrix4d Hlist;
	Matrix4d diag=Matrix4d::Identity();
	diag(3,3)=0;
	for(int i=0;i<num;i++){
        P[i]=P[i].inverse().eval();
		Hlist=P[i].transpose()*diag*P[i];
		for(int j=0;j<4;j++){
			for(int k=0;k<4;k++){
                tripletListMat.push_back(T(tetList[4*i+j],tetList[4*i+k],Hlist(j,k)));
			}
		}
	}
    // set hard constraint
    for(int i=0;i<numPrb;i++){
        tripletListMat.push_back(T(numPts+num+i,closestPts[i],1));
        tripletListMat.push_back(T(closestPts[i],numPts+num+i,1));
    }
    SpMat H(numPts+num+numPrb, numPts+num+numPrb), G(numPts+num+numPrb,numPrb);
    H.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    // factorise
    SparseLU<SpMat> weightSolver;
    //weightSolver.isSymmetric(true);
    weightSolver.compute(H);
    if(weightSolver.info() != Success){
        //        std::string error_mes = solver.lastErrorMessage();
        MGlobal::displayInfo("Cleanup the mesh first: Mesh menu => Cleanup => Remove zero edges, faces");
    }
    // RHS
    tripletListMat.clear();
    tripletListMat.reserve(numPrb);
    for(int i=0;i<numPrb;i++){
        tripletListMat.push_back(T( numPts+num+i, i, probeWeight[i]));
    }
    G.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    // solve
    SpMat Sol = weightSolver.solve(G);
    for (int i=0;i<numPrb; i++){
        for(int j=0;j<numPts; j++){
            wr[j][i] = ws[j][i] = wl[j][i] = Sol.coeff(j,i);
        }
    }
}


// create attr
MStatus probeDeformerNode::initialize()
{
    MFnTypedAttribute tAttr;
    MFnNumericAttribute nAttr;
    MFnEnumAttribute eAttr;
    MFnMatrixAttribute mAttr;
   	MRampAttribute rAttr;

    // this attr will be dirtied when weight recomputation is needed
    aComputeWeight = nAttr.create( "computeWeight", "computeWeight", MFnNumericData::kBoolean, true );
    nAttr.setStorable(false);
    nAttr.setKeyable(false);
    addAttribute( aComputeWeight );

    aMatrix = mAttr.create("probeMatrix", "pm");
    mAttr.setStorable(false);
    mAttr.setHidden(true);
    mAttr.setArray(true);
    mAttr.setUsesArrayDataBuilder(true);
    addAttribute(aMatrix);
    attributeAffects( aMatrix, outputGeom );

    aInitMatrix = mAttr.create("initProbeMatrix", "ipm");
    mAttr.setHidden(true);
    mAttr.setArray(true);
    mAttr.setStorable(true);
    mAttr.setUsesArrayDataBuilder(true);
    addAttribute(aInitMatrix);
    attributeAffects( aInitMatrix, outputGeom );

    aBlendMode = eAttr.create( "blendMode", "bm", BM_SRL );
    eAttr.addField( "expSO+expSym", BM_SRL );
    eAttr.addField( "expSE+expSym", BM_SES );
    eAttr.addField( "logmatrix3", BM_LOG3 );
    eAttr.addField( "logmatrix4", BM_LOG4 );
    eAttr.addField( "quat+linear", BM_QSL );
    eAttr.addField( "linear", BM_AFF );
    eAttr.addField( "off", BM_OFF );
    eAttr.setStorable(true);
    addAttribute( aBlendMode );
    attributeAffects( aBlendMode, outputGeom );

	aRotationConsistency = nAttr.create( "rotationConsistency", "rc", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aRotationConsistency );
    attributeAffects( aRotationConsistency, outputGeom );

	aFrechetSum = nAttr.create( "frechetSum", "fs", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aFrechetSum );
    attributeAffects( aFrechetSum, outputGeom );

	aWorldMode = nAttr.create( "worldMode", "wrldmd", MFnNumericData::kBoolean, true );
    nAttr.setStorable(true);
    addAttribute( aWorldMode );
    attributeAffects( aWorldMode, outputGeom );
    attributeAffects( aWorldMode, aComputeWeight );
    
    aWeightMode = eAttr.create( "weightMode", "wtm", WM_INV_DISTANCE );
    eAttr.addField( "inverse", WM_INV_DISTANCE );
    eAttr.addField( "cutoff", WM_CUTOFF_DISTANCE );
    eAttr.addField( "draw", WM_DRAW );
    eAttr.addField( "harmonic", WM_HARMONIC);
    eAttr.setStorable(true);
    addAttribute( aWeightMode );
    attributeAffects( aWeightMode, outputGeom );
    attributeAffects( aWeightMode, aComputeWeight );

	aMaxDist = nAttr.create("maxDistance", "md", MFnNumericData::kDouble, 8.0);
    nAttr.setMin( 0.001 );
    nAttr.setStorable(true);
	addAttribute( aMaxDist );
	attributeAffects( aMaxDist, outputGeom );
	attributeAffects( aMaxDist, aComputeWeight );

	aNormExponent = nAttr.create("normExponent", "ne", MFnNumericData::kDouble, 1.0);
    nAttr.setStorable(true);
	addAttribute( aNormExponent );
	attributeAffects( aNormExponent, outputGeom );
	attributeAffects( aNormExponent, aComputeWeight );
    
    aProbeWeight = nAttr.create("probeWeight", "prw", MFnNumericData::kDouble, 1.0);
    nAttr.setArray(true);
    nAttr.setStorable(true);
    nAttr.setUsesArrayDataBuilder(true);
    addAttribute(aProbeWeight);
	attributeAffects( aProbeWeight, outputGeom );
	attributeAffects( aProbeWeight, aComputeWeight);

    aVisualisationMode = eAttr.create( "visualisationMode", "vm", VM_OFF );
    eAttr.addField( "off", VM_OFF );
    eAttr.addField( "effect", VM_EFFECT );
    eAttr.addField( "stiffness", VM_STIFFNESS );
    eAttr.setStorable(true);
    addAttribute( aVisualisationMode );
    attributeAffects( aVisualisationMode, outputGeom );
    
    aVisualisationMultiplier = nAttr.create("visualisationMultiplier", "vmp", MFnNumericData::kDouble, 1.0);
    nAttr.setStorable(true);
	addAttribute( aVisualisationMultiplier );
	attributeAffects( aVisualisationMultiplier, outputGeom );

	//ramp
    aWeightCurveR = rAttr.createCurveRamp( "weightCurveRotation", "wcr" );
    addAttribute( aWeightCurveR );
	attributeAffects( aWeightCurveR, outputGeom );
    aWeightCurveS = rAttr.createCurveRamp( "weightCurveShear", "wcs" );
    addAttribute( aWeightCurveS );
	attributeAffects( aWeightCurveS, outputGeom );
    aWeightCurveL = rAttr.createCurveRamp( "weightCurveTranslation", "wcl" );
    addAttribute( aWeightCurveL );
	attributeAffects( aWeightCurveL, outputGeom );

    return MS::kSuccess;
}

// create ramp attributes
MStatus probeDeformerNode::accessoryNodeSetup(MDagModifier& cmd)
{
	MStatus stat;
	MObject thisNode = thisMObject();
	MRampAttribute rWeightCurveR( thisNode, probeDeformerNode::aWeightCurveR, &stat );
	MRampAttribute rWeightCurveS( thisNode, probeDeformerNode::aWeightCurveS, &stat );
	MRampAttribute rWeightCurveL( thisNode, probeDeformerNode::aWeightCurveL, &stat );
	
	MFloatArray a1,b1;// position, value
	MIntArray c1;// interpolation

	a1.append(float(0.0));
    a1.append(float(1.0));
    
    b1.append(float(1.0));
    b1.append(float(0.0));
    
	c1.append(MRampAttribute::kSmooth);
	c1.append(MRampAttribute::kSmooth);
    
    rWeightCurveR.addEntries(a1,b1,c1);
    rWeightCurveS.addEntries(a1,b1,c1);
    rWeightCurveL.addEntries(a1,b1,c1);
    return stat;
}


// this deformer also changes colours
void probeDeformerNode::postConstructor(){
	setDeformationDetails(kDeformsColors);
}

// plugin initializer
MStatus initializePlugin( MObject obj )
{
    MStatus status;
    MFnPlugin plugin( obj, "Shizuo KAJI", "0.1", "Any");
    status = plugin.registerNode( probeDeformerNode::nodeName, probeDeformerNode::id, probeDeformerNode::creator, probeDeformerNode::initialize, MPxNode::kDeformerNode );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    return status;
}
 
MStatus uninitializePlugin( MObject obj )
{
    MStatus   status;
    MFnPlugin plugin( obj );
    status = plugin.deregisterNode( probeDeformerNode::id );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    return status;
}
